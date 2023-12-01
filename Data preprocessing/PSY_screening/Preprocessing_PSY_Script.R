# Preprocessing Script for the PSY Cohort
# Raw data files are not provided (Ethical permission limitation)
# Raw data files could be provided based on authorized reasonable request to a corresponding author 
#  -> if approved by an ethical review board of Uppsala University
# Note 1: This script is not expected to be executed since raw files and paths are not included
# Note 2: Computer-specific file paths are shown as "..."
# Note 3: The equal sign = was used as an assignment operator for typing/productivity reasons
# Note 4: In many cases loops were deliberately used instead of apply functions to enable better control of the variables

Working_directory = "..." # Replace with an appropriate directory
setwd(Working_directory)

# Importing required packages
library(limma)
library(FlowSorted.Blood.450k)
library(minfi)
library(wateRmelon)
library(lumi)
library(lmtest)
library(sva)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylation450kmanifest)
library(stringi)
library(stringr)
library(data.table)

# Not IN operator
'%!in%' = function(x,y){!('%in%'(x,y))}

# The data was available from 2 independent batches

################### Importing data ###################
# Reading targets files
Targets_SCR_1 = read.table("...", sep="\t", header=TRUE)
Targets_SCR_2 = read.table("...", sep="\t", header=TRUE)

# Adjusting targets files
# 1
colnames(Targets_SCR_1)[1] = "Basename"
Targets_SCR_1$Basename = as.character(Targets_SCR_1$Basename)
Targets_SCR_1$Sentrix.Barcode = as.character(Targets_SCR_1$Sentrix.Barcode)
Targets_SCR_1$Sample.Section = as.character(Targets_SCR_1$Sample.Section)
Targets_SCR_1$barcode = paste0(Targets_SCR_1$Sentrix.Barcode,"_", Targets_SCR_1$Sample.Section)
# 2
colnames(Targets_SCR_2)[1]="Basename"
Targets_SCR_2$Basename=as.character(Targets_SCR_2$Basename)
Targets_SCR_2$Sentrix.Barcode=as.character(Targets_SCR_2$Sentrix.Barcode)
Targets_SCR_2$Sample.Section=as.character(Targets_SCR_2$Sample.Section)
Targets_SCR_2$barcode = paste0(Targets_SCR_2$Sentrix.Barcode,"_", Targets_SCR_2$Sample.Section)

# Preparing IDAT paths 
IDAT_path_SCR1 = "/home/aleksandr/Desktop/WORK/Preprocessing_SCR_Combinded/FTOPSY450K/FTOPSY450KIDAT3/MF-0305_131107_IDAT/MF-0305_131107_IDAT"
IDAT_path_SCR2 = "/home/aleksandr/Desktop/WORK/Preprocessing_SCR_Combinded/OB-0737_150326_ResultReport/OB-0737_150326_IDAT"

# Getting RG sets
Mset.raw_SCR1 = read.metharray.exp(base = IDAT_path_SCR1, force = T, recursive = TRUE)
Mset.raw_SCR2 = read.metharray.exp(base = IDAT_path_SCR2, force = T, recursive = TRUE)

# Adjusting corresponding target files
Targets_SCR_1 = Targets_SCR_1[match(Mset.raw_SCR1@colData@rownames, Targets_SCR_1$barcode),]
Targets_SCR_2 = Targets_SCR_2[match(Mset.raw_SCR2@colData@rownames, Targets_SCR_2$barcode),]

# Making combined target file 
Targets_SCR_combined = rbind(Targets_SCR_1, Targets_SCR_2)
Targets_SCR_combined$Study = sapply(Targets_SCR_combined$barcode, function(x){
  if (x %in% Targets_SCR_1$barcode){return("SCR1")} else {return("SCR2")}
})

# Making combined RG set
Mset.raw_SCR_combined = combineArrays(Mset.raw_SCR1, Mset.raw_SCR2, outType = c("IlluminaHumanMethylation450k"))

# Check for order of targets
Mset.raw_SCR_combined@colData@rownames == Targets_SCR_combined$barcode # order is matching

# Getting annotation
annot = as.data.frame(getAnnotation(Mset.raw_SCR_combined))
write.csv(annot, "Annot_Illumina_450_K.csv")


################### Background correction and filtering ###################
# Background correction with dye-bias normalization ("noob" method)
Mset.A = preprocessNoob(Mset.raw_SCR_combined)

# Background correction with Funnorm (may be required)
Mset.B = preprocessFunnorm(Mset.raw_SCR_combined)
save(Mset.B, file = "Mset_Funnorm_normalized_bumphunt_combined_SCR_PSY.rda")

# Saving beta values after preprocessing with "noob"
Beta_Noob_Normalized = getBeta(Mset.A, offset = 100) # Offset of 100 is Illumina's default
write.table(Beta_Noob_Normalized,"Beta_Noob_Normalized_combined_SCR_PSY.txt", sep = "\t")

# Obtaining beta value
betQN = betaqn(Mset.A) # from the wateRmelon package
betQN[betQN == 1] = 1 - min(betQN)

# Normalization for beta
annot$Type = ifelse(annot$Type == "I", 1,2)
betQN.BMIQ = BMIQ(betQN, as.numeric(annot[rownames(betQN), "Type"])) # from the wateRmelon package
betas = as.data.frame(betQN.BMIQ$nbeta)
pvals = detectionP(Mset.raw_SCR_combined) # getting detection p-values

# First of all, filter the samples (we have to keep only the samples that are of good quality)
cols_to_select_samples = colSums(pvals <= 0.00005) >= 0.75 * nrow(pvals) #Select samples with more than 75% of probes with detection P-value less than 0.00005

# Then, filter the probes (we have to keep only the probes that are of good quality)
lines_to_select_pvals = rownames(pvals)[rowSums(pvals <= 0.01) >= 0.75 * ncol(pvals)] #Select probes with more than 75% of samples with detection P-value less than 0.01

# Remove sites on X and Y chromosomes
indexes_to_select_X = which(annot[,'chr'] != 'chrX')
indexes_to_select_Y = which(annot[,'chr'] != 'chrY')
indexes_to_select_XY = intersect(indexes_to_select_X,indexes_to_select_Y)
lines_to_select_XY = rownames(annot[indexes_to_select_XY,])

# Remove probes that do not target CpGs (starting with "ch")
lines_to_select_ch = rownames(annot)[which(substr(rownames(annot),1,2) == 'cg')]

# Remove sites with SNPs within probe with maf >5% 
indexes_to_select_SNP = which(is.na(annot[,'Probe_maf']) == T | annot[,'Probe_maf'] < 0.05)
lines_to_select_SNP = rownames(annot[indexes_to_select_SNP,])

# Remove sites with SNPs within the SBE and CpG
indexes_to_select_SNP_SBE_CPG = which(is.na(annot[,'CpG_rs']) == T & is.na(annot[,'SBE_rs']) == T)
lines_to_select_SNP_SBE_CPG = rownames(annot[indexes_to_select_SNP_SBE_CPG,])

# Inspecting removed probes
test = annot[indexes_to_select_SNP_SBE_CPG,]
table(test$CpG_rs) # questionable CpGs have been removed
table(test$SBE_rs) # questionable probes with SBE SNP have been removed
rm(test)

# Intersect all these probes that were selected, to only obtain the probes that pass all these filtering steps
lines_to_select = intersect(lines_to_select_XY,lines_to_select_pvals)
lines_to_select = intersect(lines_to_select,lines_to_select_SNP)
lines_to_select = intersect(lines_to_select,lines_to_select_ch)
lines_to_select = intersect(lines_to_select,lines_to_select_SNP_SBE_CPG)

# Get the beta-value matrix without probes that didn't pass the quality control and samples that didn't
betas2 = betas[lines_to_select,cols_to_select_samples]

# Additional filtering for Illumina 450K
Cross_reactive_probes_YI_AN_CHEN = read.csv("48639-non-specific-probes-Illumina450k.csv") #https://www.tandfonline.com/doi/full/10.4161/epi.23470; #https://github.com/sirselim/illumina450k_filtering
betas2 = betas2[rownames(betas2) %!in% Cross_reactive_probes_YI_AN_CHEN$TargetID,]
Cross_reactive_probes_MILES_BENTON = as.data.frame(fread("HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt", header = FALSE)) #https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0569-x; #https://github.com/sirselim/illumina450k_filtering
betas2 = betas2[rownames(betas2) %!in% Cross_reactive_probes_MILES_BENTON$V1,]

# Removing sites with missing beta values
betas2 = na.omit(betas2)

# Transform beta to M values
M.val = beta2m(betas2)

# Batch effect
unadj.M = ComBat(M.val, batch = Targets_SCR_combined$Sentrix.Barcode)
unadj.beta = m2beta(unadj.M)
write.table(unadj.beta, "unadjusted beta_combined_SCR_PSY.txt", sep= "\t")

################### Cell-type adjustment ###################
# cell blood estimates
wbcset = estimateCellCounts(Mset.raw_SCR_combined, returnAll = TRUE)
setbun = as.data.frame(wbcset$counts)
write.table(setbun, file = "wbcset.txt", sep = "\t")

# adjusting beta values
beta_matrix=as.data.frame(unadj.beta)

y = matrix(nrow = nrow(beta_matrix), ncol = ncol(beta_matrix))
rownames(y) = rownames(beta_matrix)
colnames(y) = colnames(beta_matrix)

i=0
while (i < nrow(beta_matrix))  
{
  i = i + 1
  CpG = as.numeric(as.character(beta_matrix[i,]))
  beta.lm = lm(CpG~CD8T+CD4T+NK+Bcell+Mono+Gran, data = setbun)
  adj.beta = mean(CpG) + as.numeric(as.character(beta.lm$residuals))
  y[i,] = adj.beta
}

x=y

# Final adjustment of extreme values
y[y>1] = NaN
y[y<0] = NaN
min.y = min(y,na.rm = T)
max.y = max(y,na.rm = T)
x[x>1] = max.y # Should be max.y
x[x<0] = min.y # Should be min.y
min.x = min(x)

adj.betas = x
adj.betas = as.data.frame(adj.betas)
write.table(adj.betas, file="adjusted_blood_beta_combined_SCR_PSY.txt", sep = "\t")

M_adj = beta2m(adj.betas)
write.table(M_adj, "M_values_adj_combined_SCR_PSY.txt", sep = "\t")

# Selecting required M-values (runs contained samples from several cohorts)
Illum_target = Targets_SCR_combined
Illum_target = Illum_target[stri_detect_fixed(Illum_target$Basename, pattern = "PSY"),]
M_adj_SCR_combined = M_adj[,Illum_target$barcode]
adj.betas_SCR_combined = adj.betas[,Illum_target$barcode]
Illum_target$Participant = sapply(Illum_target$Basename, function(x){
  if (stri_detect_fixed(x, pattern = ".")){
    x = unlist(stri_split_fixed(x, pattern = "."))
  }
  x = unlist(stri_split_fixed(x, pattern = "_"))
  x = x[2]
  x = str_trim(x)
  return(x)
})
New_col_names = sapply(colnames(M_adj_SCR_combined), function(x){
  index = which(Illum_target$barcode == x)
  index = index[1]
  ID = Illum_target$Participant[index]
  return(ID)
})
which(duplicated(New_col_names)) # Some samples were analyzed in both batches

# Selecting only first instances of duplicated samples
Columns_to_select = !duplicated(New_col_names)
M_adj_SCR_combined = M_adj_SCR_combined[,Columns_to_select]
adj.betas_SCR_combined = adj.betas_SCR_combined[,Columns_to_select]

# Writing output
write.table(M_adj_SCR_combined,"M_values_adj_221_PSY_screen_comb_prepr.txt", sep = "\t")
write.table(adj.betas_SCR_combined,"Beta_values_adj_221_PSY_screen_comb_prepr.txt", sep = "\t")
write.csv(Illum_target, "SAMPLE_ID_FOR_PSY_screen_comb_prepr.csv", row.names = FALSE)








