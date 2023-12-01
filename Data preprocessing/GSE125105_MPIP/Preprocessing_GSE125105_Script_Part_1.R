# Preprocessing Script for the GSE125105
# Raw data files could be obtained from corresponding repository at GSE125105 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125105
# Note 1: This script is not expected to be executed since some of raw files are not provided (Contact authors for authorized request if needed)
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


################### Importing data ###################
# Preparing target file
All_files = list.files("...")
All_files = All_files[stri_detect_fixed(All_files, pattern ="Red.idat")]
Full_file_name = stri_replace_all_fixed(All_files, pattern = "_Red.idat", replacement = "")
Participant = sapply(Full_file_name, function(x){
  x = unlist(stri_split_fixed(x, pattern = "_"))
  x = x[1]
  return(x)
})
Sentrix.Barcode = sapply(Full_file_name, function(x){
  x = unlist(stri_split_fixed(x, pattern = "_"))
  x = x[2]
  return(x)
})
Sample.Section = sapply(Full_file_name, function(x){
  x = unlist(stri_split_fixed(x, pattern = "_"))
  x = x[3]
  return(x)
})
barcode = sapply(Full_file_name, function(x){
  x = unlist(stri_split_fixed(x, pattern = "_"))
  x = x[-1]
  x = paste0(x, collapse = "_")
  return(x)
})
Targets_GSE125105 = data.frame(Full_file_name = Full_file_name, Participant = Participant, Sentrix.Barcode = Sentrix.Barcode, Sample.Section = Sample.Section, barcode = barcode)
rm(list = c("All_files", "Full_file_name","Participant", "Sentrix.Barcode", "Sample.Section", "barcode"))

# Reading IDAT files
IDAT_path = "..."
Mset.raw = read.metharray.exp(base = IDAT_path, force=T, recursive = TRUE)

# Adjusting corresponding target file
Targets_GSE125105 = Targets_GSE125105[match(Mset.raw@colData@rownames, Targets_GSE125105$Full_file_name),]
Mset.raw@colData@rownames == Targets_GSE125105$Full_file_name

# Getting annotation
annot = as.data.frame(getAnnotation(Mset.raw))
write.csv(annot, "Annot_Illumina_450_K.csv")

################### Background correction and filtering ###################
# Background correction with dye-bias normalization ("noob" method)
Mset.A = preprocessNoob(Mset.raw)

# Background correction with Funnorm (may be required)
Mset.B = preprocessFunnorm(Mset.raw)
save(Mset.B, file = "Mset_Funnorm_normalized_bumphunt_combined_GSE125105.rda")

# Saving beta values after preprocessing with "noob"
Beta_Noob_Normalized = getBeta(Mset.A, offset = 100) # Offset of 100 is Illumina's default
write.table(Beta_Noob_Normalized, "Beta_Noob_Normalized_combined_GSE125105.txt", sep = "\t")

# Obtaining beta value
betQN = betaqn(Mset.A) # from the wateRmelon package
betQN[betQN == 1] = 1 - min(betQN)

# Normalization for beta
annot$Type = ifelse(annot$Type == "I", 1,2)
betQN.BMIQ = BMIQ(betQN, as.numeric(annot[rownames(betQN), "Type"])) # from the wateRmelon package
betas = as.data.frame(betQN.BMIQ$nbeta)
pvals = detectionP(Mset.raw) # getting detection p-values

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
lines_to_select_ch = rownames(annot)[which(substr(rownames(annot), 1, 2) == 'cg')]

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
unadj.M = ComBat(M.val, batch = Targets_GSE125105$Sentrix.Barcode)
unadj.beta = m2beta(unadj.M)
write.table(unadj.beta, "unadjusted beta_combined_GSE125105.txt", sep = "\t")
write.table(unadj.M, "unadjusted Mval_combined_GSE125105.txt", sep = "\t")

################### Cell-type adjustment ###################
# cell blood estimates
wbcset = estimateCellCounts(Mset.raw, returnAll = TRUE)
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
M_adj = beta2m(adj.betas)
Targets_GSE125105$Full_file_name == colnames(M_adj)

# Writing output
write.table(adj.betas,file="adjusted_blood_beta_combined_GSE125105.txt", sep = "\t")
write.table(M_adj,"M_values_adj_combined_GSE125105.txt", sep = "\t")
write.csv(Targets_GSE125105, "Targets_GSE125105.csv", row.names = FALSE)