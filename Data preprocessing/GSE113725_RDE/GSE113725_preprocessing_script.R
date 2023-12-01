Working_directory = "/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE113725_RDE" # Replace with an appropriate directory
setwd(Working_directory)

# Setting options
getOption("scipen") # default number notation is 0
options(scipen=999)
options(stringsAsFactors = FALSE)

# Importing packages (not all of them are required)
library(stringr)
library(dplyr)
library(openxlsx)
library(grid)
library(gdata)
library(magrittr)
library(data.table)
library(XML)
library(RCurl)
library(stringi)
library(httr)
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
library(GEOquery)
library(meffil)

################### Defining functions ###################

# NOT IN operator
'%!in%' = function(x,y){!('%in%'(x,y))}

# A function to quickly read data
smart_fread = function(x){
  x = as.data.frame(fread(x, nThread = 14))
  rownames(x) = x$V1
  x$V1 = NULL
  return(x)
}

################### Importing data ###################
# Obtaining data from GEO
gse = getGEO("GSE113725")
gse = gse[[1]]
GSE113725_phenotypes = gse@phenoData@data

# Getting suppl. files for GSE113725
getGEOSuppFiles("GSE113725")
files = list.files("GSE113725")
files = files[stri_detect_fixed(files, pattern = ".gz")]
files = paste0("GSE113725", "/", files)
lapply(files, gunzip)

# Reading files
GSE113725_methyl = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE113725_RDE/GSE113725/GSE113725_methylatedIntensities.csv")
GSE113725_non_methyl = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE113725_RDE/GSE113725/GSE113725_unmethylatedIntensities.csv")
GSE113725_det_p = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE113725_RDE/GSE113725/GSE113725_detectionP.csv")

all(colnames(GSE113725_methyl) == colnames(GSE113725_non_methyl)) # TRUE
all(colnames(GSE113725_methyl) == colnames(GSE113725_det_p)) # TRUE
all(rownames(GSE113725_methyl) == rownames(GSE113725_non_methyl)) # TRUE
all(rownames(GSE113725_methyl) == rownames(GSE113725_det_p)) # TRUE

colnames(GSE113725_methyl) = paste0(colnames(GSE113725_methyl), " Methylated signal")
colnames(GSE113725_non_methyl) = paste0(colnames(GSE113725_non_methyl), " Unmethylated signal")
merged_signals = cbind(GSE113725_non_methyl, GSE113725_methyl)
fwrite(merged_signals, "merged_signals.csv", row.names = TRUE, sep = ",")

# Making RAW GenomicMethylSet 
GSE113725_raw = readGEORawFile(filename = "merged_signals.csv",
                               Uname = "Unmethylated signal",
                               Mname = "Methylated signal",
                               sep = ",")

Beta_from_signals_raw = getBeta(GSE113725_raw, offset = 100)
fwrite(Beta_from_signals_raw, "Beta_from_signals_raw_GSE113725.csv", row.names = TRUE, sep = ",")


################### Normalization ###################

# Performing quantile normalization as in other cohorts
# Obtaining beta value
betQN = betaqn(Beta_from_signals_raw) # from the wateRmelon package
betQN[betQN == 1] = 1 - min(betQN)

# Getting annotation and matching the rows
annot = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot = as.data.frame(annot)
annot$Type = ifelse(annot$Type == "I",1,2)
annot = annot[rownames(betQN),]
all(rownames(betQN) == rownames(annot)) # rows are fully matching

# Beta-Mixture Quantile (BMIQ) Normalization to correct for probe type bias
betQN.BMIQ = BMIQ(betQN, annot$Type)
betas = as.data.frame(betQN.BMIQ$nbeta)

# Need to obtain detection P-vals
pvals = GSE113725_det_p
pvals = pvals[rownames(betas),]

# Check
all(colnames(betas) == colnames(pvals)) # TRUE
all(rownames(betas) == rownames(pvals)) # TRUE

################### Filtering ###################

# First of all, filter the samples (we have to keep only the samples that are of good quality)
cols_to_select_samples = colSums(pvals <= 0.00005) >= 0.75 * nrow(pvals) # Select samples with more than 75% of probes with detection P-value less than 0.00005

# Then, filter the probes (we have to keep only the probes that are of good quality)
lines_to_select_pvals = rownames(pvals)[rowSums(pvals <= 0.01) >= 0.75 * ncol(pvals)] # Select probes with more than 75% of samples with detection P-value less than 0.01

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
Cross_reactive_probes_YI_AN_CHEN = read.csv("/home/aleksandr/Desktop/WORK/Bad probes illumina 450K/48639-non-specific-probes-Illumina450k.csv") #https://www.tandfonline.com/doi/full/10.4161/epi.23470; #https://github.com/sirselim/illumina450k_filtering
betas2 = betas2[rownames(betas2) %!in% Cross_reactive_probes_YI_AN_CHEN$TargetID,]
Cross_reactive_probes_MILES_BENTON = as.data.frame(fread("/home/aleksandr/Desktop/WORK/Bad probes illumina 450K/HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt", header = FALSE)) #https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0569-x; #https://github.com/sirselim/illumina450k_filtering
betas2 = betas2[rownames(betas2) %!in% Cross_reactive_probes_MILES_BENTON$V1,]

# Removing sites with missing beta values
betas2 = na.omit(betas2)

# Transform beta to M values
M.val = beta2m(betas2) # 371308

################### Phenotype curation ###################

# Check
all(GSE113725_phenotypes$title %in% colnames(betas2)) # TRUE
all(GSE113725_phenotypes$title == colnames(betas2)) # TRUE
GSE113725_phenotypes$Batch = sapply(GSE113725_phenotypes$title, function(x){
  x = stri_split_fixed(x, pattern = "_")
  x = unlist(x)
  x = x[1]
  return(x)
})
table(GSE113725_phenotypes$Batch)
colnames(GSE113725_phenotypes) = stri_replace_all_fixed(colnames(GSE113725_phenotypes), pattern = ":ch1", replacement = "")

# Batch effect
unadj.beta = ComBat(betas2, batch = GSE113725_phenotypes$Batch)
unadj.M = beta2m(unadj.beta) # Produced NAs

# NA stats
NA_in_betas = apply(unadj.M, 1, function(x) any(is.na(x)))
table(NA_in_betas) # 144 CpGs have missing values -> removed
Clean_values = !NA_in_betas

# Selecting clean CpGs
unadj.beta = unadj.beta[Clean_values,]
unadj.M = unadj.M[Clean_values,]

write.table(unadj.beta, "unadjusted_beta_GSE113725.txt", sep = "\t")
write.table(unadj.M, "unadjusted_Mval_GSE113725.txt", sep = "\t")

# Cell proportions, EXPERIMENTAL
Cell_counts_blood = meffil.estimate.cell.counts.from.betas(Beta_from_signals_raw, cell.type.reference = "blood gse35069")
all(rownames(Cell_counts_blood) == GSE113725_phenotypes$title) # Everything is matching
setbun = as.data.frame(Cell_counts_blood)
write.table(setbun, file = "wbcset.txt", sep = "\t")

# Adjusting beta values
beta_matrix=as.data.frame(unadj.beta)

y = apply(beta_matrix, 1, function(x){
  CpG = as.numeric(x)
  beta.lm = lm(CpG~CD8T+CD4T+NK+Bcell+Mono+Gran, data = setbun)
  adj.beta = mean(CpG) + as.numeric(as.character(beta.lm$residuals))
  return(adj.beta)
})
y = t(y)
all(rownames(y) == rownames(beta_matrix))
colnames(y) = colnames(beta_matrix)
str(y)
gc()
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
all(GSE113725_phenotypes$title == colnames(M_adj)) # TRUE


# Writing output
fwrite(adj.betas,file="adjusted_blood_beta_GSE113725_full.txt", sep = "\t", row.names = TRUE) # 371164
fwrite(M_adj,"M_values_adj_GSE113725_full.txt", sep = "\t", row.names = TRUE) # 371164
write.csv(GSE113725_phenotypes, "Phenos_GSE113725_full.csv", row.names = FALSE)


################### Selecting values ###################

GSE113725_phenotypes_selected = GSE113725_phenotypes[GSE113725_phenotypes$groupid %in% c("2", "4"),]
GSE113725_phenotypes_selected$Disease = ifelse(GSE113725_phenotypes_selected$groupid == "2", "Case", "Control")
M_adj_selected = M_adj[,GSE113725_phenotypes_selected$title]

fwrite(M_adj_selected,"M_values_adj_GSE113725_selected.txt", sep = "\t", row.names = TRUE)
write.csv(GSE113725_phenotypes_selected, "Phenos_GSE113725_selected.csv", row.names = FALSE)


