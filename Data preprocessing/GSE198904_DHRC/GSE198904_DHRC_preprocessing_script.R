Working_directory = "/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE198904_DHRC" # Replace with an appropriate directory
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
gse = getGEO("GSE198904")
gse = gse[[1]]
GSE198904_phenotypes = gse@phenoData@data
colnames(GSE198904_phenotypes) = stri_replace_all_fixed(colnames(GSE198904_phenotypes), pattern = ":ch1", replacement = "")
GSE198904_phenotypes_cohort_2 = GSE198904_phenotypes[GSE198904_phenotypes$cohort == "Cohort 2",]
Duplicated_subjects = as.data.frame(table(GSE198904_phenotypes_cohort_2$subject))
# Many subjects are duplicated -> this should be addressed

# Supplementary files for GSE198904 were downloaded manually and placed in /home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE198904_supplementary
# Reading files
GSE198904_DHRC_signals = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE198904_supplementary/GSE198904_DHRC_Matrix_signal_intensities.txt")
GSE198904_DHRC_p = GSE198904_DHRC_signals[,seq(from = 1, to = ncol(GSE198904_DHRC_signals)-2, by=3)]
fwrite(GSE198904_DHRC_signals, "GSE198904_DHRC_signals.csv", row.names = TRUE, sep = ";")
rm(GSE198904_DHRC_signals)
gc()
# "blood idoloptimized epic" will be used for cell counts

# Raw MethylSet
GSE198904_DHRC_raw = readGEORawFile(filename = "GSE198904_DHRC_signals.csv",
                               Uname = "Unmethylated_Signal",
                               Mname = "Methylated_Signal",
                               array = "IlluminaHumanMethylationEPIC",
                               annotation = "ilm10b2.hg19",
                               sep = ";")

Beta_from_signals_raw = getBeta(GSE198904_DHRC_raw, offset = 100)
fwrite(as.data.frame(Beta_from_signals_raw), "Beta_from_signals_raw_GSE198904_DHRC.csv", row.names = TRUE, sep = ",")
rm(GSE198904_DHRC_raw)
gc()

# Performing quantile normalization as in other cohorts
# Obtaining beta value
betQN = betaqn(Beta_from_signals_raw) # from the wateRmelon package
betQN[betQN == 1] = 1 - min(betQN)

# Getting annotation and matching the rows
annot = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
annot = as.data.frame(annot)
annot$Type = ifelse(annot$Type == "I",1,2)
annot = annot[rownames(betQN),]
all(rownames(betQN) == rownames(annot)) # rows are fully matching
gc()

# Beta-Mixture Quantile (BMIQ) Normalization to correct for probe type bias
betQN.BMIQ = BMIQ(betQN, annot$Type)
betas = as.data.frame(betQN.BMIQ$nbeta)
rm(betQN)
gc()

# Need to obtain detection P-vals
pvals = GSE198904_DHRC_p
pvals = pvals[rownames(betas),]

# Check
all(colnames(betas) == colnames(pvals)) # TRUE
all(rownames(betas) == rownames(pvals)) # TRUE
rm(GSE198904_DHRC_p)
gc()

################### Filtering ###################

# First of all, filter the samples (we have to keep only the samples that are of good quality)
cols_to_select_samples=colSums(pvals<=0.00005)>=0.75*nrow(pvals) #Select samples with more than 75% of probes with detection P-value less than 0.00005

# Then, filter the probes (we have to keep only the probes that are of good quality)
lines_to_select_pvals=rownames(pvals)[rowSums(pvals<=0.01)>=0.75*ncol(pvals)] #Select probes with more than 75% of samples with detection P-value less than 0.01

# Remove sites on X and Y chromosomes
indexes_to_select_X=which(annot[,'chr']!='chrX')
indexes_to_select_Y=which(annot[,'chr']!='chrY')
indexes_to_select_XY=intersect(indexes_to_select_X,indexes_to_select_Y)
lines_to_select_XY=rownames(annot[indexes_to_select_XY,])

# Remove probes that do not target CpGs (starting with "ch")
lines_to_select_ch=rownames(annot)[which(substr(rownames(annot),1,2)=='cg')]

# Remove sites with SNPs within probe with maf >5% 
indexes_to_select_SNP=which(is.na(annot[,'Probe_maf'])==T | annot[,'Probe_maf']<0.05)
lines_to_select_SNP=rownames(annot[indexes_to_select_SNP,])

# Remove sites with SNPs within the SBE and CpG
indexes_to_select_SNP_SBE_CPG = which(is.na(annot[,'CpG_rs'])==T & is.na(annot[,'SBE_rs'])==T )
lines_to_select_SNP_SBE_CPG =rownames(annot[indexes_to_select_SNP_SBE_CPG,])
#
test = annot[indexes_to_select_SNP_SBE_CPG,]
table(test$CpG_rs) # questionable CpGs have been removed
table(test$SBE_rs) # questionable probes with SBE SNP have been removed
rm(test)


# Intersect all these probes that were selected, to only obtain the probes that pass all these filterings
lines_to_select=intersect(lines_to_select_XY,lines_to_select_pvals)
lines_to_select=intersect(lines_to_select,lines_to_select_SNP)
lines_to_select=intersect(lines_to_select,lines_to_select_ch)
lines_to_select = intersect(lines_to_select,lines_to_select_SNP_SBE_CPG)

# Get the beta-value matrix without probes that didn't pass the quality control and samples that didn't
betas2=betas[lines_to_select,cols_to_select_samples]

# Additional filtering for Illumina 450K
Cross_reactive_probes_YI_AN_CHEN = read.csv("/home/aleksandr/Desktop/WORK/Bad probes illumina 450K/48639-non-specific-probes-Illumina450k.csv") #https://www.tandfonline.com/doi/full/10.4161/epi.23470; #https://github.com/sirselim/illumina450k_filtering
betas2 = betas2[rownames(betas2) %!in% Cross_reactive_probes_YI_AN_CHEN$TargetID,]
Cross_reactive_probes_MILES_BENTON = as.data.frame(fread("/home/aleksandr/Desktop/WORK/Bad probes illumina 450K/HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt", header = FALSE)) #https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0569-x; #https://github.com/sirselim/illumina450k_filtering
betas2 = betas2[rownames(betas2) %!in% Cross_reactive_probes_MILES_BENTON$V1,]

# Removing EPIC bad probes
Files_EPIC_bad_probes = list.files("/home/aleksandr/Desktop/WORK/Files with pverlapping probes EPIC") #https://github.com/sirselim/illumina450k_filtering; #https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1#Sec22
Files_EPIC_bad_probes = Files_EPIC_bad_probes[stri_detect_fixed(Files_EPIC_bad_probes, pattern = ".csv")]
Files_EPIC_bad_probes = paste0("/home/aleksandr/Desktop/WORK/Files with pverlapping probes EPIC/", Files_EPIC_bad_probes)
Files_EPIC_bad_probes = lapply(Files_EPIC_bad_probes, function(x) as.data.frame(fread(x)))
Bad_probes_Epic = lapply(Files_EPIC_bad_probes, function(x) x[,1])
Bad_probes_Epic = unlist(Bad_probes_Epic)
Bad_probes_Epic = unique(Bad_probes_Epic)
betas2 = betas2[rownames(betas2) %!in% Bad_probes_Epic,]

# Remove sites with missing beta values
betas2=na.omit(betas2)
# Transform beta to M values
M.val=beta2m(betas2)

rm(betas)
rm(Files_EPIC_bad_probes)
rm(Cross_reactive_probes_YI_AN_CHEN)
rm(Cross_reactive_probes_MILES_BENTON)
rm(Bad_probes_Epic)
rm(betQN.BMIQ)
rm(pvals)
gc()

################### Phenotype curation ###################

# Check
all(GSE198904_phenotypes_cohort_2$title == colnames(betas2)) # Not matching
all(GSE198904_phenotypes_cohort_2$title %in% colnames(betas2)) # All present
nrow(GSE198904_phenotypes_cohort_2) # 432
ncol(betas2) # 432

# Match the order
betas2 = betas2[,GSE198904_phenotypes_cohort_2$title]
M.val = M.val[,GSE198904_phenotypes_cohort_2$title]
all(GSE198904_phenotypes_cohort_2$title == colnames(betas2)) # TRUE
all(GSE198904_phenotypes_cohort_2$title == colnames(M.val)) # TRUE

# Obtain batch variable
GSE198904_phenotypes_cohort_2$Batch = sapply(GSE198904_phenotypes_cohort_2$title, function(x){
  x = stri_split_fixed(x, pattern = "_")
  x = unlist(x)
  x = x[1]
  return(x)
})
table(GSE198904_phenotypes_cohort_2$Batch)
gc()

################### Save files before corrections ###################
fwrite(GSE198904_phenotypes_cohort_2, "GSE198904_phenotypes_cohort_2.csv", row.names = TRUE, sep = ";")
gc()
fwrite(betas2, "betas2.csv", row.names = TRUE, sep = ";")
gc()
fwrite(M.val, "Mval.csv", row.names = TRUE, sep = ";")
gc()
