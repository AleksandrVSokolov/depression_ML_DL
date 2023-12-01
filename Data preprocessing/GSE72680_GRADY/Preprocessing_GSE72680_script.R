# Preprocessing Script for the GSE72680
# Raw data files could be obtained from corresponding repository at GSE72680 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72680
# Note 1: To execute the script several raw files used have to be downloaded from GEO
# Note 2: Computer-specific file paths are shown as "..."
# Note 3: The equal sign = was used as an assignment operator as authors don't buy the idea of using <- for typing/productivity reasons
# Note 4: In many cases loops were deliberately used instead of apply functions to enable better control of the variables (even though loops in R are slow and computationally inefficient)

Working_directory = "/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE72680_GRADY"# Replace with an appropriate directory
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
gse = getGEO("GSE72680")
gse = gse[[1]]
experiment = gse@experimentData
GSE72680_phenotypes = gse@phenoData@data
# Phenotype data appeared to be incorrect since columns were skewed, continuing only with Beta values

# Working only with BETA values
Raw_betas = smart_fread("/home/aleksandr/Desktop/WORK/open_access_cohorts/E-GEOD-72680-Methyl-Mental-Trauma/GSE_version/GSE72680_beta_values.txt") # Downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72680
pvals = Raw_betas[,seq(from = 2, to = ncol(Raw_betas), by = 2)]
Raw_betas = Raw_betas[,seq(from = 1, to = ncol(Raw_betas)-1, by = 2)]
rm(Raw_betas)
rm(pvals)
# Dimensions may not be correct. Raw RGSet usually has more values. The resulting output was = 380475 CpGs

# Importing betas from signals
RAW_geo_Signals = readGEORawFile(filename = "/home/aleksandr/Desktop/WORK/open_access_cohorts/E-GEOD-72680-Methyl-Mental-Trauma/GSE_version/GSE72680_signals.txt", sep = "\t")
# Default Uname = "Unmethylated signal" and Mname = "Methylated signal" work

# Saving beta values after signals
Beta_from_signals_raw = getBeta(RAW_geo_Signals, offset = 100) # Offset of 100 is Illumina's default
write.table(Beta_from_signals_raw, "Beta_from_signals_raw_GSE72680.txt", sep = "\t")

# Performing quantile normalization as in other cohorts
# Obtaining beta value
betQN = betaqn(Beta_from_signals_raw) # from the wateRmelon package
betQN[betQN == 1] = 1 - min(betQN)

# Getting annotation and matching the rows
annot = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot = as.data.frame(annot)
annot$Type = ifelse(annot$Type == "I",1,2)
annot = annot[rownames(betQN),]
all(rownames(betQN) == rownames(annot)) # rows are not fully matching
betQN = betQN[rownames(betQN) %in% rownames(annot),]
annot = annot[rownames(betQN),]
all(rownames(betQN) == rownames(annot)) # now all rows are matching

# Beta-Mixture Quantile (BMIQ) Normalization to correct for probe type bias
betQN.BMIQ = BMIQ(betQN, annot$Type)
betas = as.data.frame(betQN.BMIQ$nbeta)

# Need to obtain detection P-vals
RAW_geo_Signals_data = smart_fread("/home/aleksandr/Desktop/WORK/open_access_cohorts/E-GEOD-72680-Methyl-Mental-Trauma/GSE_version/GSE72680_signals.txt")
pvals = RAW_geo_Signals_data[, seq(from=1, to=ncol(RAW_geo_Signals_data), by = 3)]
pvals = pvals[rownames(betas),]

# Renaming columns for pvals
colnames(pvals) = stri_replace_all_fixed(colnames(pvals), pattern = "_Detection PVal", replacement = "")

# Check
all(rownames(pvals) == rownames(betas)) # TRUE
all(colnames(pvals) == colnames(betas)) # TRUE
all(rownames(betas) == rownames(annot)) # TRUE
################### Filtering probes ###################

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
M.val = beta2m(betas2) # 380726

################### Getting phenotypes ###################

# Reading data
E_GEOD_72680_Phenotypes = as.data.frame(fread("/home/aleksandr/Desktop/WORK/open_access_cohorts/E-GEOD-72680-Methyl-Mental-Trauma/GSE_version/data", skip = 30, nThread = 14)) # Replace with an appropriate path
# the data file was downloaded from GEO (series matrix file)
# structure is complex and should be reformatted

# Preparing scaffold data dataframe with phenotypes
E_GEOD_72680_Phenotypes_Scaffold = E_GEOD_72680_Phenotypes[-(8:36),]
E_GEOD_72680_Phenotypes_Scaffold = t(E_GEOD_72680_Phenotypes_Scaffold)
E_GEOD_72680_Phenotypes_GEO_ID = rownames(E_GEOD_72680_Phenotypes_Scaffold)[-1]
colnames(E_GEOD_72680_Phenotypes_Scaffold) = E_GEOD_72680_Phenotypes_Scaffold[1,]
E_GEOD_72680_Phenotypes_Scaffold = E_GEOD_72680_Phenotypes_Scaffold[-1,]
E_GEOD_72680_Phenotypes_Scaffold = cbind(E_GEOD_72680_Phenotypes_GEO_ID, E_GEOD_72680_Phenotypes_Scaffold)
E_GEOD_72680_Phenotypes_Scaffold = as.data.frame(E_GEOD_72680_Phenotypes_Scaffold)

# Getting characteristics
Characteristics_matrix = E_GEOD_72680_Phenotypes[8:36, ]
Characteristics_matrix = apply(Characteristics_matrix[,-1], 2, function(x){
  All_traits = x
  All_traits = All_traits[!is.na(All_traits)]
  All_traits = All_traits[All_traits != ""]
  All_traits = sapply(All_traits, function(z){
    z = unlist(stri_split_fixed(z, pattern = ": "))
    z = z[1]
    return(z)})
  return(All_traits)
})
Characteristics_matrix = unlist(Characteristics_matrix)
names(Characteristics_matrix) = NULL
Characteristics_matrix = unique(Characteristics_matrix)
TMP_DF_E_GEOD_72680 = lapply(Characteristics_matrix, function(x){
  x = rep(NA, times = 422)
  return(x)
}) 
names(TMP_DF_E_GEOD_72680) = Characteristics_matrix
TMP_DF_E_GEOD_72680 = do.call(cbind, TMP_DF_E_GEOD_72680)
TMP_DF_E_GEOD_72680 = as.data.frame(TMP_DF_E_GEOD_72680)
colnames(TMP_DF_E_GEOD_72680) = stri_replace_all_fixed(Characteristics_matrix, pattern = " ", replacement = "_")
TMP_DF_E_GEOD_72680 = cbind(E_GEOD_72680_Phenotypes_GEO_ID, TMP_DF_E_GEOD_72680)

# Looping through participants and obtaining associated traits
i = 1 
while(i <= nrow(TMP_DF_E_GEOD_72680)){
  Current_participant = TMP_DF_E_GEOD_72680$E_GEOD_72680_Phenotypes_GEO_ID[i]
  Subset_data = E_GEOD_72680_Phenotypes[8:36, Current_participant]
  for (trait in Characteristics_matrix){
    if (any(stri_detect_fixed(str = Subset_data, pattern = trait))){
      Current_char = Subset_data[stri_detect_fixed(str = Subset_data, pattern = trait)]
      Value_to_return = unlist(stri_split_fixed(Current_char, pattern = ": "))
      Value_to_return = Value_to_return[2]
      TMP_DF_E_GEOD_72680[i,stri_replace_all_fixed(trait, pattern = " ", replacement = "_")] = Value_to_return
    }
  }
  rm(list = c("Current_participant", "Subset_data", "trait", "Current_char", "Value_to_return"))
  print(i)
  i = i + 1
}
E_GEOD_72680_Phenotypes = cbind(TMP_DF_E_GEOD_72680, E_GEOD_72680_Phenotypes_Scaffold[,-1])

# Removing temporary files
rm(list = c("E_GEOD_72680_Phenotypes_Scaffold", "TMP_DF_E_GEOD_72680", 
            "Characteristics_matrix", "E_GEOD_72680_Phenotypes_GEO_ID"))

# Matching order of M-values, Beta values and phenotypes
M.val = M.val[, E_GEOD_72680_Phenotypes$`!Sample_description.1`]
betas2 = betas2[, E_GEOD_72680_Phenotypes$`!Sample_description.1`]

# Check
colnames(M.val) == E_GEOD_72680_Phenotypes$`!Sample_description.1`
colnames(betas2) == E_GEOD_72680_Phenotypes$`!Sample_description.1`

# Inspecting the data and preparing variables
E_GEOD_72680_Phenotypes$Batch = sapply(E_GEOD_72680_Phenotypes$`!Sample_description.1`, function(x){
  x = unlist(stri_split_fixed(x, pattern = "_"))
  x = x[1]
  return(x)
})
E_GEOD_72680_Phenotypes$Batch = factor(E_GEOD_72680_Phenotypes$Batch)
str(E_GEOD_72680_Phenotypes)
E_GEOD_72680_Phenotypes$age = as.numeric(E_GEOD_72680_Phenotypes$age)
E_GEOD_72680_Phenotypes$Sex = as.factor(E_GEOD_72680_Phenotypes$Sex)

# Fixing child abuse
E_GEOD_72680_Phenotypes$`childhood_sexual/physical_abuse_moderate_to_extreme`
E_GEOD_72680_Phenotypes$`childhood_sexual/physical_abuse_moderate_to_extreme` = ifelse(E_GEOD_72680_Phenotypes$`childhood_sexual/physical_abuse_moderate_to_extreme` == "--", NA,
                                                                                       E_GEOD_72680_Phenotypes$`childhood_sexual/physical_abuse_moderate_to_extreme`)
E_GEOD_72680_Phenotypes$`childhood_sexual/physical_abuse_moderate_to_extreme` = factor(E_GEOD_72680_Phenotypes$`childhood_sexual/physical_abuse_moderate_to_extreme`, 
                                                                                       levels = c("Yes", "No"))
# Fixing treatment_for_depression
E_GEOD_72680_Phenotypes$treatment_for_depression
E_GEOD_72680_Phenotypes$treatment_for_depression = ifelse(E_GEOD_72680_Phenotypes$treatment_for_depression == "--", NA,
                                                          E_GEOD_72680_Phenotypes$treatment_for_depression)
E_GEOD_72680_Phenotypes$treatment_for_depression = factor(E_GEOD_72680_Phenotypes$treatment_for_depression, levels = c("Yes", "No"))

# Fixing bipolar disorder
E_GEOD_72680_Phenotypes$treatment_for_bipolar_disorder
E_GEOD_72680_Phenotypes$treatment_for_bipolar_disorder = ifelse(E_GEOD_72680_Phenotypes$treatment_for_bipolar_disorder == "--", NA,
                                                                E_GEOD_72680_Phenotypes$treatment_for_bipolar_disorder)
E_GEOD_72680_Phenotypes$treatment_for_bipolar_disorder = factor(E_GEOD_72680_Phenotypes$treatment_for_bipolar_disorder, levels = c("Yes", "No"))

# Fixing posttraumatic_stress_disorder 
E_GEOD_72680_Phenotypes$treatment_for_posttraumatic_stress_disorder
E_GEOD_72680_Phenotypes$treatment_for_posttraumatic_stress_disorder = ifelse(E_GEOD_72680_Phenotypes$treatment_for_posttraumatic_stress_disorder == "--", NA,
                                                                             E_GEOD_72680_Phenotypes$treatment_for_posttraumatic_stress_disorder)
E_GEOD_72680_Phenotypes$treatment_for_posttraumatic_stress_disorder = factor(E_GEOD_72680_Phenotypes$treatment_for_posttraumatic_stress_disorder, levels = c("Yes", "No"))

# Fixing anxiety
E_GEOD_72680_Phenotypes$treatment_for_anxiety_disorder
E_GEOD_72680_Phenotypes$treatment_for_anxiety_disorder = ifelse(E_GEOD_72680_Phenotypes$treatment_for_anxiety_disorder == "--", NA,
                                                                E_GEOD_72680_Phenotypes$treatment_for_anxiety_disorder)
E_GEOD_72680_Phenotypes$treatment_for_anxiety_disorder = factor(E_GEOD_72680_Phenotypes$treatment_for_anxiety_disorder, levels = c("Yes", "No"))

# Numerical scores
E_GEOD_72680_Phenotypes[,9:21] = apply(E_GEOD_72680_Phenotypes[,9:21],2, function(x){
  x = str_trim(x)
  x = as.numeric(x)
  return(x)
})
E_GEOD_72680_Phenotypes$`race/ethnicity` = as.factor(E_GEOD_72680_Phenotypes$`race/ethnicity`)
E_GEOD_72680_Phenotypes[,24:30] = apply(E_GEOD_72680_Phenotypes[,24:30],2, function(x){
  x = str_trim(x)
  x = as.numeric(x)
  return(x)
})

# BDI score categorization
# https://www.ismanet.org/doctoryourspirit/pdfs/Beck-Depression-Inventory-BDI.pdf
# https://www.pnas.org/content/pnas/suppl/2019/05/20/1816847116.DCSupplemental/pnas.1816847116.sapp.pdf
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2905604/

# Inspecting
E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score
summary(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score)

# Making categories based on thresholds
# Standard
E_GEOD_72680_Phenotypes$Beck_depression_binary = sapply(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score, function(x){
  if (is.na(x)){
    return(NA)
  }
  if (x>=19){
    x = "High"
  } else {
    x = "Low"
  }
  return(x)
})
E_GEOD_72680_Phenotypes$Beck_depression_binary = factor(E_GEOD_72680_Phenotypes$Beck_depression_binary, levels = c("High","Low"))
table(E_GEOD_72680_Phenotypes$Beck_depression_binary)

# All categories
E_GEOD_72680_Phenotypes$Beck_depression_full = sapply(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score, function(x){
  if (is.na(x)){
    return(NA)
  }
  if (x %in% 0:10){
    return("Normal")
  }
  if (x %in% 11:16){
    return("Mild.mood.dist")
  }
  if (x %in% 17:20){
    return("Borderline.clin.depr")
  }
  if (x %in% 21:30){
    return("Mod.depr")
  }
  if (x %in% 31:40){
    return("Sev.depr")
  }
  if (x > 40){
    return("Extr.depr")
  }
})
table(E_GEOD_72680_Phenotypes$Beck_depression_full)
E_GEOD_72680_Phenotypes$Beck_depression_full = ordered(E_GEOD_72680_Phenotypes$Beck_depression_full, levels = c("Normal",
                                                                                                                "Mild.mood.dist",
                                                                                                                "Borderline.clin.depr",
                                                                                                                "Mod.depr",
                                                                                                                "Sev.depr",
                                                                                                                "Extr.depr"))
# Strict threshold (without borderline depression)
E_GEOD_72680_Phenotypes$Beck_depression_binary_strict = sapply(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score, function(x){
  if (is.na(x)){
    return(NA)
  }
  if (x>=21){
    x = "High"
  } else {
    x = "Low"
  }
  return(x)
})
E_GEOD_72680_Phenotypes$Beck_depression_binary_strict = factor(E_GEOD_72680_Phenotypes$Beck_depression_binary_strict, levels = c("High","Low"))
table(E_GEOD_72680_Phenotypes$Beck_depression_binary_strict)

# Threshold for severe depression
E_GEOD_72680_Phenotypes$Beck_depression_binary_severe = sapply(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score, function(x){
  if (is.na(x)){
    return(NA)
  }
  if (x>=31){
    x = "High"
  } else {
    x = "Low"
  }
  return(x)
})
E_GEOD_72680_Phenotypes$Beck_depression_binary_severe = factor(E_GEOD_72680_Phenotypes$Beck_depression_binary_severe, levels = c("High","Low"))
table(E_GEOD_72680_Phenotypes$Beck_depression_binary_severe)

# Threshold for extreme depression
E_GEOD_72680_Phenotypes$Beck_depression_binary_extreme = sapply(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score, function(x){
  if (is.na(x)){
    return(NA)
  }
  if (x>=41){
    x = "High"
  } else {
    x = "Low"
  }
  return(x)
})
E_GEOD_72680_Phenotypes$Beck_depression_binary_extreme = factor(E_GEOD_72680_Phenotypes$Beck_depression_binary_extreme, levels = c("High","Low"))
table(E_GEOD_72680_Phenotypes$Beck_depression_binary_extreme)

# Making composite depression categorization based on BDI category (strict) and different logic of interpreting missing values
E_GEOD_72680_Phenotypes$Composite_depression = NA

for (i in 1:nrow(E_GEOD_72680_Phenotypes)){
  if (is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i]) & is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i])){
    E_GEOD_72680_Phenotypes$Composite_depression[i] = NA
  } else if (!is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i]) & !is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i])){
    if (E_GEOD_72680_Phenotypes$treatment_for_depression[i] == "Yes"){
      E_GEOD_72680_Phenotypes$Composite_depression[i] = "Depressed"
    } else if (E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i]>=21){
      E_GEOD_72680_Phenotypes$Composite_depression[i] = "Depressed"
    } else {
      E_GEOD_72680_Phenotypes$Composite_depression[i] = "Normal"
    }
  } else if (!is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i]) & is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i])){
    if (E_GEOD_72680_Phenotypes$treatment_for_depression[i] == "Yes"){
      E_GEOD_72680_Phenotypes$Composite_depression[i] = "Depressed"
    } else {
      E_GEOD_72680_Phenotypes$Composite_depression[i] = "Normal"
    }
  } else if (is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i]) & !is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i])){
    if (E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i]>=21){
      E_GEOD_72680_Phenotypes$Composite_depression[i] = "Depressed"
    } else {
      E_GEOD_72680_Phenotypes$Composite_depression[i] = "Normal"
    }
  }
}

E_GEOD_72680_Phenotypes$Composite_depression_NA = NA

for (i in 1:nrow(E_GEOD_72680_Phenotypes)){
  if (is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i]) & is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i])){
    E_GEOD_72680_Phenotypes$Composite_depression_NA[i] = NA
  } else if (!is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i])){
    if (E_GEOD_72680_Phenotypes$treatment_for_depression[i] == "Yes"){
      E_GEOD_72680_Phenotypes$Composite_depression_NA[i] = "Depressed"
    } else if (E_GEOD_72680_Phenotypes$treatment_for_depression[i] == "No" & is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i])){
      E_GEOD_72680_Phenotypes$Composite_depression_NA[i] = NA
    } else if (E_GEOD_72680_Phenotypes$treatment_for_depression[i] == "No" & E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i]>=21){
      E_GEOD_72680_Phenotypes$Composite_depression_NA[i] = "Depressed"
    } else if (E_GEOD_72680_Phenotypes$treatment_for_depression[i] == "No" & E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i]<21){
      E_GEOD_72680_Phenotypes$Composite_depression_NA[i] = "Normal"
    }
  } else if (E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i]>=21){
    E_GEOD_72680_Phenotypes$Composite_depression_NA[i] = "Depressed"
  } else {
    E_GEOD_72680_Phenotypes$Composite_depression_NA[i] = "Normal"
  }
}

E_GEOD_72680_Phenotypes$Composite_depression_Strict = NA

for (i in 1:nrow(E_GEOD_72680_Phenotypes)){
  if (is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i]) | is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i])){
    E_GEOD_72680_Phenotypes$Composite_depression_Strict[i] = NA
  } else if (E_GEOD_72680_Phenotypes$treatment_for_depression[i] == "Yes" & E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i]>=21){
    E_GEOD_72680_Phenotypes$Composite_depression_Strict[i] = "Depressed"
  } else {
    E_GEOD_72680_Phenotypes$Composite_depression_Strict[i] = "Normal"
  }
}

E_GEOD_72680_Phenotypes$Composite_depression_NA_full = NA

for (i in 1:nrow(E_GEOD_72680_Phenotypes)){
  if (is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i]) & is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i])){
    E_GEOD_72680_Phenotypes$Composite_depression_NA_full[i] = NA
  } else if (!is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i]) & E_GEOD_72680_Phenotypes$treatment_for_depression[i] == "Yes"){
    E_GEOD_72680_Phenotypes$Composite_depression_NA_full[i] = "Depressed"
  } else if (!is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i]) & E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i]>=21) {
    E_GEOD_72680_Phenotypes$Composite_depression_NA_full[i] = "Depressed"
  } else if (is.na(E_GEOD_72680_Phenotypes$treatment_for_depression[i]) & E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i] < 21) {
    E_GEOD_72680_Phenotypes$Composite_depression_NA_full[i] =  NA
  } else if (is.na(E_GEOD_72680_Phenotypes$beck_depression_inventory_total_score[i]) & E_GEOD_72680_Phenotypes$treatment_for_depression[i] == "No"){
    E_GEOD_72680_Phenotypes$Composite_depression_NA_full[i] =  NA
  } else {
    E_GEOD_72680_Phenotypes$Composite_depression_NA_full[i] = "Normal"
  }
}

table(E_GEOD_72680_Phenotypes$Composite_depression, exclude = NULL)
table(E_GEOD_72680_Phenotypes$Composite_depression_NA, exclude = NULL)
table(E_GEOD_72680_Phenotypes$Composite_depression_Strict, exclude = NULL)
table(E_GEOD_72680_Phenotypes$Composite_depression_NA_full, exclude = NULL)

# Preparing factors for composite scores
E_GEOD_72680_Phenotypes$Composite_depression = factor(E_GEOD_72680_Phenotypes$Composite_depression, levels = c("Depressed", "Normal"))
E_GEOD_72680_Phenotypes$Composite_depression_NA = factor(E_GEOD_72680_Phenotypes$Composite_depression_NA, levels = c("Depressed", "Normal"))
E_GEOD_72680_Phenotypes$Composite_depression_Strict = factor(E_GEOD_72680_Phenotypes$Composite_depression_Strict, levels = c("Depressed", "Normal"))
E_GEOD_72680_Phenotypes$Composite_depression_NA_full = factor(E_GEOD_72680_Phenotypes$Composite_depression_NA_full, levels = c("Depressed", "Normal"))

# Fixing columns of the pheno data
colnames(E_GEOD_72680_Phenotypes) = stri_replace_all_fixed(colnames(E_GEOD_72680_Phenotypes), pattern = "!", replacement = "")
colnames(E_GEOD_72680_Phenotypes) = stri_replace_all_fixed(colnames(E_GEOD_72680_Phenotypes), pattern = "-", replacement = "_")
colnames(E_GEOD_72680_Phenotypes) = stri_replace_all_fixed(colnames(E_GEOD_72680_Phenotypes), pattern = "/", replacement = "_")

################### Final processing ###################

# Sample_description.1 initial column
# E_GEOD_72680_Phenotypes$Batch new batch column
# Batch effect
#unadj.M = ComBat(M.val, batch = E_GEOD_72680_Phenotypes$Batch) # Code produces error Error in while (change > conv) {: missing value where TRUE/FALSE needed

# Check
all(colnames(betas2) == E_GEOD_72680_Phenotypes$Sample_description.1)

unadj.beta = ComBat(betas2, batch = E_GEOD_72680_Phenotypes$Batch)
unadj.M = beta2m(unadj.beta) # Produced NAs

# NA stats
NA_in_betas = apply(unadj.M, 1, function(x) any(is.na(x)))
table(NA_in_betas) # 34039 CpGs have missing values -> removed
Clean_values = !NA_in_betas

# Selecting clean CpGs
unadj.beta = unadj.beta[Clean_values,]
unadj.M = unadj.M[Clean_values,]

write.table(unadj.beta, "unadjusted_beta_GSE72680.txt", sep = "\t")
write.table(unadj.M, "unadjusted_Mval_GSE72680.txt", sep = "\t")

################### Cell proportion adjustment ###################

# Cell proportions, EXPERIMENTAL
Cell_counts_blood = meffil.estimate.cell.counts.from.betas(Beta_from_signals_raw, cell.type.reference = "blood gse35069")
all(rownames(Cell_counts_blood) == E_GEOD_72680_Phenotypes$Sample_description.1) # Everything is matching
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
all(E_GEOD_72680_Phenotypes$Sample_description.1 == colnames(M_adj))

# Writing output
fwrite(adj.betas,file="adjusted_blood_beta_combined_GSE72680.txt", sep = "\t", row.names = TRUE) # 346687
fwrite(M_adj,"M_values_adj_combined_GSE72680.txt", sep = "\t", row.names = TRUE) # 346687
write.csv(E_GEOD_72680_Phenotypes, "Phenos_GSE72680.csv", row.names = FALSE)

