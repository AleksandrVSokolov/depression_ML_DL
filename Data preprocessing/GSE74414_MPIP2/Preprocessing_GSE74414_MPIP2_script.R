Working_directory = "/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE74414_MPIP2" # Replace with an appropriate directory
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
preprocess_GSE_var = function(x){
  x = sapply(x, function(z){
    z = stri_split_fixed(z, pattern = ":")
    z = unlist(z)
    z = z[2]
    z = str_trim(z)
    return(z)
  })
  return(x)
}
geo_data_table_gsm_extract = function(GEO_ID){
  Data_GEO = getGEOfile(
    GEO = GEO_ID,
    destdir = getwd(), amount = "full"
  )
  gunzip(Data_GEO)
  Data_GEO = stri_replace_all_fixed(Data_GEO, pattern = ".gz", replacement = "") #useful!
  Data_full = parseGEO(Data_GEO)
  file.remove(Data_GEO)
  GSMS = Data_full@gsms
  Full_GSMS = lapply(GSMS, function(x) x@dataTable@table)
  Full_Columns = lapply(GSMS, function(x) x@dataTable@columns)
  Full_Columns = Full_Columns[[1]]
  # for (i in 1:length(Full_GSMS)){
  #   colnames(Full_GSMS[[i]]) = paste0(names(Full_GSMS)[i], "_", colnames(Full_GSMS[[i]]))
  # }
  Full_GSMS = do.call(cbind, Full_GSMS)
  # colnames(Full_GSMS) = sapply(colnames(Full_GSMS), function(x){
  #   x = unlist(stri_split_fixed(x, pattern = "."))[2]
  #   return(x)
  # })
  Output = list()
  Output[[1]] = Full_Columns
  Output[[2]] = Full_GSMS
  return(Output)
}
geo_data_table_gsm_extract_phenotypes = function(GEO_ID){
  Data_GEO = getGEOfile(
    GEO = GEO_ID,
    destdir = getwd(), amount = "full"
  )
  gunzip(Data_GEO, remove = TRUE, overwrite = TRUE)
  Data_GEO = stri_replace_all_fixed(Data_GEO, pattern = ".gz", replacement = "") #useful!
  Data_full = parseGEO(Data_GEO)
  GSMS = Data_full@gsms
  file.remove(Data_GEO)
  Full_GSMS = lapply(GSMS, function(x) x@header)
  Colnames_GSMS = names(Full_GSMS[[1]])
  Full_GSMS_DFs = lapply(Full_GSMS, function(x){
    x = lapply(x, function(j){
      if (length(j) > 1){
        j = paste0(j, collapse = "___")
      }
      return(j)
    })
    x = data.frame(x)
    char_vector = x$characteristics_ch1
    char_vector = unlist(stri_split_fixed(char_vector, pattern =  "___"))
    vars_GSM = preprocess_GSE_var(char_vector)
    names_GSM = preprocess_GSE_var_names(char_vector)
    vars_GSM = as.data.frame(vars_GSM)
    vars_GSM = t(vars_GSM)
    vars_GSM = as.data.frame(vars_GSM)
    colnames(vars_GSM) = names_GSM
    x$characteristics_ch1 = NULL
    x$data_row_count = NULL
    x = x[1,]
    x = cbind(x, vars_GSM)
    rownames(x) = NULL
    return(x)
  })
  Data_pheno_GSMS = as.data.frame(do.call(rbind, Full_GSMS_DFs))
  return(Data_pheno_GSMS)
}

GSE74414 = getGEO("GSE74414")
GSE74414_data = GSE74414[[1]]
GSE74414_phenotypes = pData(GSE74414_data)
GSE74414_phenotypes_baseline = GSE74414_phenotypes[GSE74414_phenotypes$`treatment:ch1` == "baseline", ]
colnames(GSE74414_phenotypes_baseline) = stri_replace_all_fixed(colnames(GSE74414_phenotypes_baseline), pattern = ":ch1", replacement = "")
colnames(GSE74414_phenotypes_baseline) = stri_replace_all_fixed(colnames(GSE74414_phenotypes_baseline), pattern = " ", replacement = "_")
pheno_Janine_GSE74414 = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE74414_MPIP2/Data_from_Janine/GSE74414_baseline_pData_107.txt")
pheno_Janine_GSE74414$Sex_2 = ifelse(pheno_Janine_GSE74414$Sex == 2,"male", "female")
pheno_Janine_GSE74414$Status_2 = ifelse(pheno_Janine_GSE74414$Status == 1,"major depressive disorder", "control")
pheno_Janine_GSE74414$Batch = sapply(pheno_Janine_GSE74414$Meth_ID, function(x){
  x = unlist(stri_split_fixed(x, pattern = "_"))
  x = x[1]
  return(x)
})
N_batches = length(unique(pheno_Janine_GSE74414$Batch))

################### Matching participants with Janine's data ###################
GSE74414_phenotypes_baseline$Barcode = NA
# DOES NOT WORK PROPERLY as some participants are not matched by age...
for (i in 1:nrow(GSE74414_phenotypes_baseline)){
  current_sample = GSE74414_phenotypes_baseline[i,]
  
  # match
  match = pheno_Janine_GSE74414[pheno_Janine_GSE74414$age == current_sample$age,]
  match = match[match$Sex_2 == current_sample$gender,]
  match = match[match$Status_2 == current_sample$disease_state,]
  
  if (nrow(match) == 0){
    GSE74414_phenotypes_baseline$Barcode[i] = NA
    next
  }
  
  if (nrow(match) == 1){
    GSE74414_phenotypes_baseline$Barcode[i] = match$Meth_ID
    next
  } else {
    match = match[match$CD8T == current_sample$cd8_t_cells,]
  }
  
  if (nrow(match) == 1){
    GSE74414_phenotypes_baseline$Barcode[i] = match$Meth_ID
    next
  } else {
    match = pheno_Janine_GSE74414[pheno_Janine_GSE74414$age == current_sample$age,]
    match = match[match$Sex_2 == current_sample$gender,]
    match = match[match$Status_2 == current_sample$disease_state,]
    matrix_1 = match[,4:9]
    matrix_0 = current_sample[,c("cd8_t_cells", "cd4_t_cells", "natural_killer_cells", "b_cells", "monocytes", "granulocytes")]
    matrix_0 = t(apply(matrix_0, 1, as.numeric))
    colnames(matrix_0) = colnames(matrix_1)
    matrix_2 = rbind(matrix_1, matrix_0)
    
    distance_probes = dist(matrix_2)
    distance_probes_selected = as.matrix(distance_probes)
    distance_probes_selected = distance_probes_selected[current_sample$geo_accession, ]
    distance_probes_selected = distance_probes_selected[names(distance_probes_selected) != current_sample$geo_accession]
    best_match = names(distance_probes_selected)[which.min(distance_probes_selected)]
    match = match[best_match,]
    GSE74414_phenotypes_baseline$Barcode[i] = match$Meth_ID
  }
}
which(is.na(GSE74414_phenotypes_baseline$Barcode)) # 15 25 34 58 67
# Need exact matching from Janine or some other method to adjust for batches

################### Getting supplementary data ###################
# Getting suppl. files for GSE113725
getGEOSuppFiles("GSE74414")
files = list.files("GSE74414")
files = files[stri_detect_fixed(files, pattern = ".gz")]
files = paste0("GSE74414", "/", files)
lapply(files, gunzip)

# Reading files
GSE74414_methyl = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE74414_MPIP2/GSE74414/GSE74414_all_raw_methylated_signals.txt")
GSE74414_non_methyl = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE74414_MPIP2/GSE74414/GSE74414_all_raw_unmethylated_signals.txt")
GSE74414_det_p = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE74414_MPIP2/GSE74414/GSE74414_all_detection.p.values.txt")
all(colnames(GSE74414_methyl) == colnames(GSE74414_non_methyl)) # TRUE
all(colnames(GSE74414_methyl) == colnames(GSE74414_det_p))  # TRUE

# Resetting row names
rownames(GSE74414_methyl) = GSE74414_methyl$ID_REF
rownames(GSE74414_non_methyl) = GSE74414_non_methyl$ID_REF
rownames(GSE74414_det_p) = GSE74414_det_p$ID_REF
GSE74414_methyl$ID_REF = NULL
GSE74414_non_methyl$ID_REF = NULL
GSE74414_det_p$ID_REF = NULL
all(rownames(GSE74414_methyl) == rownames(GSE74414_non_methyl)) # TRUE
all(rownames(GSE74414_methyl) == rownames(GSE74414_det_p))  # TRUE
all(colnames(GSE74414_methyl) == colnames(GSE74414_non_methyl)) # TRUE
all(colnames(GSE74414_methyl) == colnames(GSE74414_det_p))  # TRUE

# Remove first rows
GSE74414_methyl = GSE74414_methyl[-1,]
GSE74414_non_methyl = GSE74414_non_methyl[-1,]
GSE74414_det_p = GSE74414_det_p[-1,]
ref_rownames = rownames(GSE74414_methyl)
ref_colnames = colnames(GSE74414_methyl)

# Set as numeric
GSE74414_methyl = as.data.frame(apply(GSE74414_methyl, 2, as.numeric))
GSE74414_non_methyl = as.data.frame(apply(GSE74414_non_methyl, 2, as.numeric))
GSE74414_det_p = as.data.frame(apply(GSE74414_det_p, 2, as.numeric))
rownames(GSE74414_methyl) = ref_rownames
rownames(GSE74414_non_methyl) = ref_rownames
rownames(GSE74414_det_p) = ref_rownames

colnames(GSE74414_methyl) = paste0(colnames(GSE74414_methyl), " Methylated signal")
colnames(GSE74414_non_methyl) = paste0(colnames(GSE74414_non_methyl), " Unmethylated signal")
merged_signals = cbind(GSE74414_non_methyl, GSE74414_methyl)
fwrite(merged_signals, "merged_signals.csv", row.names = TRUE, sep = ",")

# Making RAW GenomicMethylSet 
GSE74414_raw = readGEORawFile(filename = "merged_signals.csv",
                               Uname = "Unmethylated signal",
                               Mname = "Methylated signal",
                               sep = ",")

Beta_from_signals_raw = getBeta(GSE74414_raw, offset = 100)
fwrite(Beta_from_signals_raw, "Beta_from_signals_raw_GSE74414.csv", row.names = TRUE, sep = ",")

# Check
all(rownames(Beta_from_signals_raw) == rownames(GSE74414_methyl)) # Not matching ...

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
pvals = GSE74414_det_p
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
M.val = beta2m(betas2) # 380768

# Check
all(GSE74414_phenotypes$geo_accession %in% colnames(betas2)) # TRUE
all(GSE74414_phenotypes$geo_accession == colnames(betas2)) # FALSE
colnames(GSE74414_phenotypes) = stri_replace_all_fixed(colnames(GSE74414_phenotypes), pattern = ":ch1", replacement = "")
colnames(GSE74414_phenotypes) = stri_replace_all_fixed(colnames(GSE74414_phenotypes), pattern = " ", replacement = "_")

# Reordering
betas2 = betas2[,GSE74414_phenotypes$geo_accession]
M.val = M.val[,GSE74414_phenotypes$geo_accession]

# Check
all(GSE74414_phenotypes$geo_accession == colnames(betas2)) # TRUE
all(GSE74414_phenotypes$geo_accession == colnames(M.val)) # TRUE

################### Select baseline ###################
GSE74414_phenotypes_baseline_final = GSE74414_phenotypes[GSE74414_phenotypes$treatment == "baseline",]
betas2 = betas2[,GSE74414_phenotypes_baseline_final$geo_accession]
M.val = M.val[,GSE74414_phenotypes_baseline_final$geo_accession]

# Check
all(GSE74414_phenotypes_baseline_final$geo_accession == colnames(betas2)) # TRUE
all(GSE74414_phenotypes_baseline_final$geo_accession == colnames(M.val)) # TRUE

# Surrogate variables for batch
GSE74414_phenotypes_baseline_final$age = as.numeric(GSE74414_phenotypes_baseline_final$age)
GSE74414_phenotypes_baseline_final$gender = factor(GSE74414_phenotypes_baseline_final$gender, levels = c("male", "female"))
GSE74414_phenotypes_baseline_final$disease_state = factor(GSE74414_phenotypes_baseline_final$disease_state, levels = c("control", "major depressive disorder"))

mod = model.matrix(~ disease_state + age + gender , data=GSE74414_phenotypes_baseline_final)
n.sv = num.sv(betas2, mod, method="leek")
n.sv # 0 surrogate variable
# num.sv(M.val, mod, method="leek") -> produces error: Error in eigen(t(dats) %*% dats) : infinite or missing values in 'x'
# since it does not produce any significate variables, the batch effect might be little
par(las = 1) 
boxplot(betas2, horizontal = TRUE) # Only one sample seems to be a bit off...
par(las = 1) 
boxplot(M.val, horizontal = TRUE) # Only one  sample seems to be a bit off...

write.table(betas2, "unadjusted_beta_GSE74414.txt", sep = "\t")
write.table(M.val, "unadjusted_Mval_GSE74414.txt", sep = "\t")

# Cell proportions, EXPERIMENTAL
Cell_counts_blood = meffil.estimate.cell.counts.from.betas(Beta_from_signals_raw, cell.type.reference = "blood gse35069")
all(rownames(Cell_counts_blood) == GSE74414_phenotypes_baseline_final$geo_accession) # Different lengths
setbun = as.data.frame(Cell_counts_blood)
setbun = setbun[GSE74414_phenotypes_baseline_final$geo_accession,]
all(rownames(setbun) == colnames(betas2)) # TRUE
all(rownames(setbun) == GSE74414_phenotypes_baseline_final$geo_accession) # TRUE

write.table(setbun, file = "wbcset.txt", sep = "\t")

# Matching things
unadj.beta = betas2

################### Cell type curation ###################

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
all(GSE74414_phenotypes_baseline_final$geo_accession == colnames(M_adj)) # TRUE

# Check
all(rownames(adj.betas) == rownames(M_adj)) # TRUE

# Writing output
fwrite(adj.betas,file="adjusted_blood_beta_GSE74414_baseline.txt", sep = "\t", row.names = TRUE) # 380768
fwrite(M_adj,"M_values_adj_GSE74414_baseline.txt", sep = "\t", row.names = TRUE) # 380768
write.csv(GSE74414_phenotypes_baseline_final, "Phenos_GSE74414_baseline.csv", row.names = FALSE)

