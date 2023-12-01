
library(data.table)
library(limma)
library(parallel)
options(stringsAsFactors = FALSE)

writeLines("Running limma-based feature selection")

#########  Defining reading function ######### 
core_number = detectCores()

if (core_number <= 4){
  core_number = 2
} else {
  core_number = core_number - 3
}

smart_fread = function(x, ...){
  x = as.data.frame(fread(x, nThread = core_number, ...))
  return(x)
}

######### Getting parameters ######### 
ARGS = commandArgs(trailingOnly = TRUE)

folder_main = ARGS[1]
data_pheno_name = ARGS[2]
data_methylation_name = ARGS[3]
output_path = ARGS[4]

######### Setting working directory #########
setwd(folder_main)

######### Importing files ######### 
pheno_df = smart_fread(data_pheno_name)
mval_df = smart_fread(data_methylation_name)
mval_df = t(mval_df)

######### Preparing covariates #########
pheno_df$Sex = factor(pheno_df$Sex, levels = c("Male","Female"))
pheno_df$Depression = factor(pheno_df$Depression, levels = c("Control", "Case"))
pheno_df$Study = factor(pheno_df$Study)

######### Limma steps #########
# Design.matrix = model.matrix(~ Depression + Sex + Age + Study, data = pheno_df)
Design.matrix = model.matrix(~ Depression + Sex + Age + Study, data = pheno_df)
fit = lmFit(mval_df, Design.matrix)
fitE = eBayes(fit)
Top_table = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf)
Top_table$CpG = rownames(Top_table)

######### Writing result #########
fwrite(Top_table, file=output_path, sep = ",") 


