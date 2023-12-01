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
library(GEOquery)
library(meffil)

################### Defining functions ###################

# NOT IN operator
'%!in%' = function(x,y){!('%in%'(x,y))}

# A function to quickly read data
smart_fread = function(x){
  x = as.data.frame(fread(x, nThread = 10))
  rownames(x) = x$V1
  x$V1 = NULL
  return(x)
}

################### Importing data ###################
# Obtaining betas 2
betas2 = smart_fread("betas2.csv")
gc()
M.val = smart_fread("Mval.csv")
gc()
GSE198904_phenotypes_cohort_2 = smart_fread("GSE198904_phenotypes_cohort_2.csv")
gc()
Beta_from_signals_raw = smart_fread("Beta_from_signals_raw_GSE198904_DHRC.csv")
gc()

# Check 
all(colnames(betas2) == GSE198904_phenotypes_cohort_2$title) # TRUE
all(colnames(betas2) == colnames(M.val)) # TRUE
all(colnames(betas2) == colnames(Beta_from_signals_raw)) # FALSE, column order is not matching!

# Resetting order for Beta_from_signals_raw
Beta_from_signals_raw = Beta_from_signals_raw[, GSE198904_phenotypes_cohort_2$title]
all(colnames(betas2) == colnames(Beta_from_signals_raw)) # TRUE
GSE198904_phenotypes_cohort_2$Batch = factor(GSE198904_phenotypes_cohort_2$Batch)

################### Batch correction with combat ###################
# Batch effect
unadj.beta = ComBat(betas2, batch = GSE198904_phenotypes_cohort_2$Batch)
unadj.M = beta2m(unadj.beta) # Produced NAs

# NA stats
any(is.na(unadj.beta))
any(is.na(unadj.M))

# NA stats
NA_in_betas = apply(unadj.M, 1, function(x) any(is.na(x)))
table(NA_in_betas) # 5761 CpGs have missing values -> removed
Clean_values = !NA_in_betas

# Selecting clean CpGs
unadj.beta = unadj.beta[Clean_values,]
unadj.M = unadj.M[Clean_values,]

write.table(unadj.beta, "unadjusted_beta_GSE198904_DHRC.txt", sep = "\t")
write.table(unadj.M, "unadjusted_Mval_GSE198904_DHRC.txt", sep = "\t")
gc()

################### Cell counts ###################
# Cell proportions, EXPERIMENTAL
Cell_counts_blood = meffil.estimate.cell.counts.from.betas(as.matrix(Beta_from_signals_raw), cell.type.reference = "blood gse35069")
all(rownames(Cell_counts_blood) == GSE198904_phenotypes_cohort_2$title) # Everything is matching
setbun = as.data.frame(Cell_counts_blood)
write.table(setbun, file = "wbcset.txt", sep = "\t")
gc()


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
all(GSE198904_phenotypes_cohort_2$title == colnames(M_adj)) # TRUE

# Removal to free memory
rm(unadj.beta)
rm(unadj.M)
rm(M.val)
rm(betas2)
rm(x)
rm(y)
rm(beta_matrix)
gc()

# Writing output
fwrite(adj.betas,file="adj_beta_GSE198904_cohort2_full.txt", sep = "\t", row.names = TRUE) # 673011
fwrite(M_adj,"M_values_adj_GSE198904_cohort2_full.txt", sep = "\t", row.names = TRUE) # 673011
write.csv(GSE198904_phenotypes_cohort_2, "GSE198904_phenotypes_cohort_2_full.csv", row.names = FALSE)
gc()

###################### Selecting non-duplicated participants ###################### 
GSE198904_phenotypes_cohort_2_selected = arrange(GSE198904_phenotypes_cohort_2, geo_accession)
adj.betas_selected = adj.betas[,GSE198904_phenotypes_cohort_2_selected$title]
M_adj_selected = M_adj[,GSE198904_phenotypes_cohort_2_selected$title]

# Inspecting duplicated participants
duplicated_participants = duplicated(GSE198904_phenotypes_cohort_2_selected$subject)
duplicated_participants_names = GSE198904_phenotypes_cohort_2_selected$subject[duplicated_participants]
unique_participants = !duplicated_participants

# Selecting unique participants
GSE198904_phenotypes_cohort_2_selected = GSE198904_phenotypes_cohort_2_selected[unique_participants,]
adj.betas_selected = adj.betas_selected[,unique_participants]
M_adj_selected = M_adj_selected[,unique_participants]

# Check
all(unique(GSE198904_phenotypes_cohort_2$subject) %in% GSE198904_phenotypes_cohort_2_selected$subject)
all(colnames(adj.betas_selected) == GSE198904_phenotypes_cohort_2_selected$title) # TRUE
all(colnames(M_adj_selected) == GSE198904_phenotypes_cohort_2_selected$title) # TRUE

# Writing output
fwrite(adj.betas_selected,file="adj_beta_GSE198904_cohort2_selected.txt", sep = "\t", row.names = TRUE) # 673011
gc()
fwrite(M_adj_selected,"M_values_adj_GSE198904_cohort2_selected.txt", sep = "\t", row.names = TRUE) # 673011
gc()
write.csv(GSE198904_phenotypes_cohort_2_selected, "GSE198904_phenotypes_cohort_2_selected.csv", row.names = FALSE)
gc()