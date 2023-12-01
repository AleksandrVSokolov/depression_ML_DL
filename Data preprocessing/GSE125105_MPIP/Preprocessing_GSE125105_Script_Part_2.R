setwd("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE125105_MPIP")
options(stringsAsFactors = FALSE)


library(stringr)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(grid)
library(gdata)
library(RColorBrewer)
library(networkD3)
library(webshot)
library(htmlwidgets)
library(magrittr)
library(igraph)
library(visNetwork)
library(data.table)
library(XML)
library(rvest)
library(RCurl)
library(HGNChelper)
library(stringi)
library(httr)
library(lubridate)
library(rjson)
library(rtracklayer)
library(rstudioapi)
library(tidyr)
library(Gviz)
library(limma)
library(FactoMineR)
library(ggthemes)
library(igraph)
library(RSelenium)
library(lumi)
library(outliers)
library(svglite)
library(scatterplot3d)
library(sva)
library(jsonlite)
library(ggrepel)
library(parallel)
library(GEOquery)
library(BiocManager)

# In case of problems
library(settings)
reset_options()
options(help_type = "text")

###################### Functions ###################### 
'%!in%' = function(x,y){!('%in%'(x,y))}
smart_fread = function(x){
  x = as.data.frame(fread(x, nThread = 14))
  rownames(x) = x$V1
  x$V1 = NULL
  return(x)
}

################### Importing and preparing GSE125105 data ###################

# Reading data and performing curation of series matrix file
GSE125105_pheno_init = readLines("/home/aleksandr/Desktop/WORK/open_access_cohorts/GSE125105_Max_Plank_Depression/data_GSE125105") 
# the data file was downloaded from GEO (series matrix file)
GSE125105_pheno_init_small = substr(x = GSE125105_pheno_init, 1, 200)
GSE125105_pheno_selected = GSE125105_pheno_init[c(27:46,48,49,50,53,55,56,57,58,59,60,61,62,63,67)]
GSE125105_pheno_selected = stri_split_fixed(GSE125105_pheno_selected, pattern = "\t")
Names_list_GSE125105_pheno_selected = lapply(GSE125105_pheno_selected, function(x){
  Name = x[[1]]
  if (Name == "!Sample_characteristics_ch1"){
    Name_2 = x[[2]]
    Name_2 = unlist(stri_split_fixed(Name_2, pattern = ":"))
    Name_2 = Name_2[[1]]
    Name = Name_2
  }
  Name = stri_replace_all_fixed(Name, pattern = '"', replacement = "")
  Name = stri_replace_all_fixed(Name, pattern = '!', replacement = "")
  Name = stri_replace_all_fixed(Name, pattern = '-', replacement = "_")
  return(Name)
})
Names_list_GSE125105_pheno_selected = unlist(Names_list_GSE125105_pheno_selected)
GSE125105_pheno_selected = lapply(GSE125105_pheno_selected, function(x){
  x = x[-1]
  x = sapply(x, function(z){
    z = unlist(stri_split_fixed(z, pattern = ":"))
    if (length(z)>1){
      z = z[2]
    }
    return(z)
  })
  x = stri_replace_all_fixed(x, pattern = '"', replacement = "")
  return(x)
})
GSE125105_pheno_selected = do.call(cbind,GSE125105_pheno_selected)
GSE125105_pheno_selected = as.data.frame(GSE125105_pheno_selected)
Names_list_GSE125105_pheno_selected[duplicated(Names_list_GSE125105_pheno_selected)] = paste0(Names_list_GSE125105_pheno_selected[duplicated(Names_list_GSE125105_pheno_selected)], "_2")
colnames(GSE125105_pheno_selected) = Names_list_GSE125105_pheno_selected

# Reading targets file
GSE125105_target = as.data.frame(fread("/home/aleksandr/Desktop/WORK/open_access_cohorts/GSE125105_Max_Plank_Depression/Targets_GSE125105.csv"))
GSE125105_pheno_full = inner_join(GSE125105_target, GSE125105_pheno_selected, by = c("Participant" = "ID_REF"))

# Reading M-values
GSE125105_M_val_corrected = as.data.frame(fread("/home/aleksandr/Desktop/WORK/open_access_cohorts/GSE125105_Max_Plank_Depression/M_values_adj_combined_GSE125105.txt", 
                                                nThread = 10))
rownames(GSE125105_M_val_corrected) = GSE125105_M_val_corrected$V1
GSE125105_M_val_corrected$V1 = NULL

# Attaching calculated blood data (not really needed since values have been adjusted with regression)
GSE125105_blood = as.data.frame(fread("/home/aleksandr/Desktop/WORK/open_access_cohorts/GSE125105_Max_Plank_Depression/wbcset.txt", 
                                      nThread = 10))
GSE125105_pheno_full = inner_join(GSE125105_pheno_full, GSE125105_blood, by = c("Full_file_name" = "V1"))

# Check
all(colnames(GSE125105_M_val_corrected) == GSE125105_pheno_full$Full_file_name)
rm(list = c("GSE125105_pheno_init",
            "GSE125105_pheno_init_small",
            "GSE125105_pheno_selected", 
            "Names_list_GSE125105_pheno_selected", 
            "GSE125105_target", 
            "GSE125105_blood"))

################### Saving outputs ###################
fwrite(GSE125105_M_val_corrected, file = "GSE125105_M_val_corrected.csv",  row.names = TRUE)
fwrite(GSE125105_pheno_full, file = "GSE125105_pheno_full.csv", row.names = TRUE)




