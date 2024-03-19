Working_directory = "/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Main_analysis" # Replace with an appropriate path
setwd(Working_directory)

# Setting options
getOption("scipen") # Default number notation is 0
options(scipen=999)
options(stringsAsFactors = FALSE)

################### Package import ###################
# Importing packages (Note: Not all of them may be required)
library(fun)
library(stringr)
library(dplyr)
library(ggplot2)
library(openxlsx)
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
library(RSelenium)
library(lumi)
library(outliers)
library(svglite)
library(scatterplot3d)
library(sva)
library(jsonlite)
library(ggrepel)
library(parallel)
library(bacon)
library(gridExtra)
library(grid)
library(ggplotify)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(chromoMap)
library(RIdeogram)
library(ggVennDiagram)
library(FactoMineR)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(missMethyl)
library(preprocessCore)
library(metafor)
################### Defining functions ###################

# NOT IN operator
'%!in%' = function(x,y){!('%in%'(x,y))}

# A function to convert list to data frame
list_to_df = function(data_list){
  if (length(data_list) > 1){
    data_list = do.call(rbind, data_list)
  } else {
    data_list = data_list[[1]]
  }
  return(data_list)
}

# A function to replace multiple patterns by multiple replacements in a string
multiple_stri_replacer = function(string, pattern_vector, replacement_vector){
  
  # Pattern_vector and replacement_vector should have the same length
  for (i in 1:length(pattern_vector)){
    string = stri_replace_all_fixed(str = string, pattern = pattern_vector[i], replacement = replacement_vector[i])
  }
  return(string)
}

# A function to read text files fast; uses data.table::fread
smart_fread = function(x, ...){
  x = as.data.frame(fread(x, nThread = 10, ...))
  if ("V1" %in% colnames(x)){
    rownames(x) = x$V1
    x$V1 = NULL
  }
  return(x)
}

# A function to detect at least one pattern in a string
multiple_stri_detector = function(string, pattern_vector){
  output_list = list()
  for (i in 1:length(pattern_vector)){
    output_list[[i]] = stri_detect_fixed(str = string, pattern = pattern_vector[i])
  }
  output_list = do.call(rbind, output_list)
  apply(output_list, 2, any)
}

# A function to expand a data frame where several columns contain condensed cells with a specified separator
multiple_expander = function(df, cols_to_expand, pattern){
  #
  orig_colnames = colnames(df)
  df_modif = df[, cols_to_expand, drop = FALSE]
  df_const = df[,-cols_to_expand, drop = FALSE]
  orig_colnames_modif = colnames(df_modif)
  
  # Running expansion
  df_list = list()
  for (i in 1:nrow(df_const)){
    print(i)
    curr_df_const = df_const[i,, drop = FALSE]
    curr_df_modif = df_modif[i,, drop = FALSE]
    curr_df_modif = apply(curr_df_modif, 2, function(x) unlist(stri_split_fixed(x, pattern = pattern)), simplify = FALSE)
    
    if (length(cols_to_expand) > 1){
      curr_df_modif = do.call(cbind, curr_df_modif)
    } else {
      curr_df_modif = unlist(curr_df_modif)
    }
    
    if (is.matrix(curr_df_modif)){
      for (b in 1:nrow(curr_df_modif)){
        curr_df_const[b, ] = curr_df_const[1, ]
      }
    } else {
      for (b in 1:length(curr_df_modif)){
        curr_df_const[b, ] = curr_df_const[1, ]
      }
    }
    
    curr_df_combined = cbind(curr_df_const, curr_df_modif)
    
    if (length(cols_to_expand) == 1){
      colnames(curr_df_combined)[ncol(curr_df_combined)] = orig_colnames[cols_to_expand]
    }
    df_list[[i]] = curr_df_combined
  }
  
  if (length(df_list) > 1){
    df_list = do.call(rbind, df_list)
  } else {
    df_list = df_list[[1]]
  }
  
  df_list = df_list[, orig_colnames]
  return(df_list)
}

# A function to check gene symbols
# Importing NIH dataset
Homo_Sapiens_Gene_info_NIH = smart_fread("/home/aleksandr/Desktop/WORK/UCSC_ID_MAP/Homo_sapiens.gene_info") # https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/ (replace with an appropriate path)
Homo_Sapiens_Gene_info_NIH_expanded = multiple_expander(df = Homo_Sapiens_Gene_info_NIH, cols_to_expand = 5, pattern = "|")
check_gene_symbol_NIH = function(PRF_gene_symbols, PRF_ref_NIH_expanded, PRF_replace_NA_with_old = FALSE){
  PRF_gene_symbols_check = lapply(PRF_gene_symbols, function(x){
    if (x %in% PRF_ref_NIH_expanded$Symbol_from_nomenclature_authority){
      Curr_gene = x
      Approved = TRUE
      Suggested.Symbol = x
    } else if (x %in% PRF_ref_NIH_expanded$Symbol){
      RRF_df = PRF_ref_NIH_expanded[PRF_ref_NIH_expanded$Symbol == x,]
      Curr_gene = x
      Approved = FALSE
      Suggested.Symbol = RRF_df[,"Symbol_from_nomenclature_authority"]
      Suggested.Symbol = unique(Suggested.Symbol)
      Suggested.Symbol = Suggested.Symbol[Suggested.Symbol != "-"]
      if (length(Suggested.Symbol) >= 1){
        Suggested.Symbol = Suggested.Symbol[1]
        if (Suggested.Symbol == x){
          Approved = TRUE
        }
      } else {
        Suggested.Symbol = RRF_df[,"Symbol"]
        Suggested.Symbol = unique(Suggested.Symbol)
        Suggested.Symbol = Suggested.Symbol[1]
      }
    } else if (x %in% PRF_ref_NIH_expanded$Synonyms){
      RRF_df = PRF_ref_NIH_expanded[PRF_ref_NIH_expanded$Synonyms == x,]
      Curr_gene = x
      Approved = FALSE
      Suggested.Symbol = RRF_df[,"Symbol_from_nomenclature_authority"]
      Suggested.Symbol = unique(Suggested.Symbol)
      Suggested.Symbol = Suggested.Symbol[Suggested.Symbol != "-"]
      if (length(Suggested.Symbol) >= 1){
        Suggested.Symbol = Suggested.Symbol[1]
      } else {
        Suggested.Symbol = RRF_df[,"Symbol"]
        Suggested.Symbol = unique(Suggested.Symbol)
        Suggested.Symbol = Suggested.Symbol[1]
      }
    } else {
      Curr_gene = x
      Approved = FALSE
      Suggested.Symbol = NA
      if (PRF_replace_NA_with_old){
        Suggested.Symbol = x
      }
    }
    Dataset = data.frame(x = Curr_gene, Approved = Approved, Suggested.Symbol = Suggested.Symbol)
  })
  PRF_gene_symbols_check = list_to_df(PRF_gene_symbols_check)
  return(PRF_gene_symbols_check)
}

# A function to create demographic characteristics for a dataset
characterize_dataset_generelized_two_subgroups = function(dataset,
                                                          study_char,
                                                          contrast_col_number,
                                                          contrast_vector,
                                                          participants_col_number,
                                                          model_covariates_col_vector,
                                                          columns_to_characterise_vector,
                                                          Remove_NA_predictors = TRUE,
                                                          drop_P = TRUE,
                                                          simplif_P = 3){
  
  # Running Raw Statistics
  total_number = nrow(dataset)
  total_number_message = paste0("Initial dataset includes ", total_number, " participants")
  writeLines(total_number_message)
  first_row_output = data.frame(Name = total_number_message,Category.1 = "", Category.2 = "", P.value = NA)
  
  # Removing missing data
  colnames(dataset)[participants_col_number] = "Participant_ID_FUN"
  indeces = c(participants_col_number, contrast_col_number, model_covariates_col_vector)
  indeces = unique(indeces)
  Model_df = dataset[,indeces]
  if (Remove_NA_predictors){
    
    writeLines("Removing participants with missing data")
    Participants_missing_data = apply(Model_df, 1, function(x){
      x = as.character(x)
      
      if (any(is.na(x))){
        return(TRUE)
      } else{
        return(FALSE)
      }
      
    })
    Participants_missing_data = Model_df$Participant_ID_FUN[Participants_missing_data]
    
    if (length(Participants_missing_data) < 1){
      excluded_number = 0
      
    } else {
      
      excluded_number = length(Participants_missing_data)
    }
    
    Participants_excluded_df = dataset[dataset$Participant_ID_FUN %in% Participants_missing_data, ]
    dataset = dataset[dataset$Participant_ID_FUN %!in% Participants_missing_data, ]
    excluded_message = paste0("Participants with missing data excluded: ", excluded_number,"\n",
                              "Resulting number of participants: ", nrow(dataset))
    second_row_output = data.frame(Name = excluded_message,Category.1 = "", Category.2 = "", P.value = NA)
    writeLines(excluded_message)
    
  } else {
    
    writeLines("Exclusion of participants was not performed")
    second_row_output = data.frame(Name = "Exclusion of participants was not performed", Category.1 = "", Category.2 = "", P.value = NA)
    Participants_excluded_df = NA
  }
  
  main_indeces = c(contrast_col_number, columns_to_characterise_vector)
  main_indeces = unique(main_indeces)
  Main_dataset = list()
  
  
  # Making categories
  Initial_names = colnames(dataset)[main_indeces]
  colnames(dataset)[contrast_col_number] = "variable_to_split"
  levels_var = dataset[,"variable_to_split"]
  levels_var = levels(levels_var)
  names(levels_var) = NULL
  Category_1 = levels_var[1]
  Category_2 = levels_var[2]
  
  for (i in 1:length(main_indeces)){
    curr_variable = dataset[, main_indeces[i]]
    if (is.factor(curr_variable)){
      Characterised_var = characterize_categorical_variable(df = dataset,
                                                            Category_1 = Category_1, 
                                                            Category_2 = Category_2, 
                                                            Variable_name = colnames(dataset)[main_indeces[i]], 
                                                            keep_missing = TRUE)
    } else if (is.character(curr_variable)){
      stop("All columns should be either a Factor or Numeric")
    } else if (is.numeric(curr_variable)){
      Characterised_var = characterize_numeric_variable(df = dataset, 
                                                        Category_1 = Category_1, 
                                                        Category_2 = Category_2, 
                                                        Variable_name = colnames(dataset)[main_indeces[i]], 
                                                        Mention_NAs = TRUE)
    }
    Current_DF = data.frame(Name = Initial_names[i], 
                            Category.1 = Characterised_var[[1]], 
                            Category.2 = Characterised_var[[2]], 
                            P.value = Characterised_var[[3]])
    Main_dataset[[i]] = Current_DF
  }
  Main_dataset = list_to_df(data_list = Main_dataset)
  Main_dataset$P.value = round(Main_dataset$P.value, digits = simplif_P)
  
  # Compiling full dataframe
  header = data.frame(Name = study_char, Category.1 = "", Category.2 = "", P.value = NA)
  Full_table = rbind(header, first_row_output, second_row_output, Main_dataset)
  
  if(drop_P){
    Full_table$P.value = NULL
  }
  
  output = list()
  output$Table = Full_table
  output$Excluded = Participants_excluded_df
  return(output)
}
# Several helper functions for characterization (used internally in characterize_dataset_generelized_two_subgroups)
characterize_categorical_variable = function(df, Category_1, Category_2, Variable_name, P.val.valid = TRUE, keep_missing = FALSE){
  Total_variable = df[, Variable_name]
  
  if (!is.ordered(Variable_name)){
    df[, Variable_name] = ordered(Total_variable, levels = sort(unique(Total_variable)))
    Total_variable = df[, Variable_name]
  }
  
  Contingency_table = table(df[,Variable_name], df[,"variable_to_split"])
  
  if (keep_missing){
    Contingency_table = table(df[,Variable_name], df[,"variable_to_split"], exclude = NULL)
  }
  
  Contingency_table_percents_1 = (Contingency_table[,Category_1]/(sum(Contingency_table[,Category_1])) * 100) %>% round(., digits = 1)
  Contingency_table_percents_2 = (Contingency_table[,Category_2]/(sum(Contingency_table[,Category_2])) * 100) %>% round(., digits = 1)
  Report_1 = paste0(names(Contingency_table_percents_1), ": ")
  Report_1 = paste0(Report_1, Contingency_table[,Category_1], " (",Contingency_table_percents_1, "%)", collapse = "\n")
  Report_2 = paste0(names(Contingency_table_percents_2), ": ")
  Report_2 = paste0(Report_2, Contingency_table[,Category_2], " (",Contingency_table_percents_2, "%)", collapse = "\n")
  
  Contingency_table_pval = table(df[,Variable_name], df[,"variable_to_split"])
  if (P.val.valid){
    tryCatch({
      Var.pval <<- chisq.test(Contingency_table_pval)
      Var.pval <<- Var.pval$p.value}, warning = function(w){
        message = paste0(Variable_name, " has small counts in some groups, using P.val from Fisher's Exact Test")
        writeLines(message)
        Var.pval <<- fisher.test(Contingency_table_pval)
        Var.pval <<- Var.pval$p.value
      })
  } else {
    Var.pval = "Not valid"
  }
  
  Output_list = list()
  Output_list[[1]] = Report_1
  Output_list[[2]] = Report_2
  Output_list[[3]] = Var.pval
  return(Output_list)
}

characterize_numeric_variable = function(df, Category_1, Category_2, Variable_name, Mention_NAs = FALSE){
  Total_variable = df[, Variable_name]
  Variable_1 = df[df$variable_to_split == Category_1, Variable_name]
  Variable_1_report = describe_vector_numeric(Variable_1, Mention_NAs = Mention_NAs)
  Variable_2 = df[df$variable_to_split == Category_2, Variable_name]
  Variable_2_report = describe_vector_numeric(Variable_2, Mention_NAs = Mention_NAs)
  Shapiro_normality_check = shapiro.test(Total_variable)
  Shapiro_normality_check = Shapiro_normality_check$p.value
  if (Shapiro_normality_check < 0.05){
    message = paste0(Variable_name, " is not normally distributed. Using P-val from Mann Whitney U Test")
    writeLines(message)
    formula_test = paste0(Variable_name, "~ variable_to_split")
    Total_pval = wilcox.test(formula = as.formula(formula_test), data = df)
    Total_pval = Total_pval$p.value
  } else {
    message = paste0(Variable_name, " is normally distributed. Using P-val from T Test")
    writeLines(message)
    formula_test = paste0(Variable_name, "~ variable_to_split")
    Var_test_check = var.test(formula	= as.formula(formula_test), data = df)
    Var_test_check = Var_test_check$p.value 
    if (Var_test_check >= 0.05){
      Total_pval = t.test(formula = as.formula(formula_test), data = df, var.equal = TRUE)
    } else  {
      Total_pval = t.test(formula = as.formula(formula_test), data = df, var.equal = FALSE)
    }
    Total_pval = Total_pval$p.value
  }
  Output_list = list()
  Output_list[[1]] = Variable_1_report
  Output_list[[2]] = Variable_2_report
  Output_list[[3]] = Total_pval
  return(Output_list)
}

describe_vector_numeric = function(x, Mention_NAs, show_as_percent = FALSE, Describe_min_max = TRUE){
  Mean_var = mean(x, na.rm = TRUE) %>% round(., digits = 2)
  SD_var = sd(x, na.rm = TRUE) %>% round(., digits = 2)
  Var_report = paste0(Mean_var, " ± ", SD_var)
  if (Mention_NAs){
    Missing_count = length(x[is.na(x)])
    Missing_val_percent = (Missing_count/length(x) * 100) %>% round(., digits = 2)
    if (Missing_count > 0){
      Var_report = paste0(Var_report, "\nMissing val: ", Missing_count, " (", Missing_val_percent, "%)")
    }
  }
  if (show_as_percent){
    Mean_var = Mean_var * 100
    SD_var = SD_var * 100
    Var_report = paste0(Mean_var, " ± ", SD_var, "%")
  }
  if (Describe_min_max){
    Min_var = min(x, na.rm = TRUE) %>% round(., digits = 2)
    Max_var = max(x, na.rm = TRUE) %>% round(., digits = 2)
    if (show_as_percent){
      Min_var = Min_var * 100
      Max_var = Max_var * 100
      Var_report = paste0(Var_report, "\n", "Min: ", Min_var, "%, Max: ",Max_var, "%")
    } else {
      Var_report = paste0(Var_report, "\n", "Min: ", Min_var, ", Max: ",Max_var)
    }
  }
  return(Var_report)
}

plot_cpg = function(CpG_name, 
                    B_values,
                    Phenotypes_df, 
                    X_lab = "Cohorts", 
                    Y_lab = "Beta value: ", 
                    file_path = NA,
                    text_x_vjust = 0.4){
  
  Bvals = B_values[CpG_name, ]
  Bvals = as.numeric(as.character(Bvals))
  tmp_df = Phenotypes_df
  tmp_df$CpG = Bvals
  plot = ggplot(data = tmp_df, aes(x = Study, y = CpG, fill = Depression)) +
    geom_boxplot(width = 0.3) +
    scale_fill_manual(values = c("#FC8D62","#66C2A5"))+
    labs(x = X_lab, y = paste0(Y_lab, CpG_name), fill = "Depression status") + 
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17), # Customizing the plot title
      axis.title = element_text(face = "bold", size = 14), # Customizing the axes titles
      legend.title = element_text(face = "bold", size = 14), # Customizing the legend title
      axis.text = element_text(size = 13, colour = "black"), # Customizing the axes text
      legend.text = element_text(size = 13, colour = "black"), # Customizing the legend text
      panel.background = element_blank(), # Removing ugly gray background
      axis.line = element_line(size = 1), # Adding axis lines
      axis.text.x = element_text(size = 13, colour = "black", angle = 45, vjust = text_x_vjust),
      panel.grid.major.y = element_line(size = 0.5, linetype = 3, color =  "black") # Modifying horizontal lines in the plot
    )
  
  if (!is.na(file_path)){
    ggsave(plot = plot, filename = file_path, width = 10, height = 10, dpi = 300)
  } else {
    return(plot)
  }
}

# A function to make Venn diagrams from named lists (uses ggplot2 and ggVennDiagram)
make_Venn_digram_list = function(named_list, plot_full_path = NULL, ...){
  
  # Customizable Venn diagram
  venn = Venn(named_list)
  data = process_data(venn)
  data@region$full_lable = sapply(data@region$count, function(x){
    Number = x
    Percent = x/sum(data@region$count)
    Percent = Percent*100
    Percent = round(Percent, digits = 1)
    Label_full = paste0(Number, "\n","(", Percent, "%)")
    return(Label_full)
  })
  
  # to see available shapes plot_shapes()
  plot = ggplot() +
    
    # 1. region count layer
    geom_sf(aes(fill = id), data = venn_region(data), ...) +
    
    # 2. set edge layer
    geom_sf(color="black", size = 0.5, data = venn_setedge(data), show.legend = FALSE, ...) +
    
    # 3. set label layer
    geom_sf_text(aes(label = name), data = venn_setlabel(data), ...) +
    
    # 4. region label layer
    geom_sf_label(aes(label = full_lable), data = venn_region(data), ...) +
    scale_fill_brewer(...) +
    scale_x_continuous(expand = expansion(mult = .2)) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.position = "none",
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      plot.background = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(1,1,1,1, "cm")
    )
  print(plot)
  if (is.character(plot_full_path)){
    pdf(plot_full_path, width = 10, height = 10)
    print(plot)
    dev.off()
  }
}


# This function creates the Manhattan plot and saves it in the specified folder (used internally in the function above)
make_manhattan_plot_toptable_generalized = function(Toptable_df, 
                                                    calculate_cumul_pos = TRUE, 
                                                    plot_name, 
                                                    plot_folder = FALSE, 
                                                    target_genes = FALSE){
  Toptable_df = arrange(Toptable_df, CpG_cumul_pos)
  
  if (calculate_cumul_pos){
    writeLines("Calculating CpGs cumul. pos.")
    Toptable_df$CpG_cumul_pos = mapply(function(x,y){
      output = calculate_cumulative_coordinates(chrom = x, pos = y)
      return(output)
    }, Toptable_df$Chromosome, Toptable_df$Position_Hg19)
    Toptable_df = arrange(Toptable_df, CpG_cumul_pos)
  }
  
  Toptable_df$Chromosome = stri_replace_all_fixed(Toptable_df$Chromosome, pattern = "chr", replacement = "")
  Toptable_df$Chromosome = as.numeric(Toptable_df$Chromosome)
  
  axisdf = list()
  for (i in 1:length(unique(Toptable_df$Chromosome))){
    Chrom = unique(Toptable_df$Chromosome)[i]
    Current_df = Toptable_df[Toptable_df$Chromosome == Chrom,]
    Current_df$CpG_cumul_pos = as.double(Current_df$CpG_cumul_pos)
    center = (max(Current_df$CpG_cumul_pos) + min(Current_df$CpG_cumul_pos))/2
    data_small = data.frame(Chromosome = Chrom, center = center)
    axisdf[[i]] = data_small
  }
  axisdf = list_to_df(axisdf)
  
  # Adding highlight for CpGs based on the gene
  if (is.character(target_genes)){
    writeLines("Adding highlight")
    Toptable_df = mutate(Toptable_df, is.highlight = ifelse(multiple_stri_detector(CpG_gene, pattern_vector = target_genes), 
                                                            "yes", "no"))
  }
  
  # Modification of CpGs DF
  Pval_treshold = Toptable_df[Toptable_df$adj.P.Val < 0.05,]
  Pval_treshold = arrange(Pval_treshold, -P.Value)
  Pval_treshold = Pval_treshold$P.Value[1]
  Pval_treshold = -log10(Pval_treshold)
  Toptable_df = mutate(Toptable_df, is_annotate= ifelse(-log10(Toptable_df$P.Value)>=Pval_treshold, "yes", "no"))
  
  # Make the plot
  plot = ggplot(Toptable_df, aes(x=CpG_cumul_pos, y=-log10(P.Value))) +
    
    # Show all points
    geom_point( aes(color=as.factor(Chromosome)), alpha=0.6, size=4) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$Chromosome, breaks = axisdf$center) +
    scale_y_continuous(expand = c(0, 1)) +     # Remove space between plot area and x axis
    
    # Add p-val line
    geom_hline(yintercept=Pval_treshold, linetype="dashed", 
               color = "red", size=0.5)
  
  if (is.character(target_genes)){
    # Add highlighted points
    plot = plot + geom_point(data=filter(Toptable_df, is.highlight =="yes"), color="orange", size = 6, alpha= 0.8)
  }
  
  plot = plot + 
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data = subset(Toptable_df, is_annotate=="yes"), aes(label = CpG), size = 8, force = 10, 
                     max.overlaps = 50) +
    
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 24),
      axis.title.x = element_text(size = 28, face = "bold"),
      axis.title.y = element_text(size = 28, face = "bold")
    ) +
    labs(x = "Chromosome")
  
  # Saving the plot
  writeLines("Saving plot")
  if (is.character(plot_folder)){
    plot_name = paste0(plot_folder, "/", plot_name)
  }
  png(plot_name, width = 4096, height = 2160, units = "px")
  print(plot)
  dev.off()
  
}

# This function calculates cumulated genomic coordinates for CpGs (used internally in the function test_CpGs_limma_generalized)
Chromosome_coord_table = smart_fread("Chromosome_coord_table.csv") # helper dataset for the function (replace with an appropriate path)
calculate_cumulative_coordinates = function(Coordinates = NULL, chrom = NULL, pos = NULL, pos_full = FALSE){
  if (!is.null(chrom) & !is.null(pos)){
    pos = as.double(pos)
    Curr_chrom = stri_replace_all_fixed(chrom, pattern = 'chr', replacement = '')
    if (Curr_chrom != "1"){
      Index = which(Chromosome_coord_table$Chromosome == Curr_chrom)
      Curr_table = Chromosome_coord_table[1:(Index-1), ]
      Combined_length = sum(Curr_table$`Total length (bp)`)
      pos = pos + Combined_length
    }
    if (pos_full){
      pos = paste0(chrom, ":", pos)
    }
    return(pos)
  }
  
  if (!is.null(Coordinates)){
    Coordinates = unlist(stri_split_fixed(Coordinates, pattern = c(':')))
    Coordinates = unlist(stri_split_fixed(Coordinates, pattern = c('-')))
    Curr_chrom = Coordinates[1]
    Curr_chrom = stri_replace_all_fixed(Curr_chrom, pattern = 'chr', replacement = '')
    if (Curr_chrom == "1"){
      Start = Coordinates[2]
      End = Coordinates[3]
    } else {
      Start = as.numeric(Coordinates[2])
      End = as.numeric(Coordinates[3])
      Index = which(Chromosome_coord_table$Chromosome == Curr_chrom)
      Curr_table = Chromosome_coord_table[1:(Index-1), ]
      Combined_length = sum(Curr_table$`Total length (bp)`)
      Start = Start + Combined_length
      End = End + Combined_length
    }
    coord_string = paste0(Start,"-", End)
    return(coord_string)
  }
  
}

################### Importing data ###################

# PSY SCR
PSY_SCR_mval = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/PSY_screening/M_values_adj_221_PSY_screen_comb_prepr.txt")
PSY_SCR_pheno = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/PSY_screening/Screening_pheno_small.csv")
PSY_SCR_illum = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/PSY_screening/SAMPLE_ID_FOR_PSY_screen_comb_prepr.csv")
PSY_SCR_pheno = inner_join(PSY_SCR_illum, PSY_SCR_pheno, by = c("Participant" = "Code"))
PSY_SCR_pheno = PSY_SCR_pheno[PSY_SCR_pheno$barcode %in% colnames(PSY_SCR_mval),]
all(colnames(PSY_SCR_mval) == PSY_SCR_pheno$barcode) # TRUE

# PSY Recall
PSY_RC_mval = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/PSY_recall/PSY_RC_Mval_selected.csv")
PSY_RC_pheno = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/PSY_recall/PSY_RC_Pheno.csv")
all(colnames(PSY_RC_mval) == PSY_RC_mval$barcode) # TRUE

# GSE125105_MPIP
GSE125105_MPIP_mval = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE125105_MPIP/GSE125105_M_val_corrected.csv")
GSE125105_MPIP_pheno = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE125105_MPIP/GSE125105_pheno_full.csv")
all(colnames(GSE125105_MPIP_mval) == GSE125105_MPIP_pheno$Full_file_name) # TRUE

# GSE72680_GRADY
GSE72680_GRADY_mval = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE72680_GRADY/M_values_adj_combined_GSE72680.txt")
GSE72680_GRADY_pheno = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE72680_GRADY/Phenos_GSE72680.csv")
all(colnames(GSE72680_GRADY_mval) == GSE72680_GRADY_pheno$Sample_description.1) # TRUE

# GSE113725_RDE
GSE113725_RDE_mval = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE113725_RDE/M_values_adj_GSE113725_selected.txt")
GSE113725_RDE_pheno = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE113725_RDE/Phenos_GSE113725_selected.csv")
all(colnames(GSE113725_RDE_mval) == GSE113725_RDE_pheno$title) # TRUE

# GSE198904_DHRC
GSE198904_DHRC_mval = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE198904_DHRC/M_values_adj_GSE198904_cohort2_selected.txt")
GSE198904_DHRC_pheno = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE198904_DHRC/GSE198904_phenotypes_cohort_2_selected.csv")
all(colnames(GSE198904_DHRC_mval) == GSE198904_DHRC_pheno$title) # TRUE

# GSE198904_OBSERVEMDD
GSE198904_OBS_mval = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE198904_OBSERVEMDD/M_values_adj_GSE198904_cohort1_selected.txt")
GSE198904_OBS_pheno = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE198904_OBSERVEMDD/GSE198904_phenotypes_cohort_1_selected.csv")
GSE198904_OBS_pheno = GSE198904_OBS_pheno[GSE198904_OBS_pheno$title %in% colnames(GSE198904_OBS_mval),]
all(colnames(GSE198904_OBS_mval) == GSE198904_OBS_pheno$title) # TRUE

# GSE74414_MPIP2
GSE74414_MPIP2_mval = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE74414_MPIP2/M_values_adj_GSE74414_baseline.txt")
GSE74414_MPIP2_pheno = smart_fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Data preprocessing/GSE74414_MPIP2/Phenos_GSE74414_baseline.csv")
all(colnames(GSE74414_MPIP2_mval) == GSE74414_MPIP2_pheno$geo_accession) # TRUE

# Filtering GSE74414_MPIP2
GSE74414_MPIP2_pheno[GSE74414_MPIP2_pheno$title == "Baseline subject 75",] # GSM1920137 is duplicated with GSE125105_MPIP
GSE74414_MPIP2_pheno = GSE74414_MPIP2_pheno[GSE74414_MPIP2_pheno$geo_accession %!in% "GSM1920137", ]
GSE74414_MPIP2_mval = GSE74414_MPIP2_mval[, GSE74414_MPIP2_pheno$geo_accession]

################### Curating phenotypes ###################
# PSY
PSY_SCR_pheno_selected = PSY_SCR_pheno[,c("barcode","Gender","Age","DAWBA_DEPBAND")]
colnames(PSY_SCR_pheno_selected) = c("ID", "Sex", "Age", "Depression")
PSY_SCR_pheno_selected$Study = "PSY_SCR"
PSY_SCR_pheno_selected$Depression = ifelse(PSY_SCR_pheno_selected$Depression>= 4, "Case", "Control")
PSY_SCR_pheno_selected$Sex = ifelse(PSY_SCR_pheno_selected$Sex == "W", "Female", "Male")

# PSY Recall
PSY_RC_pheno_selected = PSY_RC_pheno[, c("Participant","Gender..0.male..1.female","Age","DAWBA_DEPBAND")]
colnames(PSY_RC_pheno_selected) = c("ID", "Sex", "Age", "Depression")
PSY_RC_pheno_selected$Study = "PSY_RC"
PSY_RC_pheno_selected$Depression = ifelse(PSY_RC_pheno_selected$Depression>= 4, "Case", "Control")
PSY_RC_pheno_selected$Sex = ifelse(PSY_RC_pheno_selected$Sex == 1, "Female", "Male")

# GSE125105_MPIP
GSE125105_MPIP_pheno_small = GSE125105_MPIP_pheno[,c("Full_file_name","Sex","age","diagnosis")]
colnames(GSE125105_MPIP_pheno_small) = c("ID", "Sex", "Age", "Depression")
GSE125105_MPIP_pheno_small$Study = "GSE125105_MPIP"
GSE125105_MPIP_pheno_small$Depression = ifelse(GSE125105_MPIP_pheno_small$Depression == "case", "Case", "Control")
GSE125105_MPIP_pheno_small$Sex = ifelse(GSE125105_MPIP_pheno_small$Sex == "F", "Female", "Male")

# GSE72680_GRADY
GSE72680_GRADY_pheno_small = GSE72680_GRADY_pheno[,c("Sample_description.1","Sex", "age", "Composite_depression_NA_full")]
colnames(GSE72680_GRADY_pheno_small) = c("ID", "Sex", "Age", "Depression")
GSE72680_GRADY_pheno_small = GSE72680_GRADY_pheno_small[!apply(GSE72680_GRADY_pheno_small, 1, function(x) any(is.na(x))), ]
GSE72680_GRADY_pheno_small$Study = "GSE72680_GRADY"
GSE72680_GRADY_pheno_small$Depression = ifelse(GSE72680_GRADY_pheno_small$Depression == "Depressed", "Case", "Control")
GSE72680_GRADY_mval = GSE72680_GRADY_mval[,GSE72680_GRADY_pheno_small$ID]
all(colnames(GSE72680_GRADY_mval) == GSE72680_GRADY_pheno_small$ID) # TRUE

# GSE113725_RDE
GSE113725_RDE_pheno_small = GSE113725_RDE_pheno[,c("title", "gender", "age", "Disease" )]
colnames(GSE113725_RDE_pheno_small) = c("ID", "Sex", "Age", "Depression")
GSE113725_RDE_pheno_small$Study = "GSE113725_RDE"
GSE113725_RDE_pheno_small$Sex = ifelse(GSE113725_RDE_pheno_small$Sex == "F", "Female", "Male")

# GSE198904_DHRC
GSE198904_DHRC_pheno_small = GSE198904_DHRC_pheno[,c("title", "gender", "age", "diagnosis")]
colnames(GSE198904_DHRC_pheno_small) = c("ID", "Sex", "Age", "Depression")
GSE198904_DHRC_pheno_small$Study = "GSE198904_DHRC"
GSE198904_DHRC_pheno_small$Sex = ifelse(GSE198904_DHRC_pheno_small$Sex == "F", "Female", "Male")
GSE198904_DHRC_pheno_small$Depression = ifelse(GSE198904_DHRC_pheno_small$Depression == "MDD", "Case", "Control")

# GSE198904_OBS
GSE198904_OBS_pheno_small = GSE198904_OBS_pheno[,c("title", "gender", "age", "diagnosis")]
colnames(GSE198904_OBS_pheno_small) = c("ID", "Sex", "Age", "Depression")
GSE198904_OBS_pheno_small$Study = "GSE198904_OBS"
GSE198904_OBS_pheno_small$Sex = ifelse(GSE198904_OBS_pheno_small$Sex == "F", "Female", "Male")
GSE198904_OBS_pheno_small$Depression = ifelse(GSE198904_OBS_pheno_small$Depression == "OBSERVEMDD0001", "Case", "Control")

# GSE74414_MPIP2
GSE74414_MPIP2_pheno_small = GSE74414_MPIP2_pheno[,c("geo_accession", "gender", "age", "disease_state")]
colnames(GSE74414_MPIP2_pheno_small) = c("ID", "Sex", "Age", "Depression")
GSE74414_MPIP2_pheno_small$Study = "GSE74414_MPIP2"
GSE74414_MPIP2_pheno_small$Sex = ifelse(GSE74414_MPIP2_pheno_small$Sex == "female", "Female", "Male")
GSE74414_MPIP2_pheno_small$Depression = ifelse(GSE74414_MPIP2_pheno_small$Depression == "major depressive disorder", "Case", "Control")

################### Getting combined dataset ###################
Combined_pheno = rbind(
  PSY_SCR_pheno_selected,
  PSY_RC_pheno_selected,
  GSE125105_MPIP_pheno_small,
  GSE72680_GRADY_pheno_small,
  GSE113725_RDE_pheno_small,
  GSE198904_DHRC_pheno_small,
  GSE198904_OBS_pheno_small,
  GSE74414_MPIP2_pheno_small
)
# 1945 samples in total !
table(Combined_pheno$Depression) # 1128 cases # 817 controls

################### Inspecting overlapping methylation ###################

# Selecting CpGs
CpGs_PSY_SCR = rownames(PSY_SCR_mval)
CpGs_PSY_RC = rownames(PSY_RC_mval)
CpGs_GSE125105_MPIP = rownames(GSE125105_MPIP_mval)
CpGs_GSE72680_GRADY = rownames(GSE72680_GRADY_mval)
CpGs_GSE113725_RDE = rownames(GSE113725_RDE_mval)
CpGs_GSE198904_DHRC = rownames(GSE198904_DHRC_mval)
CpGs_GSE198904_OBS = rownames(GSE198904_OBS_mval)
CpGs_GSE74414_MPIP2 = rownames(GSE74414_MPIP2_mval)

CpG_list = list(
  "PSY Screen" = CpGs_PSY_SCR,
  "PSY Recall" = CpGs_PSY_RC,
  "GSE125105 MPIP" = CpGs_GSE125105_MPIP,
  "GSE72680 GRADY" = CpGs_GSE72680_GRADY,
  "GSE113725 RDE" = CpGs_GSE113725_RDE,
  "GSE198904 DHRC" = CpGs_GSE198904_DHRC,
  "GSE198904 OBS" = CpGs_GSE198904_OBS,
  "GSE74414 MPIP2" = CpGs_GSE74414_MPIP2
)

# Number of sets is too large to visualize by Venn diagrams
Intersecting_CpGs = Reduce(intersect, CpG_list)
sapply(CpG_list, function(x) all(Intersecting_CpGs %in% x)) # All -> TRUE
# 304765 CpGs are presented in all cohorts and arrays

Combined_Mval = list(
  PSY_SCR_mval[Intersecting_CpGs,],
  PSY_RC_mval[Intersecting_CpGs,],
  GSE125105_MPIP_mval[Intersecting_CpGs,],
  GSE72680_GRADY_mval[Intersecting_CpGs,],
  GSE113725_RDE_mval[Intersecting_CpGs,],
  GSE198904_DHRC_mval[Intersecting_CpGs,],
  GSE198904_OBS_mval[Intersecting_CpGs,],
  GSE74414_MPIP2_mval[Intersecting_CpGs,]
)

# Check
sapply(Combined_Mval, function(x) all(rownames(x) == Intersecting_CpGs)) # All -> TRUE

Combined_Mval = do.call(cbind, Combined_Mval)
Combined_Mval = Combined_Mval[,Combined_pheno$ID]

# Check
all(colnames(Combined_Mval) == Combined_pheno$ID) # TRUE

################### Saving merged data ###################

if (SAVING_merged){
  fwrite(Combined_Mval,file="Combined_Mval_non_adj.csv", sep = ",", row.names = TRUE) 
  fwrite(Combined_pheno,file="Combined_pheno_init.csv", sep = ",") 
  gc()
} else {
  Combined_Mval = smart_fread("Combined_Mval_non_adj.csv")
  Combined_pheno = smart_fread("Combined_pheno_init.csv")
}

################### Exploring cohort-related batching ###################
color_vector = sapply(colnames(Combined_Mval), function(x){
  study = Combined_pheno[Combined_pheno$ID == x, "Study"]
  lookup = list(
    "GSE113725_RDE" = "#845EC2",
    "GSE125105_MPIP" = "#D65DB1",
    "GSE198904_DHRC" = "#C34A36",
    "GSE198904_OBS" = "#FF8066",
    "GSE72680_GRADY" = "#4E8397",
    "GSE74414_MPIP2" = "#FEFEDF",
    "PSY_RC" = "#FFC75F",
    "PSY_SCR" = "#ADC5CF"
  )
  color = lookup[[study]]
  return(color)
})
png(filename = "cohort_batching.png", width = 5000, height = 2000, units = "px")
boxplot(Combined_Mval, axes = FALSE, col=color_vector)
dev.off()
# Needs adjustment by cohort covariate!
gc()

################### Preprocessing with QN and Combat ###################

PERFORM_QN = FALSE

if (PERFORM_QN){
  
  Combined_Mval_QN =  as.data.frame(normalize.quantiles(as.matrix(Combined_Mval)))
  
  colnames(Combined_Mval_QN) = colnames(Combined_Mval)
  all(colnames(Combined_Mval_QN) == Combined_pheno$ID)
  
  Combined_Mval_Combat_QN = sva::ComBat(Combined_Mval_QN, batch = Combined_pheno$Study)
  gc()
  
  png(filename = "cohort_batching_QN_combat.png", width = 5000, height = 2000, units = "px")
  boxplot(Combined_Mval_Combat_QN, axes = FALSE, col=color_vector)
  dev.off()
  gc()
  
  Combined_Mval_Combat_QN = as.data.frame(Combined_Mval_Combat_QN)
  rownames(Combined_Mval_Combat_QN) = rownames(Combined_Mval)
  fwrite(Combined_Mval_Combat_QN, file="Combined_Mval_Combat_QN.csv", sep = ",", row.names = TRUE) 
  gc()
  
} else {
  print("Correction has been performed earlier -> Importing data")
  
  Combined_Mval_Combat_QN = smart_fread("Combined_Mval_Combat_QN.csv")
  
  # Check
  all(colnames(Combined_Mval_Combat_QN) == Combined_pheno$ID) # TRUE
  gc()
}

################### Performing PCA for cohorts (corrected data) ###################

PERFORM_PCA = FALSE

if (PERFORM_PCA){
  # Hypervariable CpGs
  Combined_betas_Combat = m2beta(Combined_Mval_Combat_QN)
  p97.5 = apply(Combined_betas_Combat, 1, function(x) quantile(x, probs=.975))
  p2.5 = apply(Combined_betas_Combat, 1, function(x) quantile(x, probs=.025))
  CpG_diff = p97.5-p2.5
  CpGs_hypervar = names(CpG_diff[CpG_diff > 0.2])
  
  # Performing PCA
  matrix_CpGs = as.data.frame(t(Combined_Mval_Combat_QN[CpGs_hypervar,])) # Use M-values as input!
  PCA_fit = PCA(matrix_CpGs, ncp = 10, graph = FALSE)
  PC_1 = PCA_fit$ind$coord[,1]
  PC_2 = PCA_fit$ind$coord[,2]
  
  # Plotting PCs
  Combined_pheno$PC_1 = PC_1
  Combined_pheno$PC_2 = PC_2
  plot = ggplot(data = Combined_pheno, aes(x = PC_1, y = PC_2, col = Study)) +
    geom_point() +
    theme_bw()
  ggsave(plot = plot, filename = "PCA_plotCombat.pdf", width = 10, height = 10, dpi = 300)
  
  plot = ggplot(data = Combined_pheno, aes(x = PC_1, y = PC_2, col = Depression)) +
    geom_point() +
    theme_bw()
  ggsave(plot = plot, filename = "PCA_plotCombat_vs_Depression.pdf", width = 10, height = 10, dpi = 300)
  
  # Seems like cohorts are not clustered -> PC coefficients may not be needed
}

################### Preparing Illumina annotation ###################
PREPARE_ANNOT = FALSE

if (PREPARE_ANNOT){
  ILLUM_450K_ANNOT = as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
  ILLUM_450K_ANNOT = ILLUM_450K_ANNOT[rownames(Combined_Mval_Combat_no_NA), ]
  ILLUM_450K_ANNOT$Updated_gene_name = sapply(ILLUM_450K_ANNOT$UCSC_RefGene_Name, function(x){
    if (x == ""){
      return(NA)
    }
    x = stri_split_fixed(x, pattern = ";")
    x = unlist(x)
    x = unique(x)
    x = x[x != ""]
    new_symbols = check_gene_symbol_NIH(PRF_gene_symbols = x,
                                        PRF_ref_NIH_expanded = Homo_Sapiens_Gene_info_NIH_expanded, 
                                        PRF_replace_NA_with_old = TRUE)
    new_symbols = new_symbols$Suggested.Symbol
    new_symbols = paste(new_symbols, collapse = ";")
    return(new_symbols)
  })
  ILLUM_450K_ANNOT$CpG_cumul_pos = mapply(function(x,y){
    output = calculate_cumulative_coordinates(chrom = x, pos = y)
    return(output)
  }, ILLUM_450K_ANNOT$chr, ILLUM_450K_ANNOT$pos)
  
  fwrite(ILLUM_450K_ANNOT, file="ILLUM_450K_ANNOT_fixed.csv", sep = ",", row.names = TRUE)
  
} else {
  
  ILLUM_450K_ANNOT = smart_fread("ILLUM_450K_ANNOT_fixed.csv")
}

################### Cleaning datasets ###################
Combined_pheno$Sex = factor(Combined_pheno$Sex, levels = c("Male","Female"))
Combined_pheno$Depression = factor(Combined_pheno$Depression, levels = c("Control", "Case"))
Combined_pheno$Study = factor(Combined_pheno$Study)

# Removing missing things
NA_ind = apply(Combined_pheno, 1, function(x) any(is.na(x)))
Combined_pheno[NA_ind,] # 201172200044_R08C01; 201172200025_R08C01 excluded from Control GSE198904_DHRC
Combined_pheno_no_NA = Combined_pheno[!NA_ind,]
Combined_Mval_Combat_no_NA = Combined_Mval_Combat_QN[, Combined_pheno_no_NA$ID]
Combined_Betas_Combat_no_NA = m2beta(Combined_Mval_Combat_no_NA)

Combined_Mval_no_NA = Combined_Mval[,Combined_pheno_no_NA$ID]

################### Making subsets ###################
# Generating train and test sets
CV_pheno = Combined_pheno_no_NA[Combined_pheno_no_NA$Study %in% c("GSE125105_MPIP", "GSE198904_DHRC", "GSE72680_GRADY", "PSY_SCR", "GSE113725_RDE"),]
CV_mval = Combined_Mval_Combat_no_NA[,CV_pheno$ID]

test_pheno = Combined_pheno_no_NA[Combined_pheno_no_NA$ID %!in% CV_pheno$ID,]
test_mval = Combined_Mval_Combat_no_NA[,test_pheno$ID]

# Check 
all(CV_pheno$ID == colnames(CV_mval)) # TRUE
all(test_pheno$ID == colnames(test_mval)) # TRUE

################### Saving data for ML ###################
# Preparing illumina annotation
ILLUM_450K_ANNOT_selected_ML = ILLUM_450K_ANNOT
ILLUM_450K_ANNOT_selected_ML = arrange(ILLUM_450K_ANNOT_selected_ML, CpG_cumul_pos)

# Re-arranging CpGs
CV_mval = CV_mval[rownames(ILLUM_450K_ANNOT_selected_ML),]
test_mval = test_mval[rownames(ILLUM_450K_ANNOT_selected_ML),]

# Writing data
CV_pheno = CV_pheno[,c("ID","Sex","Age","Depression","Study" )]
test_pheno = test_pheno[,c("ID","Sex","Age","Depression","Study" )]
fwrite(CV_pheno, file="ML_data/CV_pheno.csv", sep = ",", row.names = FALSE) 
fwrite(test_pheno, file="ML_data/test_pheno.csv", sep = ",", row.names = FALSE)
#
fwrite(CV_mval, file="ML_data/CV_mval.csv", sep = ",", row.names = TRUE)
fwrite(test_mval, file="ML_data/test_mval.csv", sep = ",", row.names = TRUE)
gc()

# Writing Annotation
fwrite(ILLUM_450K_ANNOT_selected_ML, file="ML_data/ILLUM_450K_ANNOT_selected_ML.csv", sep = ",", row.names = TRUE)
gc()


################### Saving RAW data for ML ###################

# Generating train and test sets
CV_pheno = Combined_pheno_no_NA[Combined_pheno_no_NA$Study %in% c("GSE125105_MPIP", "GSE198904_DHRC", "GSE72680_GRADY", "PSY_SCR", "GSE113725_RDE"),]
CV_mval_raw = Combined_Mval_no_NA[,CV_pheno$ID]

test_pheno = Combined_pheno_no_NA[Combined_pheno_no_NA$ID %!in% CV_pheno$ID,]
test_mval_raw = Combined_Mval_no_NA[,test_pheno$ID]

# Preparing illumina annotation
ILLUM_450K_ANNOT_selected_ML = ILLUM_450K_ANNOT
ILLUM_450K_ANNOT_selected_ML = arrange(ILLUM_450K_ANNOT_selected_ML, CpG_cumul_pos)

# Re-arranging CpGs
CV_mval_raw = CV_mval_raw[rownames(ILLUM_450K_ANNOT_selected_ML),]
test_mval_raw = test_mval_raw[rownames(ILLUM_450K_ANNOT_selected_ML),]

# Writing data
CV_pheno = CV_pheno[,c("ID","Sex","Age","Depression","Study" )]
test_pheno = test_pheno[,c("ID","Sex","Age","Depression","Study" )]
fwrite(CV_pheno, file="ML_data_raw/CV_pheno.csv", sep = ",", row.names = FALSE) 
fwrite(test_pheno, file="ML_data_raw/test_pheno.csv", sep = ",", row.names = FALSE)
#
fwrite(CV_mval_raw, file="ML_data_raw/CV_mval.csv", sep = ",", row.names = TRUE)
fwrite(test_mval_raw, file="ML_data_raw/test_mval.csv", sep = ",", row.names = TRUE)
gc()

# Writing Annotation
fwrite(ILLUM_450K_ANNOT_selected_ML, file="ML_data_raw/ILLUM_450K_ANNOT_selected_ML.csv", sep = ",", row.names = TRUE)
gc()

################### Plotting RAW CV and test intensities ###################
CV_mval_raw = smart_fread("ML_data_raw/CV_mval.csv")
png(filename = "CV_data_RAW_box.png", width = 5000, height = 2000, units = "px")
boxplot(CV_mval_raw, axes = FALSE)
dev.off()

test_mval_raw = smart_fread("ML_data_raw/test_mval.csv")
png(filename = "test_data_RAW_box.png", width = 5000, height = 2000, units = "px")
boxplot(test_mval_raw, axes = FALSE)
dev.off()

################### Plotting Preprocessed CV and test intensities ###################
CV_mval = smart_fread("ML_data/CV_mval.csv")
png(filename = "CV_data_harmonized_box.png", width = 5000, height = 2000, units = "px")
boxplot(CV_mval, axes = FALSE)
dev.off()

test_mval= smart_fread("ML_data/test_mval.csv")
png(filename = "test_data_harmonized_box.png", width = 5000, height = 2000, units = "px")
boxplot(test_mval, axes = FALSE)
dev.off()

# Checks
all(colnames(CV_mval) == colnames(CV_mval_raw)) # TRUE
all(colnames(test_mval) == colnames(test_mval_raw)) # TRUE
all(rownames(CV_mval) == rownames(CV_mval_raw)) # TRUE
all(rownames(test_mval) == rownames(test_mval_raw)) # TRUE


################### Differential methylation analysis with Limma (full cohort) ###################
PERFORM_DM = FALSE

if (PERFORM_DM){
  
  # Limma code
  Design.matrix = model.matrix(~ Depression + Sex + Age + Study, data = Combined_pheno_no_NA)
  fit = lmFit(Combined_Mval_Combat_no_NA, Design.matrix)
  fitE = eBayes(fit)
  Top_table = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf)
  Top_table$CpG = rownames(Top_table)
  
  # Annotating results
  ILLUM_450K_ANNOT_selected = ILLUM_450K_ANNOT[Top_table$CpG, ]
  all(rownames(ILLUM_450K_ANNOT_selected) == rownames(Top_table)) # TRUE
  Top_table$Gene = ILLUM_450K_ANNOT_selected$Updated_gene_name
  Top_table$Enchancer = ifelse(ILLUM_450K_ANNOT_selected$Enhancer== "TRUE", "Yes", "No")
  Top_table$Chromosome = ILLUM_450K_ANNOT_selected$chr
  Top_table$Pos_hg19 = ILLUM_450K_ANNOT_selected$pos
  
  Top_table_signif = Top_table[Top_table$P.Value < 0.05,]
  Top_table_signif$Cohort_agreement = unlist(invisible(mclapply(Top_table_signif$CpG, function(x){
    tmp_df = Combined_pheno_no_NA
    tmp_df$CpG = as.numeric(t(Combined_Mval_Combat_no_NA[x,]))
    CpG_stat = tmp_df %>%
      group_by(., Study, Depression) %>%
      summarize(., Median = median(CpG))
    CpG_stat = CpG_stat %>% 
      spread(Depression, Median)
    CpG_stat$Difference = CpG_stat$Case - CpG_stat$Control
    CpG_stat$Direction = ifelse(CpG_stat$Difference > 0, "Positive", "Negative")
    CpG_stat$Direction = factor(CpG_stat$Direction, levels = c("Positive", "Negative"))
    CpG_stat_freq = as.data.frame(table(CpG_stat$Direction))
    CpG_stat_freq_pos = CpG_stat_freq[CpG_stat_freq$Var1 == "Positive", "Freq"]
    CpG_stat_freq_negat = CpG_stat_freq[CpG_stat_freq$Var1 == "Negative", "Freq"]
    if (CpG_stat_freq_pos > CpG_stat_freq_negat){
      CpG_stat_str = paste0("Pos: ", CpG_stat_freq_pos, "/", CpG_stat_freq_pos+CpG_stat_freq_negat)
    } else if (CpG_stat_freq_pos < CpG_stat_freq_negat){
      CpG_stat_str = paste0("Negat: ", CpG_stat_freq_negat, "/", CpG_stat_freq_pos+CpG_stat_freq_negat)
    } else {
      CpG_stat_str = paste0((CpG_stat_freq_pos+CpG_stat_freq_negat)/2, "/", CpG_stat_freq_pos+CpG_stat_freq_negat)
    }
    return(CpG_stat_str)
  }, mc.cores = 10)))
  Top_table_signif$Agreement_index = sapply(Top_table_signif$Cohort_agreement, function(x){
    x = stri_replace_all_fixed(x, pattern = "Pos: ", replacement = "")
    x = stri_replace_all_fixed(x, pattern = "Negat: ", replacement = "")
    x = eval(parse(text = x))
    return(x)
  })
  Top_table_signif_consistent = Top_table_signif[Top_table_signif$Agreement_index == 1,]
  
  # Writing Files
  fwrite(ILLUM_450K_ANNOT_selected,file="ILLUM_450K_ANNOT_selected.csv", sep = ",") 
  fwrite(Top_table,file="Top_table.csv", sep = ",") 
  fwrite(Top_table_signif,file="Top_table_signif.csv", sep = ",")
  fwrite(Top_table_signif_consistent,file="Top_table_signif_consistent.csv", sep = ",") 
  
} else {
  
  ILLUM_450K_ANNOT_selected = smart_fread("ILLUM_450K_ANNOT_selected.csv")
  Top_table = smart_fread("Top_table.csv")
  Top_table_signif = smart_fread("Top_table_signif.csv")
  Top_table_signif_consistent = smart_fread("Top_table_signif_consistent.csv")
  
  rownames(ILLUM_450K_ANNOT_selected) = ILLUM_450K_ANNOT_selected$Name
  rownames(Top_table) = Top_table$CpG
  rownames(Top_table_signif) = Top_table_signif$CpG
  rownames(Top_table_signif_consistent) = Top_table_signif_consistent$CpG
  
}

################### Plotting Top hits ###################
top_CpGs_plots = lapply(Top_table$CpG[1:5], function(x){
  x = plot_cpg(CpG_name = x, 
               B_values = Combined_Betas_Combat_no_NA, 
               Phenotypes_df = Combined_pheno_no_NA,
               text_x_vjust = 0.5)
  return(x)
})
pdf(file = "Top_CpGs_vs_Cohorts.pdf", width = 20, height = 50)
gridExtra::grid.arrange(grobs = top_CpGs_plots, nrow = 5, ncol = 1)
dev.off()

top_CpGs_plots = lapply(Top_table_signif_consistent$CpG[1:5], function(x){
  x = plot_cpg(CpG_name = x, 
               B_values = Combined_Betas_Combat_no_NA, 
               Phenotypes_df = Combined_pheno_no_NA,
               text_x_vjust = 0.5)
  return(x)
})
pdf(file = "Top_consistent_CpGs_vs_Cohorts.pdf", width = 20, height = 50)
gridExtra::grid.arrange(grobs = top_CpGs_plots, nrow = 5, ncol = 1)
dev.off()

################### Plotting Manhattan and Volcano Plot ###################

MAKE_MAIN_PLOTS = FALSE

if (MAKE_MAIN_PLOTS){
  
  # Manhattan Plot
  Top_table_manhattan = Top_table
  Top_table_manhattan$CpG_cumul_pos = NA
  colnames(Top_table_manhattan) = stri_replace_all_fixed(colnames(Top_table_manhattan), pattern =  "Pos_hg19", replacement = "Position_Hg19")
  make_manhattan_plot_toptable_generalized(Toptable_df = Top_table_manhattan, 
                                           calculate_cumul_pos = TRUE, 
                                           plot_name = "Manhattan_all_CpGs.png")
  
  # Volcano Plot
  Top_table_vlc = Top_table_manhattan
  Pval_treshold = Top_table_manhattan[Top_table_manhattan$adj.P.Val < 0.05,]
  Pval_treshold = max(Pval_treshold$P.Value)
  Pval_treshold = -log10(Pval_treshold)
  if (is.na(Pval_treshold)){
    Pval_treshold = NULL
  }
  PREFIX_logFC_threshold = 0.05
  
  Top_table_vlc$is.highlight = sapply(Top_table_vlc$logFC, function(x){
    if (x > PREFIX_logFC_threshold){
      x = "Up"
    } else if (x < - PREFIX_logFC_threshold){
      x = "Down"
    } else {
      x = "None"
    }
  })
  
  if (!is.null(Pval_treshold)){
    Top_table_vlc = mutate(Top_table_vlc, is_annotate = ifelse(-log10(P.Value) >= Pval_treshold & is.highlight != "None", "yes", "no"))
  } else {
    Top_table_vlc$is_annotate = "no"
  }
  
  # Make the plot
  plot = ggplot(Top_table_vlc, aes(x=logFC, y=-log10(P.Value))) +
    
    # Show all points
    geom_point(aes(color= factor(is.highlight, levels = c("None", "Down", "Up"))), alpha=0.6, size=4) +
    scale_color_manual(values = c("grey", "skyblue", "red")) 
  
  # Add standard scale for plot
  plot = plot + scale_x_continuous(limits = c(-0.5,0.5))
  
  # Add pval line
  if (!is.null(Pval_treshold)){
    plot = plot + geom_hline(yintercept=Pval_treshold, linetype="dashed", 
                             color = "red", size=0.5)
  }
  
  # Add pval line 0.05
  plot = plot + geom_hline(yintercept= -log10(0.05), linetype="dashed", 
                           color = "blue", size=0.5)
  # Add logFC lines
  if (min(Top_table_vlc$logFC) < -PREFIX_logFC_threshold){
    plot = plot + geom_vline(xintercept= -PREFIX_logFC_threshold, linetype="dashed", 
                             color = "grey", size=0.5)
  }
  
  if (max(Top_table_vlc$logFC) > PREFIX_logFC_threshold){
    plot = plot + geom_vline(xintercept= PREFIX_logFC_threshold, linetype="dashed", 
                             color = "grey", size=0.5)
  }
  
  # Add label using ggrepel to avoid overlapping
  if (any(Top_table_vlc$is_annotate == "yes")){
    plot = plot + 
      geom_label_repel(data=subset(Top_table_vlc, is_annotate=="yes"), aes(label=CpG), size=4, force = 10, 
                       max.overlaps = 50)
  } 
  
  plot = plot +
    labs(x = "log2 fold change", y = "-log10 p-value") +
    # Custom the theme:
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(size = 0.1, linetype = 2, color =  "black"), # Modifying horizontal lines in the plot
      panel.background = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold")
    )
  ggsave(file = "Volcano_All_CpGs.png", plot = plot, width=2560, height=1440, units = "px", scale = 2)
}

################### Saving Top 200 ALL consistent CpGs ###################
Top_200_all_consistent_CpGs = Top_table_signif_consistent$CpG[1:200]
Top_200_all_signif_CpGs = Top_table_signif$CpG[1:200]
writeLines(text = Top_200_all_consistent_CpGs, con = "ML_data/Top_200_all_consistent_CpGs.txt")
writeLines(text = Top_200_all_consistent_CpGs, con = "ML_data_raw/Top_200_all_consistent_CpGs.txt")
make_Venn_digram_list(list(Top_200_all_consistent_CpGs = Top_200_all_consistent_CpGs, Top_200_all_signif_CpGs = Top_200_all_signif_CpGs))

################### Plotting hits ###################
plot_cpg(CpG_name = "cg18991165", 
         B_values = Combined_Betas_Combat_no_NA, 
         Phenotypes_df = Combined_pheno_no_NA,
         text_x_vjust = 0.5, file_path = "test_plot.pdf")

################### Differential methylation analysis with Limma (individual cohorts) ###################

PERFORM_individ_analysis = FALSE

if (PERFORM_individ_analysis){
  
  # Performing analysis
  studies = unique(Combined_pheno_no_NA$Study)
  Diff_Expr_Individ = lapply(studies, function(x){
    
    writeLines(paste0("Working with study: ", x))
    tmp_pheno = Combined_pheno_no_NA[Combined_pheno_no_NA$Study == x,]
    tmp_mval = Combined_Mval_no_NA[,tmp_pheno$ID]
    
    Design.matrix = model.matrix(~ Depression + Sex + Age, data = tmp_pheno)
    fit = lmFit(tmp_mval, Design.matrix)
    fitE = eBayes(fit)
    Top_table_tmp = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
    Top_table_tmp$CpG = rownames(Top_table_tmp)
    # *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975
    SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
    SE = SE[,2]
    SE = SE[Top_table_tmp$CpG]
    Top_table_tmp$SE = SE
    
    # Annotating results
    ILLUM_450K_ANNOT_tmp = ILLUM_450K_ANNOT[Top_table_tmp$CpG, ]
    Top_table_tmp$Gene = ILLUM_450K_ANNOT_tmp$Updated_gene_name
    Top_table_tmp$Enchancer = ifelse(ILLUM_450K_ANNOT_tmp$Enhancer== "TRUE", "Yes", "No")
    Top_table_tmp$Chromosome = ILLUM_450K_ANNOT_tmp$chr
    Top_table_tmp$Pos_hg19 = ILLUM_450K_ANNOT_tmp$pos
    
    Top_table_tmp$Study = x
    
    return(Top_table_tmp)
    
  })
  Diff_Expr_Individ = do.call(rbind, Diff_Expr_Individ)
  Diff_Expr_Individ_signif = Diff_Expr_Individ[Diff_Expr_Individ$P.Value < 0.05,]
  Diff_Expr_Individ_fdr = Diff_Expr_Individ[Diff_Expr_Individ$adj.P.Val < 0.05,]
  
  fwrite(Diff_Expr_Individ, "Diff_Expr_Individ.csv", sep = ",")
  
} else {
  
  Diff_Expr_Individ = smart_fread("Diff_Expr_Individ.csv")
  Diff_Expr_Individ_signif = Diff_Expr_Individ[Diff_Expr_Individ$P.Value < 0.05,]
  Diff_Expr_Individ_fdr = Diff_Expr_Individ[Diff_Expr_Individ$adj.P.Val < 0.05,]
  
}

fwrite(Diff_Expr_Individ_signif, "Diff_Expr_Individ_signif.csv", sep = ",")


PERFORM_CpG_STATs = FALSE

if (PERFORM_CpG_STATs){
  # CpG stats
  studies = unique(Diff_Expr_Individ$Study)
  CpG_stat_individ_meth = mclapply(unique(Diff_Expr_Individ_signif$CpG), function(x){
    
    Curr_df = Diff_Expr_Individ_signif[Diff_Expr_Individ_signif$CpG == x,]
    Curr_df_fdr = Diff_Expr_Individ_fdr[Diff_Expr_Individ_fdr$CpG == x,]
    
    Number_pos_sig = nrow(Curr_df[Curr_df$logFC > 0,])
    Number_negat_sig = nrow(Curr_df[Curr_df$logFC < 0,])
    Signif_score = (Number_pos_sig + Number_negat_sig)/length(studies)
    
    Number_pos_sig_fdr = nrow(Curr_df_fdr[Curr_df_fdr$logFC > 0,])
    Number_negat_sig_fdr = nrow(Curr_df_fdr[Curr_df_fdr$logFC < 0,])
    Signif_score_fdr = (Number_pos_sig_fdr + Number_negat_sig_fdr)/length(studies)
    
    stat_df = data.frame("CpG" = x,
                         "Nominal.pos" = Number_pos_sig,
                         "Nominal.neg" = Number_negat_sig,
                         "Nominal.score" = Signif_score,
                         "FDR.pos" = Number_pos_sig_fdr,
                         "FDR.neg" = Number_negat_sig_fdr,
                         "FDR.score" = Signif_score_fdr)
    
    rm(Curr_df)
    rm(Curr_df_fdr)
    
    return(stat_df)
  }, mc.cores = 5)
  CpG_stat_individ_meth = do.call(rbind, CpG_stat_individ_meth)
  
  CpG_stat_individ_meth$Nominal.agreement = mapply(function(x,y){
    
    max_ind = max(c(x,y))
    agr_ind = max_ind/(x+y)
    return(agr_ind)
    
  }, CpG_stat_individ_meth$Nominal.pos, CpG_stat_individ_meth$Nominal.neg)
  
  CpG_stat_individ_meth$Nominal.count = CpG_stat_individ_meth$Nominal.pos + CpG_stat_individ_meth$Nominal.neg
  CpG_stat_individ_meth$FDR.count = CpG_stat_individ_meth$FDR.pos + CpG_stat_individ_meth$FDR.neg
  CpG_stat_individ_meth = CpG_stat_individ_meth[,c( "CpG",
                                                    "Nominal.pos",
                                                    "Nominal.neg",
                                                    "Nominal.score",
                                                    "Nominal.agreement",
                                                    "Nominal.count",
                                                    
                                                    "FDR.pos",
                                                    "FDR.neg",
                                                    "FDR.score",
                                                    "FDR.count")]
  
  fwrite(CpG_stat_individ_meth, "CpG_stat_individ_meth.csv", sep = ",")
  
  make_Venn_digram_list(named_list = list(
    "Nominal signif. CpGs (combined data)" = Top_table_signif$CpG,
    "Nominal signif. CpGs (separate cohorts)" = CpG_stat_individ_meth$CpG
  ), plot_full_path = "Comparison_pooled_vs_individ")
  
  fwrite(CpG_stat_individ_meth, "CpG_stat_individ_meth.csv", sep = ",")
  
} else {
  
  CpG_stat_individ_meth = smart_fread("CpG_stat_individ_meth.csv")
  
}

################### Meta-analysis of individual cohorts ###################
# Example for cg01194782
CpG_meta = "cg01194782"
Diff_Expr_Individ_test = Diff_Expr_Individ[Diff_Expr_Individ$CpG == CpG_meta,]
Study_names = stri_replace_all_fixed(Diff_Expr_Individ_test$Study, pattern = "_", replacement = " ")
meta_model = rma.uni(yi = logFC, vi = SE^2, data = Diff_Expr_Individ_test, method = "SJ", weighted = TRUE)
summary(meta_model)

# Example Visualization
weights = fmtx(weights(meta_model), digits=1)
sav = forest(meta_model, slab = Study_names, ilab = weights)
k = nrow(Diff_Expr_Individ_test)
colp <- "red"
segments(coef(meta_model), 0, coef(meta_model), k, col=colp, lty="33", lwd=0.8)

# Add text
par(xpd=NA)
par(cex=sav$cex, font=2)
# Headers
text(sav$xlim[1], k+2.5, pos=4, "Cohort")
text(-0.58, k+2.5, pos=4, "Weight %")
text(0, k+2.7, "Log2FC,\n(95% CI)")
segments(sav$ilab.xpos[1]-0.22, k+2.8, sav$ilab.xpos[2]+0.13, k+2.8)
text(sav$xlim[2]-0.10, k+2.7, "Log2FC\n(95% CI)")

# Use a non-bold font for the rest of the text
par(cex=sav$cex, font=1)
text(sav$ilab.xpos[3], 0, "100.0")
text(sav$xlim[1], -0.5, pos=4, bquote(paste("Test for heterogeneity: ",
                                          tau^2, "=", .(fmtx(meta_model$tau2, digits=2)), "; ",
                                          chi^2, "=", .(fmtx(meta_model$QE, digits=2)),
                                          ", df=", .(meta_model$k - meta_model$p), ", ",
                                          .(fmtp(meta_model$QEp, digits=2, pname="P", add0=TRUE, equal=TRUE)), "; ",
                                          I^2, "=", .(round(meta_model$I2)), "%")))
title(paste0("CpG: ", CpG_meta))

# Calculate meta for all CpGs that were significant at least once
CALCULATE_META = FALSE

if (CALCULATE_META){
  
  Signif_CpG_to_test = unique(Diff_Expr_Individ_signif$CpG)
  model_outputs_meta = mclapply(Signif_CpG_to_test, function(x){
    meta_df = Diff_Expr_Individ[Diff_Expr_Individ$CpG == x,]
    meta_model = rma.uni(yi = logFC, vi = SE^2, data = meta_df, method = "SJ", weighted = TRUE)
    
    output_stats = data.frame(
      CpG = x,
      meta_LFc = meta_model$b,
      meta_se = meta_model$se,
      meta_pval = meta_model$pval,
      tau2 = meta_model$tau2,
      I2 = meta_model$I2,
      H2 = meta_model$H2,
      Q = meta_model$QE,
      Q.p = meta_model$QEp)
    
    rm(meta_df)
    rm(meta_model)
    return(output_stats)
    
  }, mc.cores = 10)
  model_outputs_meta_df = do.call(rbind, model_outputs_meta)
  model_outputs_meta_df$meta_FDR = p.adjust(model_outputs_meta_df$meta_pval, method = "fdr")
  write.csv(model_outputs_meta_df, "model_outputs_meta_df.csv", row.names = FALSE)
  
  model_outputs_meta_df_signif = model_outputs_meta_df[model_outputs_meta_df$meta_pval < 0.05,]
  
  # Grouping stats for individual analyses
  all(model_outputs_meta_df$CpG == CpG_stat_individ_meth$CpG)
  Grouped_stats_individual_studies = cbind(CpG_stat_individ_meth, model_outputs_meta_df[,-1])
  Grouped_stats_individual_studies = arrange(Grouped_stats_individual_studies, meta_pval)
  fwrite(Grouped_stats_individual_studies, "Grouped_stats_individual_studies.csv", sep = ",")
  
} else {
  model_outputs_meta_df = smart_fread("model_outputs_meta_df.csv")
  model_outputs_meta_df_signif = model_outputs_meta_df[model_outputs_meta_df$meta_pval < 0.05,]
  Grouped_stats_individual_studies = smart_fread("Grouped_stats_individual_studies.csv")
}

#


################### Venn diagram for overlaps ###################
overlapping_CpGs = list(
  "Pooled analysis (signif. nominal)" = Top_table_signif$CpG,
  "Individ cohorts (signif. nominal)" = unique(Diff_Expr_Individ_signif$CpG),
  "Meta individ cohorts (signif. nominal)" = model_outputs_meta_df_signif$CpG
)
make_Venn_digram_list(overlapping_CpGs, plot_full_path = "CpG_by_analysis.pdf", palette = 7)
intersecting_CpGs = Reduce(intersect, overlapping_CpGs)

ILLUM_450K_ANNOT_selected_intersect = ILLUM_450K_ANNOT_selected[ILLUM_450K_ANNOT_selected$Name %in% intersecting_CpGs,]
ILLUM_450K_ANNOT_selected_intersect_AHRR = ILLUM_450K_ANNOT_selected_intersect[stri_detect_fixed(ILLUM_450K_ANNOT_selected_intersect$Updated_gene_name, 
                                                                                                 pattern = "AHRR"),] # 0 CpGs

################### GO enrichment analysis ###################
enrichment_FINAL = missMethyl::gometh(sig.cpg = intersecting_CpGs, all.cpg = Top_table$CpG)
enrichment_FINAL = enrichment_FINAL[enrichment_FINAL$ONTOLOGY == "BP", ]
enrichment_FINAL = arrange(enrichment_FINAL, P.DE)
write.csv(enrichment_FINAL, "GO_BP_enrichment_overlap.csv", row.names = FALSE)


################### Regulatory chromatin enrichment ###################
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4530010/ Human 111 reference epigenomes
# https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#core_15state
# E073 Pre-frontal cortex, E062 Primary mononuclear cells (blood)
ILLUM_450K_ANNOT_selected_chrom_state = ILLUM_450K_ANNOT_selected
ENCODE_15st_cortex = fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/HGNC/15_state_model/E073_15_coreMarks_mnemonics.bed")
ENCODE_15st_PBMC = fread("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/HGNC/15_state_model/E062_15_coreMarks_mnemonics.bed")

indeces = 1:nrow(ILLUM_450K_ANNOT_selected_chrom_state)
annot_chrom_states = mclapply(indeces, function(i){
  
  cpg = ILLUM_450K_ANNOT_selected_chrom_state$Name[i]
  chrom = ILLUM_450K_ANNOT_selected_chrom_state$chr[i]
  pos = ILLUM_450K_ANNOT_selected_chrom_state$pos[i]
  
  # select target sets
  target_blood = ENCODE_15st_PBMC[ENCODE_15st_PBMC$V1 == chrom,]
  target_cortex = ENCODE_15st_cortex[ENCODE_15st_cortex$V1 == chrom,]
  
  # select specified data
  target_blood = target_blood[target_blood$V2 <= pos & target_blood$V3 >= pos,]
  target_cortex = target_cortex[target_cortex$V2 <= pos & target_cortex$V3 >= pos,]
  
  # selecting unique or the first value
  target_blood = unique(target_blood$V4)
  target_cortex = unique(target_cortex$V4)
  target_blood = target_blood[1]
  target_cortex = target_cortex[1]
  
  df = data.frame(cpg = cpg, St15_blood = target_blood, St15_cortex = target_cortex)
  
  return(df)
  
}, mc.cores = 10)
annot_chrom_states = list_to_df(annot_chrom_states)
all(ILLUM_450K_ANNOT_selected_chrom_state$Name == annot_chrom_states$cpg) # Order is matching
ILLUM_450K_ANNOT_selected_chrom_state$St15_blood = annot_chrom_states$St15_blood
ILLUM_450K_ANNOT_selected_chrom_state$St15_cortex = annot_chrom_states$St15_cortex
write.csv(ILLUM_450K_ANNOT_selected_chrom_state, "ILLUM_450K_ANNOT_selected_chrom_state.csv", row.names = FALSE)


# Enrichemnt
ChromToCpG_blood = ILLUM_450K_ANNOT_selected_chrom_state[,c("St15_blood", "Name")]
ChromToCpG_cortex = ILLUM_450K_ANNOT_selected_chrom_state[,c("St15_cortex", "Name")]
CpG_set = as.character(intersecting_CpGs)
CpG_universe = ILLUM_450K_ANNOT_selected_chrom_state$Name

# Up-methylated CpGs and Down-methylated CpGs (in cases)
CpG_set_up = CpG_set[CpG_set %in% model_outputs_meta_df_signif[model_outputs_meta_df_signif$meta_LFc > 0, "CpG"]]
CpG_set_down = CpG_set[CpG_set %in% model_outputs_meta_df_signif[model_outputs_meta_df_signif$meta_LFc < 0, "CpG"]]

# Up
Chrom_blood_enrichment = enricher(gene = CpG_set_up, 
                               TERM2GENE = ChromToCpG_blood,
                               minGSSize = 1, 
                               maxGSSize = 1e7)
Chrom_blood_enrichment_result = Chrom_blood_enrichment@result
Chrom_blood_enrichment_result$tissue = "Primary mononuclear cells from peripheral blood"
Chrom_blood_enrichment_result$change = "Increased methylation in depressed"

Chrom_cortex_enrichment = enricher(gene = CpG_set_up, 
                                  TERM2GENE = ChromToCpG_cortex,
                                  minGSSize = 1, 
                                  maxGSSize = 1e7)
Chrom_cortex_enrichment_result = Chrom_cortex_enrichment@result
Chrom_cortex_enrichment_result$tissue = "Brain Dorsolateral Prefrontal Cortex"
Chrom_cortex_enrichment_result$change = "Increased methylation in depressed"

# Down
Chrom_blood_enrichment_down = enricher(gene = CpG_set_down, 
                                  TERM2GENE = ChromToCpG_blood,
                                  minGSSize = 1, 
                                  maxGSSize = 1e7)
Chrom_blood_enrichment_result_down = Chrom_blood_enrichment_down@result
Chrom_blood_enrichment_result_down$tissue = "Primary mononuclear cells from peripheral blood"
Chrom_blood_enrichment_result_down$change = "Decreased methylation in depressed"

Chrom_cortex_enrichment_down = enricher(gene = CpG_set_down, 
                                   TERM2GENE = ChromToCpG_cortex,
                                   minGSSize = 1, 
                                   maxGSSize = 1e7)
Chrom_cortex_enrichment_result_down = Chrom_cortex_enrichment_down@result
Chrom_cortex_enrichment_result_down$tissue = "Brain Dorsolateral Prefrontal Cortex"
Chrom_cortex_enrichment_result_down$change = "Decreased methylation in depressed"

# Description for Chromatin states from https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#core_15state
content = read_html("https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#core_15state")
tables = content %>% html_table(fill = TRUE)
description_15st = tables[[1]]
description_15st$MNEMONIC = paste0(description_15st$`STATE NO.`, "_", description_15st$MNEMONIC)
# prepare enrichment results 
chrom_enrich_list = list(
  Chrom_blood_enrichment_result,
  Chrom_cortex_enrichment_result,
  Chrom_blood_enrichment_result_down,
  Chrom_cortex_enrichment_result_down
)
chrom_enrich_list = lapply(chrom_enrich_list, function(x){
  
  colnames(x) = c("State",
                  "Description",
                  "SetRatio",
                  "BgRatio",
                  "pvalue",
                  "p.adjust",
                  "qvalue",
                  "CpGs",
                  "Count",
                  "tissue",
                  "change")
  
  x$Description = sapply(x$Description, function(z){
    z = description_15st[description_15st$MNEMONIC == z, "DESCRIPTION"]
    return(z)
  })
  return(x)
})
chrom_enrich_list = list_to_df(chrom_enrich_list)


################### eQTM analysis ###################
eQTM_blood = smart_fread("BIOS_QTL_data/2015_09_02_cis_eQTMsFDR0.05-CpGLevel.txt")
eQTMs_selected = eQTM_blood[eQTM_blood$SNPName %in% intersecting_CpGs, ]
regulated_genes = eQTMs_selected$HGNCName


freq_gene = as.data.frame(table(regulated_genes))
freq_gene = arrange(freq_gene, -Freq)
regulated_genes_unique = unique(regulated_genes)
freq_gene$CpG = sapply(freq_gene$regulated_genes, function(x){
  CpGs = eQTMs_selected[eQTMs_selected$HGNCName == x,]
  CpGs = CpGs$SNPName
  CpGs = paste0(CpGs, collapse = ";")
  return(CpGs)
})

write.csv(freq_gene, "freq_gene.csv", row.names = FALSE)

# Modifying UCSC_RefGene_Name
eQTL_anno = ILLUM_450K_ANNOT_selected
eQTL_anno$UCSC_RefGene_Name = sapply(eQTL_anno$Name, function(x){
  
  if (x %!in% eQTM_blood$SNPName){
    return("")
  }
  
  cpg_eqtm = eQTM_blood[eQTM_blood$SNPName == x,]
  genes = cpg_eqtm$HGNCName
  genes = paste0(genes, collapse = ";")
  return(genes)
  
})

eQTL_anno$UCSC_RefGene_Group = sapply(eQTL_anno$UCSC_RefGene_Name, function(x){
  
  if (x==""){
    return("NA")
  }
  
  x = stri_split_fixed(x, pattern = ";")
  x = unlist(x)
  
  attachment = rep("NA", times = length(x))
  attachment = paste0(attachment, collapse = ";")
  
  return(attachment)
})

eQTL_anno$UCSC_RefGene_Accession = "NA"

enrichment_EQTM = missMethyl::gometh(sig.cpg = intersecting_CpGs, 
                                     all.cpg = Top_table$CpG,
                                     prior.prob = FALSE,
                                     anno = eQTL_anno, 
                                     fract.counts = FALSE, 
                                     sig.genes = FALSE)

enrichment_EQTM = enrichment_EQTM[enrichment_EQTM$ONTOLOGY == "BP", ]
enrichment_EQTM = arrange(enrichment_EQTM, P.DE)
write.csv(enrichment_FINAL, "GO_BP_enrichment_overlap_eQTM.csv", row.names = FALSE)

################### Saving biol. analysis ###################
# combining files
results_biology = list("Enrichment_BP_init" = enrichment_FINAL,
                       "Gene_CpG_eQTM" = eQTMs_selected,
                       "Gene_CpG_eQTM_freq" = freq_gene,
                       "Enrichment_BP_eQTM" = enrichment_EQTM,
                       "Enrichment_chrom_states" = chrom_enrich_list)
openxlsx::write.xlsx(results_biology, "results_biology.xlsx", overwrite = TRUE)

################### Generate phenotypes table for cohorts ###################
studies = unique(Combined_pheno_no_NA$Study)
table(Combined_pheno_no_NA$Study, Combined_pheno_no_NA$Depression)

charact_list = lapply(studies, function(x){
  
  tmp_df = Combined_pheno_no_NA[Combined_pheno_no_NA$Study == x,]
  charact_df = characterize_dataset_generelized_two_subgroups(dataset = tmp_df,
                                                              study_char = x,
                                                              contrast_col_number = 4,
                                                              contrast_vector = c("Control", "Case"),
                                                              participants_col_number = 1,
                                                              model_covariates_col_vector = c(2,3),
                                                              columns_to_characterise_vector = c(4, 2, 3),
                                                              Remove_NA_predictors = TRUE,
                                                              drop_P = TRUE,
                                                              simplif_P = 3)
  return(charact_df)
})

names(studies) = c("PSY (screen)", 
                   "PSY (recall)", 
                   "GSE125105 (MPIP1)", 
                   "GSE72680 (GRADY)", 
                   "GSE113725 (RDE)",
                   "GSE198904 (DHRC)",
                   "GSE198904 (OBS)",
                   "GSE74414 (MPIP2)")

charact_list = lapply(charact_list, function(study){
  study$Table$Name = as.character(study$Table$Name)
  study$Table$Name[1] = names(studies[studies == study$Table$Name[1]])
  study$Table$Name[2] = stri_replace_all_fixed(str = study$Table$Name[2], pattern = "Initial dataset includes ", replacement = "")
  study$Table$Name[4] = stri_replace_all_fixed(str = study$Table$Name[4], pattern = "Depression", replacement = "Depression status")
  
  study$Table =  study$Table[-3,]
  
  colnames(study$Table) = c("","Controls", "Cases")
  
  return(study$Table)
})

charact_df_total = characterize_dataset_generelized_two_subgroups(dataset = Combined_pheno_no_NA,
                                                                  study_char = "Total",
                                                                  contrast_col_number = 4,
                                                                  contrast_vector = c("Control", "Case"),
                                                                  participants_col_number = 1,
                                                                  model_covariates_col_vector = c(2,3),
                                                                  columns_to_characterise_vector = c(4, 2, 3),
                                                                  Remove_NA_predictors = TRUE,
                                                                  drop_P = TRUE,
                                                                  simplif_P = 3)
charact_df_total = charact_df_total$Table
charact_df_total$Name = as.character(charact_df_total$Name)
charact_df_total$Name[2] = stri_replace_all_fixed(str = charact_df_total$Name[2], pattern = "Initial dataset includes ", replacement = "")
charact_df_total$Name[4] = stri_replace_all_fixed(str = charact_df_total$Name[4], pattern = "Depression", replacement = "Depression status")
charact_df_total =  charact_df_total[-3,]
colnames(charact_df_total) = c("","Controls", "Cases")

charact_list = do.call(rbind, charact_list)
charact_list = rbind(charact_list, charact_df_total)
openxlsx::write.xlsx(charact_list, "Demographics_RAW.xlsx", overwrite = TRUE)


################### Viewing stats ###################

# Unique CpGs in the individual cohort
length(unique(Diff_Expr_Individ_signif$CpG)) # 158618

# Total FDR significant associations
nrow(Diff_Expr_Individ_fdr) # 4422

# Unique FDR significant CpGs
length(unique(Diff_Expr_Individ_fdr$CpG)) # 4418

# Number of CpGs with high significant agreements
tmp_dat = Grouped_stats_individual_studies[Grouped_stats_individual_studies$Nominal.agreement == 1,]
tmp_dat  = tmp_dat[tmp_dat$Nominal.count > 3, ]

# Number of meta-significant CpGs
tmp_dat = Grouped_stats_individual_studies[Grouped_stats_individual_studies$meta_pval < 0.05,]
table(tmp_dat$Nominal.count)

tmp_dat = tmp_dat[tmp_dat$FDR.count > 0,] # 29

# FDR chisq
freq_tab = data.frame(Total = c(304765-4418, 2451-29),SigFDR = c(4418, 29))
chisq.test(freq_tab, correct = FALSE)
# 	Pearson's Chi-squared test with Yates' continuity correction
# X-squared = 1.2101, df = 1, p-value = 0.2713

# Number of participants per subset
table(CV_pheno$Study, CV_pheno$Depression)
table(CV_pheno$Depression)
table(test_pheno$Study, test_pheno$Depression)
table(test_pheno$Depression)

# Calculate AUC gain
harm_data = openxlsx::read.xlsx("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Model_performances.xlsx", sheet = 1)
non_harm_data = openxlsx::read.xlsx("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Model_performances_rawdata.xlsx", sheet = 1)
auc_gain_CV = non_harm_data$`10.CV.AUC.(hold-out)` -   harm_data$`10.CV.AUC.(hold-out)`
mean(auc_gain_CV) # 0.1469167

harm_data_test = openxlsx::read.xlsx("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Model_performances.xlsx", sheet = 2)
non_harm_data_test = openxlsx::read.xlsx("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/Model_performances_rawdata.xlsx", sheet = 2)
harm_data_test$Model == non_harm_data_test$Model
auc_gain_test = non_harm_data_test$`AUC.(last.hold-out)`-   harm_data_test$`AUC.(last.hold-out)`
mean(auc_gain_test) # 0.1680625

################### Sensitivity analysis for meta with smoking ###################

# Generate beta values and pheno for smoking estimation
phenos_smoking = Combined_pheno_no_NA
phenos_smoking$sex = ifelse(phenos_smoking$Sex == "Male", 1,2)
phenos_smoking$Sex = NULL

library(EpiSmokEr)
smoking_datasets = list()

# Calculating smoking score using Elliott et al. 2014
# Data has from 158 to 172 of 187 CpGs, initially reported by Zellinger et al. 2013

mapping = c(
  "PSY_SCR_mval",
  "PSY_RC_mval",
  "GSE125105_MPIP_mval",
  "GSE72680_GRADY_mval",
  "GSE113725_RDE_mval",
  "GSE198904_DHRC_mval",
  "GSE198904_OBS_mval",
  "GSE74414_MPIP2_mval"
)
names(mapping) = unique(as.character(phenos_smoking$Study))

for (i in 1:length(unique(phenos_smoking$Study))){
  cohort = unique(as.character(phenos_smoking$Study))[i]
  print(cohort)
  curr_df = phenos_smoking[phenos_smoking$Study == cohort, ]
  
  # Selecting methylation dataset and beta values
  mvals = get(mapping[cohort])
  mvals = mvals[,curr_df$ID]
  print("Mapping correct:")
  print(base::all(colnames(mvals) == curr_df$ID))
  curr_beta = m2beta(mvals)
  
  rownames(curr_df) = curr_df$ID
  smoke_status = epismoker(dataset=curr_beta, samplesheet = curr_df, method = "SSc")
  smoking_datasets[[i]] = smoke_status
  
}
smoking_datasets = list_to_df(smoking_datasets)
base::all(smoking_datasets$SampleName == phenos_smoking$ID) # EpiSmokEr overwrites method all!?

detach("package:EpiSmokEr", unload=TRUE)

all(smoking_datasets$SampleName == phenos_smoking$ID) # Now all correctly works

phenos_smoking$smoking_score = smoking_datasets$smokingScore
# treshold 17.55 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3915234/
phenos_smoking$above_treshold = ifelse(phenos_smoking$smoking_score > 17.55, "Above", "Below")
table(phenos_smoking$above_treshold, phenos_smoking$Study) # Potentially only 4 participants GSE125105_MPIP are smokers based on European Threshold

# treshold 11.79 for Asian population
phenos_smoking$above_treshold_strict = ifelse(phenos_smoking$smoking_score > 11.79, "Above", "Below")
table(phenos_smoking$above_treshold_strict) # Only 44 participants are above strict threshold
44/1942 #0.02265705 Potentially only 2% of the sample are smokers

plot = ggplot(data = phenos_smoking, aes(x = Study, y = smoking_score)) +
  geom_boxplot() +
  geom_hline(yintercept = 17.55, col="red") +
  geom_hline(yintercept = 11.79, col="blue") +
  labs(y = "Estimated smoking score (Elliott H. 2014)", x = "Cohort-batch") +
  ggtitle("Estimated smoking score vs cohort-batch (cohort-level preprocessing)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
        strip.text.y = element_text(angle = 45, size = 8),
        strip.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_line(linewidth = 1, color = "black"),
        panel.grid = element_line(linewidth = 0.25, color = "black", linetype = 2),
        axis.text.y = element_text(face = "bold"))

plot
pdf(file = "Fig S5.pdf", width = 10,height = 8)
plot
dev.off()





