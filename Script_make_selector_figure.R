setwd("/home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta")


# Setting options
getOption("scipen") # Default number notation is 0
options(scipen=999)
options(stringsAsFactors = FALSE)

################### Package import ###################
# Importing packages (Note: Not all of them may be required)
library(dplyr)
library(ggplot2)
library(data.table)

################### Defining functions ###################

# NOT IN operator
'%!in%' = function(x,y){!('%in%'(x,y))}


################### Data import ###################
dataset = read.csv("merged_fold_data_10CV.csv")
dataset$X = NULL
datatset_2 = read.csv("merged_fold_data_10CV_Limma.csv")
datatset_2$X = NULL
dataset = rbind(dataset, datatset_2)
str(dataset)
dataset_cohort = dataset[dataset$Data.harmonization == "batch-level",]

dataset_cohort$Feature.selection = factor(dataset_cohort$Feature.selection, 
                                          levels = c("Top 200 Limma",
                                                     "Top 5% Hypervar. CpGs",
                                                     "Top 1% Hypervar. CpGs",
                                                     "Top 0.1% Hypervar. CpGs",
                                                     "Best ANOVA F (200 CpGs)",
                                                     "L1 Lin. SVC (200 CpGs)",
                                                     "L1 Logist. regr. (200 CpGs)",
                                                     "ExtraTrees (200 CpGs)"))

plot = ggplot(data = dataset_cohort, aes(x = Feature.selection, y = X10.CV.AUC..hold.out., fill = Feature.selection)) +
  scale_fill_brewer(palette="Dark2") +
  geom_boxplot(notchwidth = 0.1, outlier.size = 0.5, fatten=0.75) +
  facet_grid(rows = vars(Model)) +
  labs(y = "Fold-level CV AUC (hould-out, no harmonization)x30", x = "Feature selection", fill = "Feature selection") +
  ggtitle("Feature selection strategies in 10x3-fold CV on non-harmonized data") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
        strip.text.y = element_text(angle = 45, size = 8),
        strip.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_line(linewidth = 1, color = "black"),
        panel.grid = element_line(linewidth = 0.25, color = "black", linetype = 2),
        axis.text.y = element_text(face = "bold"))
plot

pdf(file = "Fig S6.pdf", width = 7, height = 20)
plot
dev.off()
