import os
import numpy
import subprocess
import datatable
import pandas

from sklearn.feature_selection import VarianceThreshold
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_classif
from sklearn.feature_selection import SelectFromModel
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import ExtraTreesClassifier



################### Defining feature selectors ###################

# Limma-based feature selector
def select_features_R_limma(train_methyl,
                            train_pheno,
                            CpG_names,
                            current_fold,
                            base_dir,
                            command = "Rscript",
                            path2script ="/HelperScriptsPY/RunLimmaMatrix.R",
                            feature_number = 50):

  # Writing Limma inputs                      
  pheno_path = base_dir + "/pheno_" + current_fold + ".csv"
  meth_path = base_dir + "/meth_" + current_fold + ".csv"
  toptable_path = base_dir + "/TopTable_" + current_fold + ".csv"

  # Saving things 
  train_pheno.to_csv(pheno_path)
  train_methyl = datatable.Frame(train_methyl)
  train_methyl.names = CpG_names
  train_methyl.to_csv(meth_path)

  # Running Limma inside R
  args = [base_dir, pheno_path, meth_path, toptable_path]
  cmd = [command, path2script] + args
  subprocess.run(cmd)

  # Removing junk
  os.remove(pheno_path)
  os.remove(meth_path)

  # Selecting features
  toptable = pandas.read_csv(toptable_path)
  top_features = toptable["CpG"].values[range(feature_number)]
     

  selection_idx = train_methyl.names
  selection_idx = [True if x in top_features else False for x in selection_idx]
  return selection_idx
  

# Select Top % hypervar CpGs
# Base function
def selector_variance(train_methyl,
                      train_pheno,
                      CpG_names,
                      current_fold,
                      base_dir,
                      percentile):
    
    current_fold = str(current_fold)
    print("Performing variance-based feature selection")
    print(f"Using percentile: {percentile}")

    variance_methyl = numpy.var(train_methyl, axis=0)
    percentile_treshold = numpy.percentile(variance_methyl, percentile)
    selector = VarianceThreshold(threshold=percentile_treshold)
    selected_data = selector.fit_transform(train_methyl)
    selection_indeces = selector.get_support()
    print("Selected features")
    print(numpy.unique(selection_indeces, return_counts=True))

    # saving selected CpGs
    CpG_names = numpy.array(CpG_names, dtype="str")
    selected_names = CpG_names[selection_indeces]
    path_CpGs = base_dir + "/selected_CpGs_" + current_fold + ".txt"
    file = open(path_CpGs,'w')
    for item in selected_names:
        file.write(item+"\n")
    file.close()

    return selection_indeces

# Top 5%
def selector_variance_5(train_methyl,
                        train_pheno,
                        CpG_names,
                        current_fold,
                        base_dir):
    
    indeces = selector_variance(train_methyl=train_methyl,
                                train_pheno=train_pheno,
                                CpG_names=CpG_names,
                                current_fold=current_fold,
                                base_dir=base_dir,
                                percentile=95)
    return indeces

# Top 1%
def selector_variance_1(train_methyl,
                        train_pheno,
                        CpG_names,
                        current_fold,
                        base_dir):
    
    indeces = selector_variance(train_methyl=train_methyl,
                                train_pheno=train_pheno,
                                CpG_names=CpG_names,
                                current_fold=current_fold,
                                base_dir=base_dir,
                                percentile=99)
    return indeces

# Top 0.1%
def selector_variance_01(train_methyl,
                        train_pheno,
                        CpG_names,
                        current_fold,
                        base_dir):
    
    indeces = selector_variance(train_methyl=train_methyl,
                                train_pheno=train_pheno,
                                CpG_names=CpG_names,
                                current_fold=current_fold,
                                base_dir=base_dir,
                                percentile=99.9)
    return indeces


# Select k-best features based on ANOVA F score
def selector_k_best_f_200(train_methyl,
                          train_pheno,
                          CpG_names,
                          current_fold,
                          base_dir):
    current_fold = str(current_fold)
    print("Performing k-best ANOVA F feature selection")

    selector = SelectKBest(f_classif, k=200)
    selected_data_F = selector.fit_transform(train_methyl, train_pheno)
    selection_indeces = selector.get_support()
    print("Selected features")
    print(numpy.unique(selection_indeces, return_counts=True))

    # saving selected CpGs
    CpG_names = numpy.array(CpG_names, dtype="str")
    selected_names = CpG_names[selection_indeces]
    path_CpGs = base_dir + "/selected_CpGs_" + current_fold + ".txt"
    file = open(path_CpGs,'w')
    for item in selected_names:
        file.write(item+"\n")
    file.close()

    return selection_indeces


# Select features based on L1-penalized models (Linear C-support vector classifier)
def selector_lin_svc_200(train_methyl,
                         train_pheno,
                         CpG_names,
                         current_fold,
                         base_dir):
    current_fold = str(current_fold)
    print("Performing L1 LSVC-based feature selection")
    print(train_methyl.shape)
    print(train_pheno.shape)
    
    lsvc = LinearSVC(C=1, penalty="l1", dual="auto", max_iter=5000).fit(train_methyl, train_pheno)
    selector = SelectFromModel(lsvc, prefit=True, max_features=200)
    selected_data_LSVC = selector.transform(train_methyl)
    selection_indeces = selector.get_support()
    print("Selected features")
    print(numpy.unique(selection_indeces, return_counts=True))

    # saving selected CpGs
    CpG_names = numpy.array(CpG_names, dtype="str")
    selected_names = CpG_names[selection_indeces]
    path_CpGs = base_dir + "/selected_CpGs_" + current_fold + ".txt"
    file = open(path_CpGs,'w')
    for item in selected_names:
        file.write(item+"\n")
    file.close()

    return selection_indeces


# Select features based on L1-penalized models (logistic regression)
def selector_log_reg_200(train_methyl,
                         train_pheno,
                         CpG_names,
                         current_fold,
                         base_dir):
    current_fold = str(current_fold)
    print("Performing L1 LogRegr-based feature selection")
    
    logregr = LogisticRegression(C=1, penalty="l1", solver="liblinear", max_iter=5000).fit(train_methyl, train_pheno)
    selector = SelectFromModel(logregr, prefit=True, max_features=200)
    selected_data_regr = selector.transform(train_methyl)
    selection_indeces = selector.get_support()
    print("Selected features")
    print(numpy.unique(selection_indeces, return_counts=True))

    # saving selected CpGs
    CpG_names = numpy.array(CpG_names, dtype="str")
    selected_names = CpG_names[selection_indeces]
    path_CpGs = base_dir + "/selected_CpGs_" + current_fold + ".txt"
    file = open(path_CpGs,'w')
    for item in selected_names:
        file.write(item+"\n")
    file.close()

    return selection_indeces


# Select features based on decision trees
def selector_trees_200(train_methyl,
                       train_pheno,
                       CpG_names,
                       current_fold,
                       base_dir):
    current_fold = str(current_fold)
    print("Performing Tree-based feature selection")
    
    trees = ExtraTreesClassifier().fit(train_methyl, train_pheno)
    selector = SelectFromModel(trees, prefit=True, max_features=200)
    selected_data_trees = selector.transform(train_methyl)
    selection_indeces = selector.get_support()
    print("Selected features")
    print(numpy.unique(selection_indeces, return_counts=True))

    # saving selected CpGs
    CpG_names = numpy.array(CpG_names, dtype="str")
    selected_names = CpG_names[selection_indeces]
    path_CpGs = base_dir + "/selected_CpGs_" + current_fold + ".txt"
    file = open(path_CpGs,'w')
    for item in selected_names:
        file.write(item+"\n")
    file.close()

    return selection_indeces