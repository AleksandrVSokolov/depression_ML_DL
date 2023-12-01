# Importing modules
import sys
import os
sys.path.append(os.getcwd())

import pandas
import os
import gc
from HelperScriptsPY import helper_functions
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score, roc_auc_score, auc, roc_curve
from sklearn.utils import shuffle
from sklearn import svm
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.linear_model import LogisticRegression


"""
# Classifiers:     
    LogisticRegr,
    SVM_c,
    SVM_c_RBF,
    DecisionTree,
    RandomForest,
    AdaBoost,
    Elastic_net_classifier,
    RidgeLogisticRegr,
    Lasso_classifier


    
Example commands:

conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta
/home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_ML.py  "RandomForest" "test_folder_ML" 3 \
    "ALL" "Mval" "Optimization" "all_cohorts" 20 "Main_analysis/ML_data/"

conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta
/home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_ML.py  "LogisticRegr" "test_folder_ML_RAW" 3 \
    "ALL" "Mval" "FINAL_Test" "all_cohorts" 200 "Main_analysis/ML_data_raw/"

conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta
/home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_ML.py  "RandomForest" "test_folder_ML" 3 \
    "Top_200_all_consistent_CpGs.txt" "Mval" "Optimization" "all_cohorts" 200 "Main_analysis/ML_data_/"

"""

#### Parsing arguments ####
ARG_LIST = sys.argv
MODEL_ARG = str(ARG_LIST[1])
FOLDER_ARG = str(ARG_LIST[2])
FOLD_ARG = int(ARG_LIST[3])
CONSIST_ARG = str(ARG_LIST[4])
VALUES_ARG = str(ARG_LIST[5])
FINAL_TEST_ARG = str(ARG_LIST[6])
COHORT_POOL_ARG = str(ARG_LIST[7])
LIMMA_FEATURE_ARG = int(ARG_LIST[8])
DATA_ARG = str(ARG_LIST[9])

print('------------------------------------------------------------------------')
print(f'Training ML model {MODEL_ARG}')
print(f'Training for n-folds {FOLD_ARG}')
print(f'Selected CpGs: {CONSIST_ARG}')
print(f'Selected values: {VALUES_ARG}')
print(f'Limma selected features: {LIMMA_FEATURE_ARG}')
print(f'Perform final test: {FINAL_TEST_ARG}')

cons_str = "_cpg_" + CONSIST_ARG
cons_str = cons_str.replace(".txt","")

# Preparing target directory
dir_path = (FOLDER_ARG + "/" + MODEL_ARG + "_fo_" + str(FOLD_ARG) + "_va_" + 
            VALUES_ARG + cons_str + "_limma_" + str(LIMMA_FEATURE_ARG) + "___" + FINAL_TEST_ARG)

# Preapring folder
os.makedirs(dir_path, exist_ok=True)



#################################################################################
##########################    Section    ########################################
#################################################################################

if COHORT_POOL_ARG == "all_cohorts":
    print("Using all cohorts")
    cohort_pool=["GSE113725_RDE",
    "GSE125105_MPIP",
    "GSE198904_DHRC",
    "GSE198904_OBS",
    "GSE72680_GRADY",
    "GSE74414_MPIP2",
    "PSY_RC",
    "PSY_SCR"]
else:
    print("Using specific cohort")
    cohort_pool = [COHORT_POOL_ARG]

methylation_data = helper_functions.load_DNA_methylation_data(select_cpg=CONSIST_ARG, 
                                                              cohort_pool=cohort_pool, 
                                                              cpg_folder=DATA_ARG)
CpG_names = methylation_data["CpG_names"]


if FINAL_TEST_ARG == "Optimization":

    print("Script is run in the model optimization mode")

    # Selecting values
    if VALUES_ARG == "Mval":
        print("Selecting M-values")
        meth_tensor = methylation_data["CV_mval_tensor"]
    elif VALUES_ARG == "Betaval":
        print("Selecting Beta-values")
        meth_tensor = methylation_data["CV_beta_tensor"] 
    else:
        raise ValueError("Wrong value parameter")
    
    depr_status_tensor = methylation_data["depr_status_tensor_CV"]

    # Print tensor shapes
    print("CV methylation shape:")
    print(meth_tensor.shape)
    print("CV depr shape:")
    print(depr_status_tensor.shape)

else:

    print("Script is run in the final test mode")

    # Selecting values
    if VALUES_ARG == "Mval":
        print("Selecting M-values")
        meth_tensor_train = methylation_data["CV_mval_tensor"]
        meth_tensor_test = methylation_data["test_mval_tensor"]
    elif VALUES_ARG == "Betaval":
        print("Selecting Beta-values")
        meth_tensor_train = methylation_data["CV_beta_tensor"]
        meth_tensor_test = methylation_data["test_beta_tensor"]
    else:
        raise ValueError("Wrong value parameter")
    
    depr_status_tensor_train = methylation_data["depr_status_tensor_CV"]
    depr_status_tensor_test = methylation_data["depr_status_tensor_test"]

    # Print tensor shapes
    print("Train methylation shape:")
    print(meth_tensor_train.shape)
    print("Test methylation shape:")
    print(meth_tensor_test.shape)
    print("Train depr shape:")
    print(depr_status_tensor_train.shape)
    print("Test depr shape:")
    print(depr_status_tensor_test.shape)


#################################################################################
##########################    Section    ########################################
#################################################################################
# Preparing ML Classifiers
print("Preparing classifiers")

# Regression-based classification
LogisticRegr = LogisticRegression(penalty=None,
                                  max_iter=5000)

# SVM-based classification
SVM_c = svm.SVC(kernel="linear", probability=True)
SVM_c_RBF = svm.SVC(kernel="rbf", probability=True)

# Decision-tree-based classification
DecisionTree = DecisionTreeClassifier()

# Random-forest and related
RandomForest = RandomForestClassifier()
AdaBoost = AdaBoostClassifier()

# Elastic net classifier
Elastic_net_classifier = LogisticRegression(penalty='elasticnet', 
                                            l1_ratio=0.5,
                                            max_iter=5000,
                                            solver='saga')

# Ridge regression classifier
RidgeLogisticRegr = LogisticRegression(penalty='l2',
                                       max_iter=5000)

# Lasso classifier
Lasso_classifier = LogisticRegression(penalty='l1',
                                      max_iter=5000,
                                      solver='saga')

# Models to test
model_list_ML = [
    LogisticRegr,
    SVM_c,
    SVM_c_RBF,
    DecisionTree,
    RandomForest,
    AdaBoost,
    Elastic_net_classifier,
    RidgeLogisticRegr,
    Lasso_classifier
]

model_names_ML = [
    "LogisticRegr",
    "SVM_c",
    "SVM_c_RBF",
    "DecisionTree",
    "RandomForest",
    "AdaBoost",
    "Elastic_net_classifier",
    "RidgeLogisticRegr",
    "Lasso_classifier"
]


# Classifier dict to send to a function
classifier_dict = {}
for i in range(len(model_names_ML)):
    classifier_dict[model_names_ML[i]] = model_list_ML[i]


#################################################################################
##########################    Section    ########################################
#################################################################################

def run_ML_model(X_train,
                 Y_train,
                 X_test,
                 Y_test,
                 current_fold,
                 curr_depr_pheno_train,
                 dir_path=dir_path,
                 classifier_dict=classifier_dict):

    # Convert fold to char
    current_fold = str(current_fold)

    # Running feature selection
    if CONSIST_ARG == "ALL":
        print("All CpGs were imported -> need feature selection")
        base_wd = os.getcwd()
        full_path_for_R = base_wd + "/" + dir_path
        full_path_R_script = base_wd + "/HelperScriptsPY/RunLimmaMatrix.R"
        selection_logical = helper_functions.select_features_R_limma(
            train_methyl = X_train,
            train_pheno = curr_depr_pheno_train,
            CpG_names = CpG_names,
            current_fold = current_fold,
            base_dir = full_path_for_R,
            path2script = full_path_R_script,
            feature_number = LIMMA_FEATURE_ARG)
        X_train = X_train[:,selection_logical]
        X_test = X_test[:,selection_logical]

        print(f"Selected training shape: {X_train.shape}")
        print(f"Selected testing shape: {X_test.shape}")

    elif CONSIST_ARG == "ALL_keep":
        print("All CpGs were imported but no selection requested")

    else:
        print("CpGs were provided from a list -> no selection with limma")


    # Preparing a model
    model_input = classifier_dict[MODEL_ARG]

    # Fitting model
    model_input.fit(X_train, Y_train)

    #### Evaluating model ####
    # train data
    # accuracy train
    classes_train_rounded = model_input.predict(X_train)
    classes_train = model_input.predict_proba(X_train)[:, 1]
    accuracy_train = accuracy_score(Y_train, classes_train_rounded)
    
    # Auc train
    auc_train = roc_auc_score(Y_train, classes_train)

    # saving figures
    fig_path_train = dir_path + "/" + "fold_" + current_fold + "_conf_matrix_train.pdf"

    helper_functions.make_confusion_matrix(y_true=Y_train,
                                    y_pred=classes_train_rounded,
                                    savefig=True,
                                    figname=fig_path_train)

    # cleaning
    gc.collect()

    # test data
    classes_test = model_input.predict_proba(X_test)[:, 1]
    classes_test_rounded = model_input.predict(X_test)
    accuracy_test = accuracy_score(Y_test, classes_test_rounded)

    # Auc test
    auc_test = roc_auc_score(Y_test, classes_test)

    # saving figures
    fig_path_test  = dir_path + "/" + "fold_" + current_fold + "_conf_matrix_test.pdf"

    helper_functions.make_confusion_matrix(y_true=Y_test,
                                        y_pred=classes_test_rounded,
                                        savefig=True,
                                        figname=fig_path_test)
    
    # Making ROC curve plot for test
    fpr_test, tpr_test, thresholds_test = roc_curve(Y_test, classes_test)
    roc_auc_plot_test = auc(fpr_test, tpr_test)

    fpr_test_pd = pandas.DataFrame(fpr_test)
    fpr_test_pd.to_csv(dir_path + "/" + current_fold + "_" + "fpr_test.csv")

    tpr_test_pd = pandas.DataFrame(tpr_test)
    tpr_test_pd.to_csv(dir_path + "/" + current_fold + "_" + "tpr_test.csv")

    # Making Plot
    plot_roc_test_path = dir_path + "/" + "fold_" + current_fold + "_roc_curve_test.pdf"
    helper_functions.plot_roc_curve(fpr=fpr_test,
                                    tpr=tpr_test,
                                    auc=roc_auc_plot_test,
                                    path=plot_roc_test_path)

    # cleaning
    gc.collect()

    # saving final evaluation
    final_performance = {
        "model": MODEL_ARG,
        "fold": current_fold,
        "total_folds": FOLD_ARG,
        "mval": VALUES_ARG,
        "cpg_type":CONSIST_ARG,
        "limma_features":LIMMA_FEATURE_ARG,
        "accuracy_train": accuracy_train,
        "accuracy_test": accuracy_test,
        "auc_train": auc_train,
        "auc_test": auc_test,
        "roc_auc_plot_test":roc_auc_plot_test
    }

    final_performance = pandas.DataFrame(final_performance, index=[0])
    perf_path = dir_path + "/" + "fold_" + current_fold + "_model_performance.csv"
    final_performance.to_csv(perf_path)
    
    # Saving predictions and labels
    Y_train = pandas.DataFrame(Y_train)
    Y_train.to_csv(dir_path + "/" + current_fold + "_" + "Y_train.csv")

    Y_test = pandas.DataFrame(Y_test)
    Y_test.to_csv(dir_path + "/" + current_fold + "_" + "Y_test.csv")

    classes_train_rounded = pandas.DataFrame(classes_train_rounded)
    classes_train_rounded.to_csv(dir_path + "/" + current_fold + "_" + "classes_train_rounded.csv")

    classes_test_rounded = pandas.DataFrame(classes_test_rounded)
    classes_test_rounded.to_csv(dir_path + "/" + current_fold + "_" + "classes_test_rounded.csv")

    classes_train = pandas.DataFrame(classes_train)
    classes_train.to_csv(dir_path + "/" + current_fold + "_" + "classes_train.csv")

    classes_test = pandas.DataFrame(classes_test)
    classes_test.to_csv(dir_path + "/" + current_fold + "_" + "classes_test.csv")

#################################################################################
##########################    Section    ########################################
#################################################################################

if FINAL_TEST_ARG == "Optimization":

    print("Script is run in the optimization mode -> analyses are in the CV")

    #### Running folds ####
    kfold = KFold(n_splits=FOLD_ARG, shuffle=True)

    fold_no = 1

    for train, test in kfold.split(meth_tensor, depr_status_tensor):
        # Printing message
        print('--------------------')
        print(f'Training for fold {fold_no} for model {MODEL_ARG}...')

        # Write phenotypes files
        curr_depr_pheno = methylation_data["CV_pheno"]
        curr_depr_pheno_train = curr_depr_pheno.iloc[train]
        curr_depr_pheno_train.to_csv(dir_path + "/" + str(fold_no) + "_" + "pheno_train.csv")

        curr_depr_pheno_test = curr_depr_pheno.iloc[test]
        curr_depr_pheno_test.to_csv(dir_path + "/" + str(fold_no) + "_" + "pheno_test.csv")

        # Run training
        run_ML_model(meth_tensor[train],
                     depr_status_tensor[train],
                     meth_tensor[test],
                     depr_status_tensor[test],
                     fold_no,
                     curr_depr_pheno_train,
                     dir_path,
                     classifier_dict)

        # Increase fold number
        fold_no = fold_no + 1

else:

    print("Script is run in the final testing mode -> CV data is used for trainining")
    print("Testing is performed on last hold-out set")

    # Shuffle samples
    curr_depr_pheno_train = methylation_data["CV_pheno"]
    curr_depr_pheno_test = methylation_data["test_pheno"]

    meth_tensor_train, depr_status_tensor_train, curr_depr_pheno_train = shuffle(meth_tensor_train, 
                                                                                 depr_status_tensor_train,
                                                                                 curr_depr_pheno_train)
    meth_tensor_test, depr_status_tensor_test, curr_depr_pheno_test = shuffle(meth_tensor_test, 
                                                                              depr_status_tensor_test,
                                                                              curr_depr_pheno_test )

    print('--------------------')
    print(f'Training for final test for model {MODEL_ARG}...')

    # Run training
    run_ML_model(meth_tensor_train,
                 depr_status_tensor_train,
                 meth_tensor_test,
                 depr_status_tensor_test,
                 "FINAL_TEST",
                 curr_depr_pheno_train,
                 dir_path,
                 classifier_dict)

    # Write phenotypes files
    curr_depr_pheno_train.to_csv(dir_path + "/" + "FINAL_TEST" + "_" + "pheno_train.csv")
    curr_depr_pheno_test.to_csv(dir_path + "/" + "FINAL_TEST" + "_" + "pheno_test.csv")



#################################################################################
##########################    Section    ########################################
#################################################################################

if FINAL_TEST_ARG == "Optimization":

    #### Averaging stats
    fold_stats = [dir_path + "/" + x for x in os.listdir(dir_path) if "_model_performance" in x]
    data_fr = [pandas.read_csv(x) for x in fold_stats]
    data_fr = pandas.concat(data_fr)

    output_stats_folds =  data_fr[["accuracy_train", "accuracy_test", "auc_train", "auc_test", "roc_auc_plot_test"]].mean()
    output_stats_folds = pandas.DataFrame(output_stats_folds)
    output_stats_folds = output_stats_folds.transpose()

    output_stats_base_info = data_fr.drop(["accuracy_train",
                                        "accuracy_test", 
                                        "auc_train", 
                                        "auc_test",
                                        "roc_auc_plot_test",
                                        "fold"] , axis=1)
    output_stats_base_info = output_stats_base_info.iloc[0]
    output_stats_base_info = pandas.DataFrame(output_stats_base_info)
    output_stats_base_info = output_stats_base_info.transpose()

    output_stats_folds = pandas.concat([output_stats_base_info, output_stats_folds], axis=1)

    average_stat_path = dir_path + "/average_stats.csv"
    output_stats_folds.to_csv(average_stat_path)
else:
    print("Script is run in the final testing mode -> no aggregation required")



#################################################################################
##########################    Section    ########################################
#################################################################################
#### Inspecting distribution of stats

if FINAL_TEST_ARG == "Optimization":

    for x in range(1,FOLD_ARG+1,1):
        helper_functions.explain_predictions_folds(dir_path=dir_path, x=x)

else:

    helper_functions.explain_predictions_folds(dir_path=dir_path, x="FINAL_TEST")
