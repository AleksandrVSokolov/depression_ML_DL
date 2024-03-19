#!/bin/bash

<<comment

conda activate tf-py38

cd /home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/

chmod +x FEATURE_RUN_10_CV_ML.sh

./FEATURE_RUN_10_CV_ML.sh

comment

declare -a subfolders=("1"
"2"
"3"
"4"
"5"
"6"
"7"
"8"
"9"
"10")

declare -a models=("LogisticRegr"
"SVM_c"
"SVM_c_RBF_15"
"DecisionTree"
"RandomForest"
"AdaBoost"
"Elastic_net_classifier"
"RidgeLogisticRegr"
"Lasso_classifier")

declare -a selectors=("selector_variance_5"
"selector_variance_1"
"selector_variance_01"
"selector_k_best_f_200"
"selector_lin_svc_200"
"selector_log_reg_200"
"selector_trees_200")

echo "Running ML models with selectors in CV"

for model in "${models[@]}"
do
    for selector in "${selectors[@]}"
    do
        for folder in "${subfolders[@]}"
        do
            /home/aleksandr/miniconda3/envs/tf-py38/bin/python  TrainCV_ML.py  "$model" "10_CV_ML_Selectors/$model/$selector/$folder" 3 \
            "$selector" "Mval" "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

        done
    done
done

for model in "${models[@]}"
do
    for selector in "${selectors[@]}"
    do
        for folder in "${subfolders[@]}"
        do
            /home/aleksandr/miniconda3/envs/tf-py38/bin/python  TrainCV_ML.py  "$model" "10_CV_ML_Selectors_RAW/$model/$selector/$folder" 3 \
            "$selector" "Mval" "Optimization" "all_cohorts" 200 "Main_analysis/ML_data_raw/"

        done
    done
done


echo "Running FINAL test for ML models with selectors"

for model in "${models[@]}"
do
    for selector in "${selectors[@]}"
    do
        /home/aleksandr/miniconda3/envs/tf-py38/bin/python  TrainCV_ML.py  "$model" "FINAL_TEST_ML_Selectors" 3 \
        "$selector" "Mval" "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data/"

    done
done

for model in "${models[@]}"
do
    for selector in "${selectors[@]}"
    do
        /home/aleksandr/miniconda3/envs/tf-py38/bin/python  TrainCV_ML.py  "$model" "FINAL_TEST_ML_Selectors_RAW" 3 \
        "$selector" "Mval" "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data_raw/"

    done
done

