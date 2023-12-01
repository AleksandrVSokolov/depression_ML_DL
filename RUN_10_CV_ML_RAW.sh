#!/bin/bash

<<comment

conda activate tf-py38

cd /home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/

chmod +x RUN_10_CV_ML_RAW.sh

./RUN_10_CV_ML_RAW.sh

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
"SVM_c_RBF"
"DecisionTree"
"RandomForest"
"AdaBoost"
"Elastic_net_classifier"
"RidgeLogisticRegr"
"Lasso_classifier")

echo "Running ML models with top 200 CpGs from Limma"

for model in "${models[@]}"
do
    for folder in "${subfolders[@]}"
    do
        /home/aleksandr/miniconda3/envs/tf-py38/bin/python  TrainCV_ML.py  "$model" "ML_CV_and_FINAL_RAW/10_CV_ML_Limma_CpGs/$model/$folder" 3 \
        "ALL" "Mval" "Optimization" "all_cohorts" 200 "Main_analysis/ML_data_raw/"

        /home/aleksandr/miniconda3/envs/tf-py38/bin/python  TrainCV_ML.py  "$model" "ML_CV_and_FINAL_RAW/10_CV_ML_Limma_CpGs/$model/$folder" 3 \
        "ALL" "Betaval" "Optimization" "all_cohorts" 200 "Main_analysis/ML_data_raw/"

    done
done

echo "Running ML models with Top 200_ALL_consistentonsistent CpGs"

for model in "${models[@]}"
do
    for folder in "${subfolders[@]}"
    do
        /home/aleksandr/miniconda3/envs/tf-py38/bin/python  TrainCV_ML.py  "$model" "ML_CV_and_FINAL_RAW/10_CV_ML_200_ALL_consistent_CpGs/$model/$folder" 3 \
        "Top_200_all_consistent_CpGs.txt" "Mval" "Optimization" "all_cohorts" 200 "Main_analysis/ML_data_raw/"

        /home/aleksandr/miniconda3/envs/tf-py38/bin/python  TrainCV_ML.py  "$model" "ML_CV_and_FINAL_RAW/10_CV_ML_200_ALL_consistent_CpGs/$model/$folder" 3 \
        "Top_200_all_consistent_CpGs.txt" "Betaval" "Optimization" "all_cohorts" 200 "Main_analysis/ML_data_raw/"

    done
done

echo "Running ML models with Top 10000 Limma CpGs and penalties (as feature selection)"

declare -a models=("Elastic_net_classifier"
"RidgeLogisticRegr"
"Lasso_classifier")

for model in "${models[@]}"
do
    for folder in "${subfolders[@]}"
    do
        /home/aleksandr/miniconda3/envs/tf-py38/bin/python  TrainCV_ML.py  "$model" "ML_CV_and_FINAL_RAW/10_CV_ML_10K_Limma/$model/$folder" 3 \
        "ALL" "Mval" "Optimization" "all_cohorts" 10000 "Main_analysis/ML_data_raw/"

        /home/aleksandr/miniconda3/envs/tf-py38/bin/python  TrainCV_ML.py  "$model" "ML_CV_and_FINAL_RAW/10_CV_ML_10K_Limma/$model/$folder" 3 \
        "ALL" "Betaval" "Optimization" "all_cohorts" 10000 "Main_analysis/ML_data_raw/"

    done
done


echo "Running FINAL test for ML models with significant and consistent CpGs"

declare -a models=("LogisticRegr"
"SVM_c"
"SVM_c_RBF"
"DecisionTree"
"RandomForest"
"AdaBoost"
"Elastic_net_classifier"
"RidgeLogisticRegr"
"Lasso_classifier")

declare -a CpGtypes=("ALL"
"Top_200_all_consistent_CpGs.txt")


for model in "${models[@]}"
do
    for CpG in "${CpGtypes[@]}"
    do
        /home/aleksandr/miniconda3/envs/tf-py38/bin/python  TrainCV_ML.py  "$model" "ML_CV_and_FINAL_RAW/FINAL_TEST_ML" 3 \
        "$CpG" "Mval" "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data_raw/"

        /home/aleksandr/miniconda3/envs/tf-py38/bin/python  TrainCV_ML.py  "$model" "ML_CV_and_FINAL_RAW/FINAL_TEST_ML" 3 \
        "$CpG" "Betaval" "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data_raw/"

    done
done

echo "Running ML models with Top 10000 Limma CpGs and penalties (as feature selection) FINAL TEST"

declare -a models=("Elastic_net_classifier"
"RidgeLogisticRegr"
"Lasso_classifier")

for model in "${models[@]}"
do
    /home/aleksandr/miniconda3/envs/tf-py38/bin/python  TrainCV_ML.py  "$model" "ML_CV_and_FINAL_RAW/FINAL_TEST_ML_10K/$model" 3 \
    "ALL" "Mval" "FINAL_TEST" "all_cohorts" 10000 "Main_analysis/ML_data_raw/"

    /home/aleksandr/miniconda3/envs/tf-py38/bin/python  TrainCV_ML.py  "$model" "ML_CV_and_FINAL_RAW/FINAL_TEST_ML_10K/$model" 3 \
    "ALL" "Betaval" "FINAL_TEST" "all_cohorts" 10000 "Main_analysis/ML_data_raw/"
done
