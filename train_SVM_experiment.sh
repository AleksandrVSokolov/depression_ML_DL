#!/bin/bash

<<comment

conda activate tf-py38

cd /home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/

chmod +x train_SVM_experiment.sh

./train_SVM_experiment.sh

comment


declare -a models=("SVM_c"
"SVM_c_1"
"SVM_c_2"
"SVM_c_3"
"SVM_c_4"
"SVM_c_5"
"SVM_c_6"
"SVM_c_RBF"
"SVM_c_RBF_1"
"SVM_c_RBF_2"
"SVM_c_RBF_3"
"SVM_c_RBF_4"
"SVM_c_RBF_5"
"SVM_c_RBF_6"
"SVM_c_RBF_7"
"SVM_c_RBF_8"
"SVM_c_RBF_9"
"SVM_c_RBF_10"
"SVM_c_RBF_11"
"SVM_c_RBF_12"
"SVM_c_RBF_13"
"SVM_c_RBF_14"
"SVM_c_RBF_15"
"SVM_c_RBF_16")

for model in "${models[@]}"
do
    /home/aleksandr/miniconda3/envs/tf-py38/bin/python  TrainCV_ML.py  "$model" "test_grid_SVM" 3 \
    "ALL" "Mval" "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

done
