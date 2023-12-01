#!/bin/bash
<<comment

conda activate tf-py38

cd /home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/

chmod +x RUN_10_CV_and_FINAL_DNN.sh

./RUN_10_CV_and_FINAL_DNN.sh

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

for folder in "${subfolders[@]}"
do
    /home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_JointAE.py  "innerDNNClassifier" "10CV_DNN/10_CV_JointAE/$folder" 2000 3 0.1 \
    "sigmoid" "sigmoid" "linear" 0.0001 "ALL" "Mval" "binary_crossentropy" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"
    
    /home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_JointAE.py  "innerDNNClassifier" "10CV_DNN/10_CV_JointAE_consistent/$folder" 2000 3 0.1 \
    "sigmoid" "sigmoid" "linear" 0.0001 "Top_200_all_consistent_CpGs.txt" "Mval" "binary_crossentropy" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

done

for folder in "${subfolders[@]}"
do
    /home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_JointVAE.py  "innerDNNClassifier" "10CV_DNN/10_CV_JointVAE/$folder" 2250 3 1 \
    "tanh" "tanh" "linear" 0.0001 "ALL" "Mval" "squared_hinge" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

    /home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_JointVAE.py  "innerDNNClassifier" "10CV_DNN/10_CV_JointVAE_consistent/$folder" 2250 3 1 \
    "tanh" "tanh" "linear" 0.0001 "Top_200_all_consistent_CpGs.txt" "Mval" "squared_hinge" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

done

for folder in "${subfolders[@]}"
do
    /home/aleksandr/miniconda3/envs/tf-py38/bin/python   TrainCV_SimpleDNN.py "DNNclassfierSimpleBatch" "10CV_DNN/10_CV_SmallDNN/$folder" 1000 3 "relu" \
    "sigmoid" 0.0001 "ALL" "Mval" "binary_crossentropy" "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

    /home/aleksandr/miniconda3/envs/tf-py38/bin/python   TrainCV_SimpleDNN.py "DNNclassfierSimpleBatch" "10CV_DNN/10_CV_SmallDNN_consistent/$folder" 1000 3 "relu" \
    "sigmoid" 0.0001 "Top_200_all_consistent_CpGs.txt" "Mval" "binary_crossentropy" "Not_save" 128 \
    "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

done

echo "Final test for DNN models"

/home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_JointAE.py  "innerDNNClassifier" "FINAL_TEST_DNN/JointAE" 2000 3 0.1 \
    "sigmoid" "sigmoid" "linear" 0.0001 "ALL" "Mval" "binary_crossentropy" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data/"

/home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_JointAE.py  "innerDNNClassifier" "FINAL_TEST_DNN/JointAE" 2000 3 0.1 \
    "sigmoid" "sigmoid" "linear" 0.0001 "Top_200_all_consistent_CpGs.txt" "Mval" "binary_crossentropy" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data/"


/home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_JointVAE.py  "innerDNNClassifier" "FINAL_TEST_DNN/JointVAE" 2250 3 1 \
    "tanh" "tanh" "linear" 0.0001 "ALL" "Mval" "squared_hinge" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data/"

/home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_JointVAE.py  "innerDNNClassifier" "FINAL_TEST_DNN/JointVAE" 2250 3 1 \
    "tanh" "tanh" "linear" 0.0001 "Top_200_all_consistent_CpGs.txt" "Mval" "squared_hinge" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data/"


/home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_SimpleDNN.py "DNNclassfierSimpleBatch" "FINAL_TEST_DNN/SmallDNN" 1000 3 "relu" \
    "sigmoid" 0.0001 "ALL" "Mval" "binary_crossentropy" "Not_save" 128 "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data/"

/home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_SimpleDNN.py "DNNclassfierSimpleBatch" "FINAL_TEST_DNN/SmallDNN" 1000 3 "relu" \
    "sigmoid" 0.0001 "Top_200_all_consistent_CpGs.txt" "Mval" "binary_crossentropy" "Not_save" 128 \
    "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data/"