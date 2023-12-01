#!/bin/bash
#SBATCH -J Bianca_DNN_RAW
#SBATCH -A sens2023584
#SBATCH -t 3-00:00:00
#SBATCH -p node
#SBATCH -N 1
#SBATCH -C gpu

cd /proj/sens2023584/Depression_methylation_multicohort_meta/

module load R_packages/4.2.1
module load python_ML_packages/3.9.5-gpu
export PYTHONPATH=/castor/project/home/aleso209/.local/lib/python3.9/site-packages/:$PYTHONPATH

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
    python3 TrainCV_JointAE.py  "innerDNNClassifier" "10CV_DNN_RAW/10_CV_JointAE/$folder" 2000 3 0.1 \
    "sigmoid" "sigmoid" "linear" 0.0001 "ALL" "Mval" "binary_crossentropy" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data_raw/"
    
    python3 TrainCV_JointAE.py  "innerDNNClassifier" "10CV_DNN_RAW/10_CV_JointAE_consistent/$folder" 2000 3 0.1 \
    "sigmoid" "sigmoid" "linear" 0.0001 "Top_200_all_consistent_CpGs.txt" "Mval" "binary_crossentropy" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data_raw/"

done

for folder in "${subfolders[@]}"
do
    python3  TrainCV_JointVAE.py  "innerDNNClassifier" "10CV_DNN_RAW/10_CV_JointVAE/$folder" 2250 3 1 \
    "tanh" "tanh" "linear" 0.0001 "ALL" "Mval" "squared_hinge" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data_raw/"

    python3 TrainCV_JointVAE.py  "innerDNNClassifier" "10CV_DNN_RAW/10_CV_JointVAE_consistent/$folder" 2250 3 1 \
    "tanh" "tanh" "linear" 0.0001 "Top_200_all_consistent_CpGs.txt" "Mval" "squared_hinge" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data_raw/"

done

for folder in "${subfolders[@]}"
do
    python3   TrainCV_SimpleDNN.py "DNNclassfierSimpleBatch" "10CV_DNN_RAW/10_CV_SmallDNN/$folder" 1000 3 "relu" \
    "sigmoid" 0.0001 "ALL" "Mval" "binary_crossentropy" "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data_raw/"

    python3   TrainCV_SimpleDNN.py "DNNclassfierSimpleBatch" "10CV_DNN_RAW/10_CV_SmallDNN_consistent/$folder" 1000 3 "relu" \
    "sigmoid" 0.0001 "Top_200_all_consistent_CpGs.txt" "Mval" "binary_crossentropy" "Not_save" 128 \
    "Optimization" "all_cohorts" 200 "Main_analysis/ML_data_raw/"

done

echo "Final test for DNN models"

python3 TrainCV_JointAE.py  "innerDNNClassifier" "FINAL_TEST_DNN_RAW/JointAE" 2000 3 0.1 \
    "sigmoid" "sigmoid" "linear" 0.0001 "ALL" "Mval" "binary_crossentropy" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data_raw/"

python3 TrainCV_JointAE.py  "innerDNNClassifier" "FINAL_TEST_DNN_RAW/JointAE" 2000 3 0.1 \
    "sigmoid" "sigmoid" "linear" 0.0001 "Top_200_all_consistent_CpGs.txt" "Mval" "binary_crossentropy" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data_raw/"


python3 TrainCV_JointVAE.py  "innerDNNClassifier" "FINAL_TEST_DNN_RAW/JointVAE" 2250 3 1 \
    "tanh" "tanh" "linear" 0.0001 "ALL" "Mval" "squared_hinge" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data_raw/"

python3 TrainCV_JointVAE.py  "innerDNNClassifier" "FINAL_TEST_DNN_RAW/JointVAE" 2250 3 1 \
    "tanh" "tanh" "linear" 0.0001 "Top_200_all_consistent_CpGs.txt" "Mval" "squared_hinge" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data_raw/"


python3 TrainCV_SimpleDNN.py "DNNclassfierSimpleBatch" "FINAL_TEST_DNN_RAW/SmallDNN" 1000 3 "relu" \
    "sigmoid" 0.0001 "ALL" "Mval" "binary_crossentropy" "Not_save" 128 "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data_raw/"

python3 TrainCV_SimpleDNN.py "DNNclassfierSimpleBatch" "FINAL_TEST_DNN_RAW/SmallDNN" 1000 3 "relu" \
    "sigmoid" 0.0001 "Top_200_all_consistent_CpGs.txt" "Mval" "binary_crossentropy" "Not_save" 128 \
    "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data_raw/"



cp -r /proj/sens2023584/Depression_methylation_multicohort_meta/10CV_DNN_RAW /proj/sens2023584/nobackup/wharf/aleso209/aleso209-sens2023584/10CV_DNN_RAW
cp -r /proj/sens2023584/Depression_methylation_multicohort_meta/FINAL_TEST_DNN_RAW /proj/sens2023584/nobackup/wharf/aleso209/aleso209-sens2023584/FINAL_TEST_DNN_RAW