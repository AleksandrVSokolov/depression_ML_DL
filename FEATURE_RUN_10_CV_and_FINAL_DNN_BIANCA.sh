#!/bin/bash
#SBATCH -J FEATURE_Bianca_DNN
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

declare -a selectors=("selector_variance_5"
"selector_variance_1"
"selector_variance_01"
"selector_k_best_f_200"
"selector_lin_svc_200"
"selector_log_reg_200"
"selector_trees_200")

for selector in "${selectors[@]}"
do
    for folder in "${subfolders[@]}"
    do
        python3 TrainCV_JointAE.py  "innerDNNClassifier" "FEATURE_10CV_DNN/10_CV_JointAE/$selector/$folder" 2000 3 0.1 \
        "sigmoid" "sigmoid" "linear" 0.0001 "$selector" "Mval" "binary_crossentropy" "mean_squared_error" 0.1 "BN" "l1_l2" \
        "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

    done
done


for selector in "${selectors[@]}"
do
    for folder in "${subfolders[@]}"
    do
        python3  TrainCV_JointVAE.py  "innerDNNClassifier" "FEATURE_10CV_DNN/10_CV_JointVAE/$selector/$folder" 2250 3 1 \
        "tanh" "tanh" "linear" 0.0001 "$selector" "Mval" "squared_hinge" "mean_squared_error" 0.1 "BN" "l1_l2" \
        "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

    done
done


for selector in "${selectors[@]}"
do
    for folder in "${subfolders[@]}"
    do
        python3   TrainCV_SimpleDNN.py "DNNclassfierSimpleBatch" "FEATURE_10CV_DNN/10_CV_SmallDNN/$selector/$folder" 1000 3 "relu" \
        "sigmoid" 0.0001 "$selector" "Mval" "binary_crossentropy" "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

    done
done

echo "Final test for DNN models"

for selector in "${selectors[@]}"
do
python3 TrainCV_JointAE.py  "innerDNNClassifier" "FEATURE_FINAL_TEST_DNN/JointAE" 2000 3 0.1 \
    "sigmoid" "sigmoid" "linear" 0.0001 "$selector" "Mval" "binary_crossentropy" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data/"

python3 TrainCV_JointVAE.py  "innerDNNClassifier" "FEATURE_FINAL_TEST_DNN/JointVAE" 2250 3 1 \
    "tanh" "tanh" "linear" 0.0001 "$selector" "Mval" "squared_hinge" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data/"

python3 TrainCV_SimpleDNN.py "DNNclassfierSimpleBatch" "FEATURE_FINAL_TEST_DNN/SmallDNN" 1000 3 "relu" \
    "sigmoid" 0.0001 "$selector" "Mval" "binary_crossentropy" "Not_save" 128 "FINAL_TEST" "all_cohorts" 200 "Main_analysis/ML_data/"
done


cp -r /proj/sens2023584/Depression_methylation_multicohort_meta/FEATURE_10CV_DNN /proj/sens2023584/nobackup/wharf/aleso209/aleso209-sens2023584/FEATURE_10CV_DNN
cp -r /proj/sens2023584/Depression_methylation_multicohort_meta/FEATURE_FINAL_TEST_DNN /proj/sens2023584/nobackup/wharf/aleso209/aleso209-sens2023584/FEATURE_FINAL_TEST_DNN