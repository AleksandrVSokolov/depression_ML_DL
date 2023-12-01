#!/bin/bash
#SBATCH -J Bianca_JointVAE
#SBATCH -A sens2023584
#SBATCH -t 3-00:00:00
#SBATCH -p node
#SBATCH -N 1
#SBATCH -C gpu

cd /proj/sens2023584/Depression_methylation_multicohort_meta/

module load R_packages/4.2.1
module load python_ML_packages/3.9.5-gpu
export PYTHONPATH=/castor/project/home/aleso209/.local/lib/python3.9/site-packages/:$PYTHONPATH


declare -a alphas=(0.5
1
3
5
10)

declare -a CpGtypes=("ALL"
"Top_200_all_consistent_CpGs.txt")

declare -a drops=(0.0
0.1)


for alpha in "${alphas[@]}"
do
   for value in "${CpGtypes[@]}"
    do
        for drop in "${drops[@]}"
        do
            python3 TrainCV_JointVAE.py  "innerDNNClassifier" "test_grid_jointVAE_BIANCA" 2000 3 "$alpha" \
            "sigmoid" "sigmoid" "sigmoid" 0.0001 "$value" "Betaval" "binary_crossentropy" "binary_crossentropy" "$drop" "BN" "l1_l2" \
                "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"
            
            python3 TrainCV_JointVAE.py  "innerDNNClassifier" "test_grid_jointVAE_BIANCA" 3000 3 "$alpha" \
            "sigmoid" "sigmoid" "sigmoid" 0.0001 "$value" "Betaval" "binary_crossentropy" "binary_crossentropy" "$drop" "BN" "l1_l2" \
                "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"


            python3 TrainCV_JointVAE.py  "innerDNNClassifier" "test_grid_jointVAE_BIANCA" 2000 3 "$alpha" \
            "relu" "sigmoid" "sigmoid" 0.0001 "$value" "Betaval" "binary_crossentropy" "binary_crossentropy" "$drop" "BN" "l1_l2" \
                "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

            python3 TrainCV_JointVAE.py  "innerDNNClassifier" "test_grid_jointVAE_BIANCA" 2000 3 "$alpha" \
            "sigmoid" "sigmoid" "linear" 0.0001 "$value" "Mval" "binary_crossentropy" "mean_squared_error" "$drop" "BN" "l1_l2" \
                "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

            python3 TrainCV_JointVAE.py  "innerDNNClassifier" "test_grid_jointVAE_BIANCA" 2000 3 "$alpha" \
            "tanh" "sigmoid" "linear" 0.0001 "$value" "Mval" "binary_crossentropy" "mean_squared_error" "$drop" "no_BN" "l1_l2" \
                "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

            python3 TrainCV_JointVAE.py  "innerDNNClassifier" "test_grid_jointVAE_BIANCA" 2000 3 "$alpha" \
            "relu" "sigmoid" "linear" 0.0001 "$value" "Mval" "binary_crossentropy" "mean_squared_error" "$drop" "BN" "l1_l2" \
                "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

            python3 TrainCV_JointVAE.py  "innerDNNClassifier" "test_grid_jointVAE_BIANCA" 2000 3 "$alpha" \
            "tanh" "tanh" "linear" 0.0001 "$value" "Mval" "squared_hinge" "mean_squared_error" "$drop" "BN" "l1_l2" \
                "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

            python3 TrainCV_JointVAE.py  "innerDNNClassifier" "test_grid_jointVAE_BIANCA" 2000 3 "$alpha" \
            "tanh" "tanh" "linear" 0.0001 "$value" "Mval" "squared_hinge" "mean_squared_error" "$drop" "no_BN" "l1_l2" \
                "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"
        
        done
    done
done

cp -r /proj/sens2023584/Depression_methylation_multicohort_meta/test_grid_jointVAE_BIANCA /proj/sens2023584/nobackup/wharf/aleso209/aleso209-sens2023584/test_grid_jointVAE_BIANCA