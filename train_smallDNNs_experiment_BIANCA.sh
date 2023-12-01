#!/bin/bash
#SBATCH -J Bianca_smallDNN
#SBATCH -A sens2023584
#SBATCH -t 3-00:00:00
#SBATCH -p node
#SBATCH -N 1
#SBATCH -C gpu

cd /proj/sens2023584/Depression_methylation_multicohort_meta/

module load R_packages/4.2.1
module load python_ML_packages/3.9.5-gpu
export PYTHONPATH=/castor/project/home/aleso209/.local/lib/python3.9/site-packages/:$PYTHONPATH

declare -a models=("DNNclassfier" 
    "DNNclassfiernoreg" 
    "DNNclassfierL1L2" 
    "DNNclassfierL1L2drop" 
    "DNNclassfierL1L2manydrops"
    "DNNclassfierL1L2manyBigdrops"
    "DNNclassfierOnlyDrops" 
    "DNNclassfierSimple" 
    "DNNclassfierSimpleBatch")

declare -a inner_act=("relu" 
    "sigmoid" 
    "tanh")


declare -a CpGtypes=("ALL"
"Top_200_all_consistent_CpGs.txt")


for model in "${models[@]}"
do
   for activation in "${inner_act[@]}"
    do

        for CpG in "${CpGtypes[@]}"
        do
            python3 TrainCV_SimpleDNN.py "$model" "test_grid_simpleDNN_BIANCA" 2000 3 "$activation" \
            "sigmoid" 0.0001 "$CpG" "Betaval" "binary_crossentropy" "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

            python3 TrainCV_SimpleDNN.py "$model" "test_grid_simpleDNN_BIANCA" 2000 3 "$activation" \
            "sigmoid" 0.0001 "$CpG" "Mval" "binary_crossentropy" "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

            python3 TrainCV_SimpleDNN.py "$model" "test_grid_simpleDNN_BIANCA" 2000 3 "$activation" \
            "sigmoid" 0.00001 "$CpG" "Betaval" "binary_crossentropy" "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

            python3 TrainCV_SimpleDNN.py "$model" "test_grid_simpleDNN_BIANCA" 2000 3 "$activation" \
            "sigmoid" 0.00001 "$CpG" "Mval" "binary_crossentropy" "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"
        done
    done
done

cp -r /proj/sens2023584/Depression_methylation_multicohort_meta/test_grid_simpleDNN_BIANCA /proj/sens2023584/nobackup/wharf/aleso209/aleso209-sens2023584/test_grid_simpleDNN_BIANCA