#!/bin/bash
<<comment

conda activate tf-py38

cd /home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta/

chmod +x epoch_num_experiment.sh

./epoch_num_experiment.sh

comment



declare -a epochs=(250
500
750
1000
1250
1500
1750
2000
2250
2500)


for ep in "${epochs[@]}"
do
    /home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_SimpleDNN.py "DNNclassfierSimpleBatch" "test_grid_epochs/smallDNNs" "$ep" 3 "relu" \
    "sigmoid" 0.0001 "ALL" "Mval" "binary_crossentropy" "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

    /home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_JointAE.py  "innerDNNClassifier" "test_grid_epochs/JointAE" "$ep" 3 0.1 \
    "sigmoid" "sigmoid" "linear" 0.0001 "ALL" "Mval" "binary_crossentropy" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

    /home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_JointVAE.py  "innerDNNClassifier" "test_grid_epochs/JointVAE" "$ep" 3 1 \
    "tanh" "tanh" "linear" 0.0001 "ALL" "Mval" "squared_hinge" "mean_squared_error" 0.1 "BN" "l1_l2" \
    "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

done