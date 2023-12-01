# Importing modules
import sys
import os
sys.path.append(os.getcwd())

import numpy
import pandas
import os
import gc
from HelperScriptsPY import helper_functions
import time
import tensorflow as tf
import matplotlib.pyplot as plt
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score, roc_auc_score, auc, roc_curve
from sklearn.utils import shuffle
from Models.ClassifierModels import innerDNNClassifier

"""
Example code:

conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta
/home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_JointAE.py  "innerDNNClassifier" "test_folder_jointAE" 100 3 0.3 \
    "sigmoid" "sigmoid" "sigmoid" 0.0001 "Top_200_all_consistent_CpGs.txt" "Betaval" "binary_crossentropy" "binary_crossentropy" 0.1 "BN" "l1_l2" \
        "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"


conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta
/home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_JointAE.py  "innerDNNClassifier" "test_folder_jointAE" 2000 3 0.3 \
    "sigmoid" "sigmoid" "sigmoid" 0.0001 "ALL" "Betaval" "binary_crossentropy" "binary_crossentropy" 0.1 "BN" "l1_l2" \
        "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/Depression_methylation_multicohort_meta
/home/aleksandr/miniconda3/envs/tf-py38/bin/python TrainCV_JointAE.py  "innerDNNClassifier" "test_folder_jointAE" 550 3 0.3 \
    "tanh" "tanh" "linear" 0.0001 "ALL" "Mval" "squared_hinge" "mean_squared_error" 0.1 "no_BN" "l1_l2" \
        "Not_save" 128 "Optimization" "all_cohorts" 200 "Main_analysis/ML_data/"

"""

#### Parsing arguments ####
ARG_LIST = sys.argv

# Model parameters
MODEL_ARG = str(ARG_LIST[1])
FOLDER_ARG = str(ARG_LIST[2])
EPOCH_ARG = int(ARG_LIST[3])
FOLD_ARG = int(ARG_LIST[4])
ALPHA_ARG = float(ARG_LIST[5])
INNER_ARG = str(ARG_LIST[6])
OUTER_CLF_ARG = str(ARG_LIST[7])
OUTER_AE_ARG = str(ARG_LIST[8])
LR_ARG = float(ARG_LIST[9])
CONSIST_ARG = str(ARG_LIST[10])
VALUES_ARG = str(ARG_LIST[11])
CLASS_LOSS_ARG = str(ARG_LIST[12])
RECON_LOSS_ARG = str(ARG_LIST[13])
DROP_RATE_ARG = float(ARG_LIST[14])
BN_ARG = str(ARG_LIST[15])
KR_ARG = str(ARG_LIST[16])

# General arguments
WEIGHT_ARG = str(ARG_LIST[17])
BATCH_ARG = int(ARG_LIST[18])
FINAL_TEST_ARG = str(ARG_LIST[19])
COHORT_POOL_ARG = str(ARG_LIST[20])
LIMMA_FEATURE_ARG = int(ARG_LIST[21])
DATA_ARG = str(ARG_LIST[22])

print('------------------------------------------------------------------------')
print(f'Training JointAE model with inner classifier {MODEL_ARG}')
print(f'Training for epochs {EPOCH_ARG} for n-folds {FOLD_ARG}')
print(f'Balancing alpha: {ALPHA_ARG}')
print(f'Inner activation: {INNER_ARG}')
print(f'Outer activation (classifier): {OUTER_CLF_ARG}')
print(f'Outer activation (autoencoder): {OUTER_AE_ARG}')
print(f'Learning rate: {LR_ARG}')
print(f'Selected CpGs: {CONSIST_ARG}')
print(f'Selected values: {VALUES_ARG}')
print(f'Classification loss: {CLASS_LOSS_ARG}')
print(f'Reconstruction loss: {RECON_LOSS_ARG}')
print(f'Drop rate: {DROP_RATE_ARG}')
print(f'Batch normalization: {BN_ARG}')
print(f'Kernel regularization: {KR_ARG}')
print(f'Saving model: {WEIGHT_ARG}')
print(f'Using batch: {BATCH_ARG}')
print(f'Limma selected features: {LIMMA_FEATURE_ARG}')
print(f'Perform final test: {FINAL_TEST_ARG}')

# Preparing target directory
al_str = str(ALPHA_ARG)
al_str = al_str.replace(".","_")
al_str = "_alpha_" + al_str  +  "_"

lr_str = str(LR_ARG)
lr_str = lr_str.replace(".","_")
lr_str = "_lr_" + lr_str  +  "_"

dr_str = str(DROP_RATE_ARG)
dr_str = dr_str.replace(".","_")
dr_str = "_drop_" + dr_str  +  "_"

cons_str = "_cpg_" + CONSIST_ARG
cons_str = cons_str.replace(".txt","")

dir_path = (FOLDER_ARG + "/" + MODEL_ARG + "_ep_" + str(EPOCH_ARG) + "_fo_" + str(FOLD_ARG) + al_str +
            "_in_" + INNER_ARG + "_ouc_" + OUTER_CLF_ARG + "_oue_" + OUTER_AE_ARG + lr_str + 
            cons_str + "_val_" + VALUES_ARG + "_closs_" + CLASS_LOSS_ARG + 
            "_rloss_" + RECON_LOSS_ARG + dr_str + "_bn_" + BN_ARG + "_kr_" + KR_ARG + "_ba_" + str(BATCH_ARG) +
            "_limma_" + str(LIMMA_FEATURE_ARG) + 
            "___" + FINAL_TEST_ARG)

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

# Classifier dict to send to a function
classifier_dict = {
    "innerDNNClassifier": innerDNNClassifier
}

# Preparing tensorflow and GPU
os.environ["TF_GPU_ALLOCATOR"] = "cuda_malloc_async"

# Setting memory for tensorflow
gpus = tf.config.list_physical_devices('GPU')
for gpu in gpus:
    tf.config.experimental.set_memory_growth(gpu, True)


#################################################################################
##########################    Section    ########################################
#################################################################################
# Defining classification loss
if CLASS_LOSS_ARG == "binary_crossentropy":
    print("Using BCE loss for classification")
    @tf.function
    def classif_loss_fn(y_true, y_pred):
        bce = tf.keras.losses.BinaryCrossentropy()
        loss = bce(y_true, y_pred)
        return loss
elif CLASS_LOSS_ARG == "binary_focal_crossentropy":
    print("Using BFCE loss for classification")
    @tf.function
    def classif_loss_fn(y_true, y_pred):
        bce = tf.keras.losses.BinaryFocalCrossentropy()
        loss = bce(y_true, y_pred)
        return loss
elif CLASS_LOSS_ARG == "hinge":
    print("Using Hinge loss for classification")
    @tf.function
    def classif_loss_fn(y_true, y_pred):
        h = tf.keras.losses.Hinge()
        loss = h(y_true, y_pred)
        return loss
elif CLASS_LOSS_ARG == "squared_hinge":
    print("Using SquaredHinge loss for classification")
    @tf.function
    def classif_loss_fn(y_true, y_pred):
        h2 = tf.keras.losses.SquaredHinge()
        loss = h2(y_true, y_pred)
        return loss

# Defining reconstruction loss
if RECON_LOSS_ARG == "mean_squared_error":
    print("Using MSE loss for reconstruction")
    recon_loss_fn = tf.keras.losses.MeanSquaredError()

elif RECON_LOSS_ARG == "binary_crossentropy":
    print("Using BCE loss for reconstruction")
    @tf.function
    def recon_loss_fn(y_true, y_pred):
        bce = tf.keras.losses.BinaryFocalCrossentropy()
        loss = bce(y_true, y_pred)
        return loss
    


#################################################################################
##########################    Section    ########################################
#################################################################################

def run_TF_model(X_train,
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

    # Encoder part
    latent_space = 32

    encoder_input = tf.keras.layers.Input(shape=X_train.shape[1])
    x = tf.keras.layers.Dense(128, activation=INNER_ARG, name="dense_1")(encoder_input)
    x = tf.keras.layers.Dense(64, activation=INNER_ARG, name="dense_2")(x)
    bottleneck = tf.keras.layers.Dense(latent_space, activation=INNER_ARG, name="bottleneck")(x)
    encoder = tf.keras.Model(encoder_input, bottleneck, name="encoder")

    # Decoder part
    output_shape = X_train.shape[1]

    latent_inputs = tf.keras.layers.Input(shape=latent_space)
    y = tf.keras.layers.Dense(64, activation=INNER_ARG, name="dense_3")(latent_inputs)
    y = tf.keras.layers.Dense(128, activation=INNER_ARG, name="dense_4")(y)
    output = tf.keras.layers.Dense(output_shape, activation=OUTER_AE_ARG, name="output")(y)
    decoder = tf.keras.Model(latent_inputs, output, name="decoder")

    # Classifier part
    if BN_ARG == "BN":
        bn_param = True
    elif BN_ARG == "no_BN":
        bn_param = False
    else:
        ValueError("Incorrect batch normalization parameter")

    if KR_ARG not in ["none", 'l1', 'l2', "l1_l2"]:
        ValueError("Incorrect kernel regularization parameter")
    elif KR_ARG == "none":
        kr_param = None
    else:
        kr_param = KR_ARG

    classifier_input = tf.keras.layers.Input(shape = latent_space)
    dnn_model_classifier = tf.keras.Model(inputs = classifier_input, 
                                          outputs = classifier_dict[MODEL_ARG](input_dim1D=classifier_input,
                                                                               BatchNorm=bn_param,
                                                                               inner_activations=INNER_ARG,
                                                                               outer_activation=OUTER_CLF_ARG,
                                                                               DropRate=DROP_RATE_ARG,
                                                                               kernel_regularizer=kr_param))
    
    class AEclassifier(tf.keras.Model):
        """
        Dense auto-encoder and classifier
        """
        def __init__(self, 
                     encoder, 
                     decoder, 
                     dnn_model_classifier, 
                     recon_loss,
                     classif_loss,
                     balance_alpha, **kwargs):
            super().__init__(**kwargs)
            self.encoder = encoder
            self.decoder = decoder
            self.classifier = dnn_model_classifier
            self.recon_loss = recon_loss
            self.classif_loss = classif_loss
            self.alpha = balance_alpha
        
        def train_step(self, data):

            # Data unpacking
            x, y = data

            with tf.GradientTape() as tape:

                encoded_space = self.encoder(x)

                reconstruction = self.decoder(encoded_space)

                y_pred = self.classifier(encoded_space)
                
                reconstruction_loss = self.recon_loss(x, reconstruction)

                classification_loss = self.classif_loss(y, y_pred)

                total_loss = self.alpha*reconstruction_loss + (1 - self.alpha)*classification_loss
                
            grads = tape.gradient(total_loss, self.trainable_weights)
            self.optimizer.apply_gradients(zip(grads, self.trainable_weights))

            return {
                "total_loss": total_loss,
                "reconstruction_loss": reconstruction_loss,
                "class_loss": classification_loss
            }
    
    # Model compilation
    joint_classif = AEclassifier(encoder,
                                 decoder, 
                                 dnn_model_classifier, 
                                 recon_loss=recon_loss_fn,
                                 classif_loss=classif_loss_fn,
                                 balance_alpha=ALPHA_ARG)
    

        
    optimizer_custom = tf.keras.optimizers.Adam(learning_rate=LR_ARG)
    joint_classif.compile(optimizer = optimizer_custom)

    # Model training
    ae_classif_history = joint_classif.fit(X_train,
                                    Y_train,
                                    epochs=EPOCH_ARG,
                                    batch_size=BATCH_ARG,
                                    shuffle=True)
    
    train_figure_path = dir_path + "/" + "fold_" + current_fold + "_training.pdf"
    train_history_path = dir_path+ "/" + "fold_" + current_fold + "_training_hist.csv"
    
    train_figure_path = dir_path + "/" + "fold_" + current_fold + "_training.pdf"
    train_history_path = dir_path+ "/" + "fold_" + current_fold + "_training_hist.csv"

    pandas.DataFrame(ae_classif_history.history).plot()
    plt.savefig(train_figure_path) 
    pandas.DataFrame(ae_classif_history.history).to_csv(train_history_path)
    
    # Saving weights
    if WEIGHT_ARG == "Save":
        print("Saving weights")
        t = time.localtime()
        current_time = time.strftime("%y-%m-%d-%H-%M", t)
        path_AE = dir_path + "/" + "Joined_AE_weights_" + current_time + "_fold_" + current_fold  + "_.h5"
        path_classifier = dir_path + "/" + "Joined_AE_classif_weights_" + current_time + "_fold_" + current_fold  + "_.h5"
        joint_classif.autoencoder.save_weights(path_AE, overwrite = True)
        joint_classif.classifier.save_weights(path_classifier, overwrite = True)

    #### Evaluating model ####

    # train data
    encoded_space_train = joint_classif.encoder.predict(X_train, batch_size=BATCH_ARG)
    reconstruction_train = joint_classif.decoder.predict(encoded_space_train, batch_size=BATCH_ARG)
    classes_train = joint_classif.classifier.predict(encoded_space_train, batch_size=BATCH_ARG)

    # reconstruction loss train
    train_reconstuction_loss = recon_loss_fn(X_train, reconstruction_train)

    # classification loss train
    train_classification_loss = classif_loss_fn(Y_train, classes_train)
    
    if CLASS_LOSS_ARG in ['hinge','squared_hinge']:
        # We used tanh in the output - > return back to sigmoid
        print("Tanh-output")
        print(classes_train[0:5])
        classes_train = helper_functions.sigmoid_from_tanh(classes_train)
        print("Probabilities")
        print(classes_train[0:5])
    
    # accuracy train
    classes_train_rounded = numpy.round(classes_train, decimals=0)
    accuracy_train = accuracy_score(Y_train, classes_train_rounded)
    
    # Auc train
    auc_train = roc_auc_score(Y_train, classes_train)

    # plot confusion
    fig_path_train = dir_path + "/" + "fold_" + current_fold + "_conf_matrix_train.pdf"
    helper_functions.make_confusion_matrix(y_true=Y_train,
                                    y_pred=classes_train_rounded,
                                    savefig=True,
                                    figname=fig_path_train)
    
    # plot reconstruction
    fig_path_recon_train = dir_path + "/" + "fold_" + current_fold + "_AE_recon_train.pdf"
    helper_functions.test_AE_performance(
        real_data=X_train,
        predicted_data=reconstruction_train,
        plot_title="Autoencoding train performance",
        plot_file= fig_path_recon_train, 
        color_hex="#e52b50")

    # cleaning
    tf.keras.backend.clear_session()
    gc.collect()

    # test data
    encoded_space_test = joint_classif.encoder.predict(X_test, batch_size=BATCH_ARG)
    reconstruction_test = joint_classif.decoder.predict(encoded_space_test, batch_size=BATCH_ARG)
    classes_test = joint_classif.classifier.predict(encoded_space_test, batch_size=BATCH_ARG)

    # reconstruction loss test
    test_reconstuction_loss = recon_loss_fn(X_test, reconstruction_test)

    # classification loss test
    test_classification_loss = classif_loss_fn(Y_test, classes_test)

    if CLASS_LOSS_ARG in ['hinge','squared_hinge']:
        # We used tanh in the output - > return back to sigmoid
        print("Tanh-output")
        print(classes_test[0:5])
        classes_test = helper_functions.sigmoid_from_tanh(classes_test)
        print("Probabilities")
        print(classes_test[0:5])

    # accuracy test
    classes_test_rounded = numpy.round(classes_test, decimals=0)
    accuracy_test = accuracy_score(Y_test, classes_test_rounded)

    # Auc test
    auc_test = roc_auc_score(Y_test, classes_test)

    # plot confusion
    fig_path_test  = dir_path + "/" + "fold_" + current_fold + "_conf_matrix_test.pdf"
    helper_functions.make_confusion_matrix(y_true=Y_test,
                                        y_pred=classes_test_rounded,
                                        savefig=True,
                                        figname=fig_path_test)
    
    # plot reconstruction
    fig_path_recon_test = dir_path + "/" + "fold_" + current_fold + "_AE_recon_test.pdf"
    helper_functions.test_AE_performance(
        real_data=X_test,
        predicted_data=reconstruction_test,
        plot_title="Autoencoding test performance",
        plot_file= fig_path_recon_test, 
        color_hex="#e52b50")
    
    
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
    del joint_classif
    tf.keras.backend.clear_session()
    gc.collect()

    # saving final evaluation
    final_performance = {
        "model": MODEL_ARG,
        "epochs": EPOCH_ARG,
        "fold": current_fold,
        "total_folds": FOLD_ARG,
        "batch_size": BATCH_ARG,
        "mval": VALUES_ARG,
        "balance_alpha": ALPHA_ARG,
        "cpg_type":CONSIST_ARG,
        "limma_features":LIMMA_FEATURE_ARG,
        "outer_cls_activ": OUTER_CLF_ARG,
        "outer_ae_activ": OUTER_AE_ARG,
        "inner_act": INNER_ARG,
        "learning_rate": LR_ARG,
        "drop_rate":DROP_RATE_ARG,
        "classif_loss_fun": CLASS_LOSS_ARG,
        "recon_loss_fun": RECON_LOSS_ARG,
        "batch_normal": BN_ARG,
        "kernel_regular": KR_ARG,
        "train_recon_loss":  float(train_reconstuction_loss),
        "test_recon_loss":  float(test_reconstuction_loss),
        "train_classification_loss":  float(train_classification_loss),
        "test_classification_loss":  float(test_classification_loss),
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
        run_TF_model(meth_tensor[train],
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
    run_TF_model(meth_tensor_train,
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

    output_stats_base_info = data_fr.drop(["train_classification_loss",
                                        "test_classification_loss",
                                        "train_recon_loss",
                                        "test_recon_loss",
                                        "accuracy_train",
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
