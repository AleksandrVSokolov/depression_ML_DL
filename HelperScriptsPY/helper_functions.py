### A bunch of helper function to use
### Storing them here so they're easily accessible.
import itertools
import zipfile
import os
import datatable
import pandas
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from sklearn.metrics import confusion_matrix, accuracy_score, precision_recall_fscore_support

# A function to make confusion matrix
def make_confusion_matrix(y_true, y_pred, classes=None, figsize=(10, 10), text_size=15, norm=False, savefig=False, figname ="confusion_matrix.png"): 
  """Makes a labelled confusion matrix comparing predictions and ground truth labels.

  If classes is passed, confusion matrix will be labelled, if not, integer class values
  will be used.

  Args:
    y_true: Array of truth labels (must be same shape as y_pred).
    y_pred: Array of predicted labels (must be same shape as y_true).
    classes: Array of class labels (e.g. string form). If `None`, integer labels are used.
    figsize: Size of output figure (default=(10, 10)).
    text_size: Size of output figure text (default=15).
    norm: normalize values or not (default=False).
    savefig: save confusion matrix to file (default=False).
  
  Returns:
    A labelled confusion matrix plot comparing y_true and y_pred.

  Example usage:
    make_confusion_matrix(y_true=test_labels, # ground truth test labels
                          y_pred=y_preds, # predicted labels
                          classes=class_names, # array of class label names
                          figsize=(15, 15),
                          text_size=10)
  """  
  # Create the confustion matrix
  cm = confusion_matrix(y_true, y_pred)
  cm_norm = cm.astype("float") / cm.sum(axis=1)[:, np.newaxis] # normalize it
  n_classes = cm.shape[0] # find the number of classes we're dealing with

  # Plot the figure and make it pretty
  fig, ax = plt.subplots(figsize=figsize)
  cax = ax.matshow(cm, cmap=plt.cm.Blues) # colors will represent how 'correct' a class is, darker == better
  fig.colorbar(cax)

  # Are there a list of classes?
  if classes:
    labels = classes
  else:
    labels = np.arange(cm.shape[0])
  
  # Label the axes
  ax.set(title="Confusion Matrix",
         xlabel="Predicted label",
         ylabel="True label",
         xticks=np.arange(n_classes), # create enough axis slots for each class
         yticks=np.arange(n_classes), 
         xticklabels=labels, # axes will labeled with class names (if they exist) or ints
         yticklabels=labels)
  
  # Make x-axis labels appear on bottom
  ax.xaxis.set_label_position("bottom")
  ax.xaxis.tick_bottom()

  # Set the threshold for different colors
  threshold = (cm.max() + cm.min()) / 2.

  # Plot the text on each cell
  for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
    if norm:
      plt.text(j, i, f"{cm[i, j]} ({cm_norm[i, j]*100:.1f}%)",
              horizontalalignment="center",
              color="white" if cm[i, j] > threshold else "black",
              size=text_size)
    else:
      plt.text(j, i, f"{cm[i, j]}",
              horizontalalignment="center",
              color="white" if cm[i, j] > threshold else "black",
              size=text_size)

  # Save the figure to the current working directory
  if savefig:
    fig.savefig(figname)
    plt.close()

# A function to test autoencoding performance
def test_AE_performance(real_data, predicted_data, plot_title, plot_file, color_hex):
    
    if predicted_data.shape[0] >= 100:
      random_participants = np.random.choice(predicted_data.shape[0], size=100, replace=False, p=None)
    else:
      random_participants = np.random.choice(predicted_data.shape[0], size=predicted_data.shape[0], replace=False, p=None)

    if predicted_data.shape[1] >= 100:
      random_CpGs = np.random.choice(predicted_data.shape[1], size=100, replace=False, p=None)
    else:
      random_CpGs = np.random.choice(predicted_data.shape[1], size=predicted_data.shape[1], replace=False, p=None)

    # Select CpGs
    x = real_data[:,random_CpGs]
    y = predicted_data[:,random_CpGs]

    # Select Participants
    x = x[random_participants,:]
    y = y[random_participants,:]

    # Flatten
    real_data = np.ravel(real_data)
    predicted_data = np.ravel(predicted_data)
    x = np.ravel(x)
    y = np.ravel(y) 

    # Adding text
    total_cor = np.corrcoef(real_data, predicted_data)[0, 1]
    sample_cor = np.corrcoef(x, y)[0, 1]
    txt = "r (total) = " + str("{:.4f}".format(total_cor)) + "\n" + "r (sample) = " + str("{:.4f}".format(sample_cor))

    fig, ax = plt.subplots()
    ax.scatter(x, y, c=color_hex, alpha=0.5, label="Predicted VS Real")
    ax.set_xlabel("Real", weight="heavy")
    ax.set_ylabel("Predicted", weight="heavy")
    ax.legend()
    plt.title(plot_title, weight="heavy")
    fig.text(.5, -.1, txt, ha='center')
    print(f"Correlation: {total_cor}")

    # Save figure
    fig.savefig(fname=plot_file, format="pdf", dpi=300., bbox_inches="tight")

# Create function to unzip a zipfile into current working directory
def unzip_data(filename):
  """
  Unzips filename into the current working directory.

  Args:
    filename (str): a filepath to a target zip folder to be unzipped.
  """
  zip_ref = zipfile.ZipFile(filename, "r")
  zip_ref.extractall()
  zip_ref.close()

# Walk through a directory and find out how many files are in each subdirectory.
def walk_through_dir(dir_path):
  """
  Walks through dir_path returning its contents.

  Args:
    dir_path (str): target directory
  
  Returns:
    A print out of:
      number of subdiretories in dir_path
      number of images (files) in each subdirectory
      name of each subdirectory
  """
  for dirpath, dirnames, filenames in os.walk(dir_path):
    print(f"There are {len(dirnames)} directories and {len(filenames)} images in '{dirpath}'.")
    
# Preprocessing data and import
def load_DNA_methylation_data(select_cpg="ALL", 
                              cohort_pool=["GSE113725_RDE",
                              "GSE125105_MPIP",
                              "GSE198904_DHRC",
                              "GSE198904_OBS",
                              "GSE72680_GRADY",
                              "GSE74414_MPIP2",
                              "PSY_RC",
                              "PSY_SCR"],
                              cpg_folder="Main_analysis/ML_data/"):

    #### Loading data ####
    # Methylation
    print("Reading data")
    CV_mval = datatable.fread(cpg_folder + "CV_mval.csv")
    CV_mval = CV_mval.to_pandas()
    CV_mval.set_index('C0', inplace=True)
    test_mval = datatable.fread(cpg_folder + "test_mval.csv")
    test_mval = test_mval.to_pandas()
    test_mval.set_index('C0', inplace=True)

    # Transposing
    CV_mval = CV_mval.transpose()
    test_mval = test_mval.transpose()

    # Phenotypes
    print("Reading Phenotypes")
    CV_pheno = pandas.read_csv(cpg_folder + "CV_pheno.csv")
    test_pheno =  pandas.read_csv(cpg_folder + "test_pheno.csv")

    # Subsetting data based on cohorts
    CV_pheno = CV_pheno[CV_pheno['Study'].isin(cohort_pool)]
    test_pheno = test_pheno[test_pheno['Study'].isin(cohort_pool)]
    CV_mval = CV_mval[CV_mval.index.isin(CV_pheno["ID"])]
    test_mval = test_mval[test_mval.index.isin(test_pheno["ID"])]


    # Annotation
    print("Reading Annotation")
    ordered_illum_annot = datatable.fread(cpg_folder + "ILLUM_450K_ANNOT_selected_ML.csv")
    ordered_illum_annot = ordered_illum_annot.to_pandas()
    ordered_illum_annot.set_index('C0', inplace=True)

    if select_cpg != "ALL":
      # CpGs
      print("Reading CpG list and selecting CpGs")
      CpG_list_path = cpg_folder + select_cpg
      CpG_list = pandas.read_csv(CpG_list_path, header=None)[0].tolist()
      # Selecting important features
      CV_mval = CV_mval.iloc[:,CV_mval.columns.isin(CpG_list)]
      test_mval = test_mval.iloc[:,test_mval.columns.isin(CpG_list)]
    else:
      print("All CpGs will be kept")

    CpG_names = CV_mval.columns.values
    CpG_names = CpG_names.tolist()

    # Preparing Phenotypes
    print("Preparing phenotypes")
    CV_pheno["depr_one_hot"] = [1. if x == "Case" else 0. for x in CV_pheno["Depression"]]
    test_pheno["depr_one_hot"] = [1. if x == "Case" else 0. for x in test_pheno["Depression"]]

    # Obtaining tensors
    depr_status_tensor_CV = CV_pheno["depr_one_hot"].to_numpy(dtype="float32")
    depr_status_tensor_test = test_pheno["depr_one_hot"].to_numpy(dtype="float32")

    CV_mval_tensor = CV_mval.to_numpy(dtype="float32")
    test_mval_tensor = test_mval.to_numpy(dtype="float32")

    # Preparing betas
    print("Preparing beta values")
    def m2Beta(m_val:float):
        beta = 2.**m_val/(2.**m_val + 1.)
        return beta
    CV_beta_tensor = m2Beta(CV_mval_tensor)
    test_beta_tensor = m2Beta(test_mval_tensor)

    output_dict = {"CV_mval": CV_mval,
               "test_mval": test_mval, 
               "CV_pheno": CV_pheno, 
               "test_pheno": test_pheno, 
               "CV_mval_tensor": CV_mval_tensor,
               "test_mval_tensor":test_mval_tensor,
               "CV_beta_tensor":CV_beta_tensor,
               "test_beta_tensor":test_beta_tensor,
               "depr_status_tensor_CV":depr_status_tensor_CV,
               "depr_status_tensor_test":depr_status_tensor_test,
               "ordered_illum_annot":ordered_illum_annot,
               "CpG_names":CpG_names}
    
    return output_dict

# Make barplots of counts vs categories
def make_stats_bars(data, category_main, sub_category, plot_path):
    data = data.groupby(category_main).value_counts([sub_category])
    data = data.reset_index()
    data = data.astype({category_main: 'str'})
    data = pandas.DataFrame(data)
    fig,ax = plt.subplots()
    ax = sns.barplot(data=data, x=category_main, y='count', hue=sub_category, palette=sns.color_palette("Set1"))
    fig.savefig(plot_path)
    plt.close()

    # Saving dataset
    file_csv = plot_path.split(".")[0]
    file_csv = file_csv + ".csv"
    data.to_csv(file_csv)

# Drop all unnamed columns
def drop_useless_column(pd_df):
    pd_df = pd_df.drop(pd_df.filter(like='Unnamed:').columns, axis=1)
    return pd_df

# Recalculate Sigmoid function (probabibilities) from tanh
def sigmoid_from_tanh(tanh):
        
    if tanh == -1.:
        return 0.
    elif tanh == 1.:
        return 1.
    
    x = np.log((1+tanh)/(1-tanh))/2
    sigmoid = 1/(1 + np.exp(-x))
    return sigmoid

sigmoid_from_tanh = np.vectorize(sigmoid_from_tanh)

def explain_predictions_folds(dir_path, x):
  
  # Reading files
  current_pheno_train = pandas.read_csv(dir_path + "/" + str(x) + "_" + "pheno_train.csv")
  current_pheno_test = pandas.read_csv(dir_path + "/" + str(x) + "_" + "pheno_test.csv")
  current_predictions_train = pandas.read_csv(dir_path + "/" + str(x) + "_" + "classes_train_rounded.csv")
  current_predictions_test = pandas.read_csv(dir_path + "/" + str(x) + "_" + "classes_test_rounded.csv")
  true_labs_train = pandas.read_csv(dir_path + "/" + str(x) + "_" + "Y_train.csv")
  true_labs_test = pandas.read_csv(dir_path + "/" + str(x) + "_" + "Y_test.csv")

  # Fixing first column
  current_predictions_train = drop_useless_column(current_predictions_train)
  current_predictions_test = drop_useless_column(current_predictions_test)
  true_labs_train = drop_useless_column(true_labs_train)
  true_labs_test = drop_useless_column(true_labs_test)
  current_pheno_train = drop_useless_column(current_pheno_train)
  current_pheno_test = drop_useless_column(current_pheno_test)

  # Renaming columns
  current_predictions_train = current_predictions_train.set_axis(["predicted_train"], axis=1)
  current_predictions_test = current_predictions_test.set_axis(["predicted_test"], axis=1)
  true_labs_train = true_labs_train.set_axis(["true_labs_train"], axis=1)
  true_labs_test = true_labs_test.set_axis(["true_labs_test"], axis=1)

  # Merging data
  tmp_df_train = pandas.concat([current_pheno_train, current_predictions_train, true_labs_train], 
                              axis=1, ignore_index=True)
  tmp_df_test = pandas.concat([current_pheno_test, current_predictions_test, true_labs_test], 
                              axis=1, ignore_index=True)
  tmp_df_train = tmp_df_train.set_axis(["ID","Sex","Age","Depression","Study","depr_one_hot",
                                      "predicted_train","true_labs_train"], axis=1)
  tmp_df_test = tmp_df_test.set_axis(["ID","Sex","Age","Depression","Study","depr_one_hot",
                                      "predicted_test","true_labs_test"], axis=1)
  
  # Plotting data
  fold_str = "fold_" + str(x)

  make_stats_bars(data=tmp_df_train, 
                  category_main="predicted_train",
                  sub_category="Study",
                  plot_path=dir_path + "/" + fold_str + "_predictions_study_train.pdf")

  make_stats_bars(data=tmp_df_train, 
                  category_main="predicted_train",
                  sub_category="Sex",
                  plot_path=dir_path + "/" + fold_str + "_predictions_gender_train.pdf")

  make_stats_bars(data=tmp_df_test, 
                  category_main="predicted_test",
                  sub_category="Study",
                  plot_path=dir_path + "/" + fold_str + "_predictions_study_test.pdf")

  make_stats_bars(data=tmp_df_test, 
                  category_main="predicted_test",
                  sub_category="Sex",
                  plot_path=dir_path + "/" + fold_str + "_predictions_gender_test.pdf")
  
  # Statisticts for false predictions
  tmp_df_train["Mismatch"] = ["False_pred" if tmp_df_train["predicted_train"].values[x] != tmp_df_train["true_labs_train"].values[x] else "True_pred"
                              for x in range(tmp_df_train.shape[0])]
  tmp_df_test["Mismatch"] = ["False_pred" if tmp_df_test["predicted_test"].values[x] != tmp_df_test["true_labs_test"].values[x] else "True_pred"
                             for x in range(tmp_df_test.shape[0])]
  
  # Write final datasets
  tmp_df_train.to_csv(dir_path + "/" + fold_str + "_final_dataset_train.csv")
  tmp_df_test.to_csv(dir_path + "/" + fold_str + "_final_dataset_test.csv")

  # Plot distributions of false predictions across cohorts and sex
  make_stats_bars(data=tmp_df_train, 
                  category_main="Mismatch",
                  sub_category="Study",
                  plot_path=dir_path + "/" + fold_str + "_false_predictions_study_train.pdf")

  make_stats_bars(data=tmp_df_train, 
                  category_main="Mismatch",
                  sub_category="Sex",
                  plot_path=dir_path + "/" + fold_str + "_false_predictions_gender_train.pdf")

  make_stats_bars(data=tmp_df_test, 
                  category_main="Mismatch",
                  sub_category="Study",
                  plot_path=dir_path + "/" + fold_str + "_false_predictions_study_test.pdf")

  make_stats_bars(data=tmp_df_test, 
                  category_main="Mismatch",
                  sub_category="Sex",
                  plot_path=dir_path + "/" + fold_str + "_false_predictions_gender_test.pdf")  
  
def plot_roc_curve(fpr, tpr, auc, path):
    plt.figure()
    plt.plot([0, 1], [0, 1], 'k--')
    plt.plot(fpr, tpr, label='Classifier (area = {:.3f})'.format(auc))
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve')
    plt.legend(loc='best')
    plt.savefig(path, dpi=300)
    plt.close()

def select_features_R_limma(train_methyl,
                            train_pheno,
                            CpG_names,
                            current_fold,
                            base_dir,
                            command = "Rscript",
                            path2script ="/HelperScriptsPY/RunLimmaMatrix.R",
                            feature_number = 50):

  # Writing Limma inputs                      
  pheno_path = base_dir + "/pheno_" + current_fold + ".csv"
  meth_path = base_dir + "/meth_" + current_fold + ".csv"
  toptable_path = base_dir + "/TopTable_" + current_fold + ".csv"

  # Saving things 
  train_pheno.to_csv(pheno_path)
  train_methyl = datatable.Frame(train_methyl)
  train_methyl.names = CpG_names
  train_methyl.to_csv(meth_path)

  # Running Limma inside R
  args = [base_dir, pheno_path, meth_path, toptable_path]
  cmd = [command, path2script] + args
  subprocess.run(cmd)

  # Removing junk
  os.remove(pheno_path)
  os.remove(meth_path)

  # Selecting features
  toptable = pandas.read_csv(toptable_path)
  top_features = toptable["CpG"].values[range(feature_number)]
     

  selection_idx = train_methyl.names
  selection_idx = [True if x in top_features else False for x in selection_idx]
  return selection_idx
  