import os
import re
import pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold
import shutil
import subprocess
import time

directory = "/home/peng.zhou-umw/project/Biomarker_Detection/Biomarker_Detection_Paper/Compare_Datasets_20230328/GSE54456_PC_Whole_revised"

dataset_name = "GSE54456_PC"

# src_path = "/home/peng.zhou-umw/project/Biomarker_Detection/Biomarker_Detection_Paper/Compare_Datasets_20230328/scripts/BMGD" + "/src"

model_list = ["svm", "knn", "logReg", "bayes"] #svm knn logReg decisionTree randomForest bayes xgboost stacking

gene_type_list = ["protein_coding"] # ["microarray", "protein_coding", "lincRNA", "CircRNA"]

# filter_info_table_name = "2_gene_filter_thresholds.csv"
outlier_remove = True

pvalue_type_list = ["padj"]
pvalue_threshold_list = [0.05]
foldchange_threshold_list = [1.0]
maxmin_remove_threshold_list = [1.0]
overlap_area_threshold_list = [0.01]






# partition_num = 10

PC_maxmin_remove_threshold = str(1.0)
linc_maxmin_remove_threshold = str(0.01)
Circ_maxmin_remove_threshold = str(0.0)
microarray_maxmin_remove_threshold = str(0.0)

#========================

Target_info = ["Balanced Accuracy 10CV"]
High_Corr_Remove_info = [1.0]# 0.8, 0.9,
curve_point_info = ["knee"] # ["highest", "knee"]
Model_Selection_Methods_info = ["Balanced Accuracy 10CV"] #["Balanced Accuracy 10CV", "F1 10CV", "AUC 10CV"]
pipeline_options = ["General"] #["General", "Fresh"]
pipeline_path = "./src/pipeline_files/ML_pipeline_CV_20220728_SHAP.py"

pseudo_foldchange = str(0.000001) # 1e-6


High_Corr_Remove_type_info = ["New"]#["Original", "New"]


#========================



# gene_filter_info_df = pd.read_csv(directory + "/" + filter_info_table_name)
#
# pvalue_type_list = gene_filter_info_df["pvalue_type"].tolist()
# pvalue_threshold_list = gene_filter_info_df["pvalue_threshold"].tolist()
# foldchange_threshold_list = gene_filter_info_df["foldchange_threshold"].tolist()
# # overlap_area_threshold_list = gene_filter_info_df["overlap_area_threshold"].tolist()
# # maxmin_remove_threshold_list = gene_filter_info_df["maxmin_remove_threshold"].tolist()

# if outlier_remove:
#     pvalue_type_list = pvalue_type_list * 2
#     pvalue_threshold_list = pvalue_threshold_list * 2
#     foldchange_threshold_list = foldchange_threshold_list * 2

folder_list = []
dataset_name_list = []

# general_output_folder = directory + "/" + dataset_name
#
# folder_list.append(general_output_folder)
#
# dataset_name_list.append(dataset_name)
#
#
# if outlier_remove:
#     outlier_remove_output_folder = directory + "/" + dataset_name + "_outlier_remove"
#
#     folder_list.append(outlier_remove_output_folder)
#
#     dataset_name_list.append(dataset_name + "_outlier_remove")

folder_list = [folder_path]
if outlier_remove:
    dataset_name_list = [dataset_name + "_outlier_remove"]
else:
    dataset_name_list = [dataset_name]


from_separation = True










# pvalue_type = "padj"
# pvalue_threshold = str(0.05)
# foldchange_threshold = str(0.585)

# pvalue_type = "padj" ##"pvalue"
# pvalue_type_dict = {}
# pvalue_type_dict["protein_coding"] = "padj"
# pvalue_type_dict["lincRNA"] = "padj"
# pvalue_type_dict["microarray"] = "padj"
# pvalue_type_dict["Circ"] = "padj"
#
#
# # pvalue_threshold = str(0.05)
# pvalue_threshold_dict = {}
# pvalue_threshold_dict["protein_coding"] = str(0.05)
# pvalue_threshold_dict["lincRNA"] = str(0.05)
# pvalue_threshold_dict["microarray"] = str(0.05)
# pvalue_threshold_dict["Circ"] = str(0.05)
#
# # foldchange_threshold = str(1.0)
# foldchange_threshold_dict = {}
# foldchange_threshold_dict["protein_coding"] = str(1.0)
# foldchange_threshold_dict["lincRNA"] = str(1.0)
# foldchange_threshold_dict["microarray"] = str(1.0)
# foldchange_threshold_dict["Circ"] = str(1.0)



# DE_type = "wilcox" ## "wilcox", "DEseq2", "limma"


# reformat_CircRNA_path = ML_path + '/reformat_CircRNA_input.txt'


### if not use the e
# linc_padj_path = "summary_step1.txt.ReadCount.lincRNA.DE.output"
# PC_padj_path = "summary_step1.txt.ReadCount.ProteinCoding.DE.output"
# linc_FPKM_path = "summary_step1.txt.FPKM.lincRNA"
# PC_FPKM_path = "summary_step1.txt.FPKM.ProteinCoding"

########### need to define & test when usage
# Circ_padj_path = output_directory + "/GSE171936.txt.edgeR.output"
# Circ_FPKM_path = output_directory + "/GSE171936_normalize.DEC.txt"
###################

current_partition = 0
# for current_partition in range(len(pvalue_type_list)):

current_parent_partition_folder_path = folder_list[current_partition]

output_directory = current_parent_partition_folder_path + "/output"

ML_path = output_directory + "/ML"

if not os.path.isdir(ML_path):
    os.mkdir(ML_path)

os.chdir(ML_path)



if from_separation:
    cvf = RepeatedStratifiedKFold(n_splits=10, n_repeats=1, random_state=7)

    # DE_output_folder_path = output_directory + "/" + DE_type + "_DE_output"
    #
    # if os.path.isdir(DE_output_folder_path):
    #     shutil.rmtree(DE_output_folder_path)
    # os.mkdir(DE_output_folder_path)

    ID_list = []
    Date_list = []
    JobID_list = []

    Train_dataset_path_list = []
    Test_dataset_path_list = []
    folder_path_list = []


    Target_list = []
    High_Corr_Remove_list = []
    High_Corr_Remove_type_list = []
    curve_point_list = []
    Dataset_list = []
    pipeline_type_list = []
    Model_Selection_Methods_list = []


    ID_num = 1

    for gene_type in gene_type_list:
        if gene_type == "protein_coding":
            output_name = "ProteinCoding"
            FPKM_name = "summary_step1.txt.FPKM.ProteinCoding"
            ReadCount_name = "summary_step1.txt.ReadCount.ProteinCoding"
        elif gene_type == "lincRNA":
            output_name = "lincRNA"
            FPKM_name = "summary_step1.txt.FPKM.lincRNA"
            ReadCount_name = "summary_step1.txt.ReadCount.lincRNA"
        elif gene_type == "microarray":
            output_name = "microarray"
            FPKM_name = "summary_step1.txt.microarray"
            ReadCount_name = "summary_step1.txt.microarray"

        FPKM_path = output_directory + "/Filtered_" + FPKM_name

        ReadCount_path = output_directory + "/Filtered_" + ReadCount_name

        folder_path = ML_path + "/" + gene_type

        if os.path.isdir(folder_path):
            shutil.rmtree(folder_path)
        os.mkdir(folder_path)

        all_df = pd.read_csv(output_directory + '/Filtered_' + gene_type + '_input.all.txt', sep="\t", index_col=0)

        condition_list = all_df["y"].tolist()
        sample_name_list = all_df.index.tolist()

        ReadCount_df = pd.read_csv(ReadCount_path, sep="\t", index_col=0)
        FPKM_df = pd.read_csv(FPKM_path, sep="\t", index_col=0)


        record_flag = 0
        for train_index, test_index in cvf.split(sample_name_list, condition_list):

            partition_folder_path = folder_path + "/CV_partition_" + str(record_flag)

            record_flag = record_flag + 1

            if os.path.isdir(partition_folder_path):
                shutil.rmtree(partition_folder_path)
            os.mkdir(partition_folder_path)

            os.chdir(partition_folder_path)


            X_train_list = []
            for index in train_index:
                X_train_list.append(sample_name_list[index])

            X_test_list = []
            for index in test_index:
                X_test_list.append(sample_name_list[index])

            FPKM_output_train_path = partition_folder_path + "/" + FPKM_name + ".train"
            FPKM_output_test_path = partition_folder_path + "/" + FPKM_name + ".test"

            FPKM_df[X_train_list].to_csv(FPKM_output_train_path, sep="\t")
            FPKM_df[X_test_list].to_csv(FPKM_output_test_path, sep="\t")

            ReadCount_output_train_path = partition_folder_path + "/" + ReadCount_name + ".train"
            ReadCount_output_test_path = partition_folder_path + "/" + ReadCount_name + ".test"

            ReadCount_df[X_train_list].to_csv(ReadCount_output_train_path, sep="\t")
            ReadCount_df[X_test_list].to_csv(ReadCount_output_test_path, sep="\t")

            # partition_folder_path + '/reformat_' + gene_type + '_input.' + dataset_type + '.txt'

            reformat_output_train_path = partition_folder_path + '/reformat_' + gene_type + '_input.train.txt'
            reformat_output_test_path = partition_folder_path + '/reformat_' + gene_type + '_input.test.txt'

            all_df.T[X_train_list].T.to_csv(reformat_output_train_path, sep="\t")
            all_df.T[X_test_list].T.to_csv(reformat_output_test_path, sep="\t")

            current_date = time.strftime('%Y%m%d', time.localtime(time.time()))
            exact_time = time.strftime('%Y-%m-%d_%H:%M', time.localtime(time.time()))

            for High_Corr_Remove_type in High_Corr_Remove_type_info:
                for Target in Target_info:
                    for Model_Selection_Method_type in Model_Selection_Methods_info:
                        for curve_point_type in curve_point_info:
                            for High_Corr_Remove_threshold in High_Corr_Remove_info:
                                for pipeline_type in pipeline_options:
                                    ID_list.append(current_date + str(ID_num).zfill(2))
                                    ID_num  = ID_num + 1
                                    Date_list.append(exact_time)
                                    JobID_list.append(" ")

                                    Train_dataset_path_list.append(reformat_output_train_path)
                                    Test_dataset_path_list.append(reformat_output_test_path)
                                    folder_path_list.append(partition_folder_path)


                                    Target_list.append(Target)
                                    High_Corr_Remove_list.append(High_Corr_Remove_threshold)
                                    High_Corr_Remove_type_list.append(High_Corr_Remove_type)
                                    curve_point_list.append(curve_point_type)
                                    Dataset_list.append(gene_type + "_" + dataset_name_list[current_partition])
                                    pipeline_type_list.append(pipeline_type)
                                    Model_Selection_Methods_list.append(Model_Selection_Method_type)


    ML_df = pd.DataFrame({"ID": ID_list, "Date": Date_list, "JobID": JobID_list, "Target": Target_list,
                          "Folder_path": folder_path_list, "TrainDataset_path": Train_dataset_path_list, "TestDataset_path": Test_dataset_path_list,
                         "High_Corr_Remove": High_Corr_Remove_list, "High_Corr_Remove_type": High_Corr_Remove_type_list,
                          "curve_point": curve_point_list, "Dataset": Dataset_list, "pipeline_type": pipeline_type_list,
                          "Model_Selection_Method":Model_Selection_Methods_list,

                          })

    ML_df.to_csv(ML_path + '/record.csv', index=False)



ML_df = pd.read_csv(ML_path + '/record.csv')

df = ML_df.copy()

error = 0

for i in range(len(df)):
    Target = df.iloc[i].at['Target']
    if Target == "Test":
        continue
    ID = df.iloc[i].at['ID']
    print(ID)

    pipeline_type_read = df.iloc[i].at['pipeline_type']

    if pipeline_type_read == "General":
        pipeline_type = "general"
    elif pipeline_type_read == "Fresh":
        pipeline_type = "fresh"

    High_Corr_Remove_type = df.iloc[i].at['High_Corr_Remove_type']

    Data_address = df.iloc[i].at['TrainDataset_path']

    Testset_address = df.iloc[i].at['TestDataset_path'] #"False"

    Dataset_type_content = df.iloc[i].at['Dataset']
    folder_path = df.iloc[i].at['Folder_path']

    if Dataset_type_content.split("_")[0] == "protein":
        gene_type = "protein_coding"
    elif Dataset_type_content.split("_")[0] == "lincRNA":
        gene_type = "lincRNA"
    elif Dataset_type_content.split("_")[0] == "CircRNA":
        gene_type = "CircRNA"
    elif Dataset_type_content.split("_")[0] == "microarray":
        gene_type = "microarray"
    else:
        print("Datatype not found")
        error = error + 1



    if High_Corr_Remove_type == "New":
        if gene_type == "protein_coding":
            linc_padj = "False"
            linc_FPKM = "False"
            PC_padj = folder_path + "/summary_step1.txt.ReadCount.ProteinCoding.train.DE.output"
            PC_FPKM = folder_path + "/summary_step1.txt.FPKM.ProteinCoding.train"
            Circ_padj = "False"
            Circ_FPKM = "False"
            microarray_padj = "False"
            microarray_FPKM = "False"
        elif gene_type == "lincRNA":
            linc_padj = folder_path + "/summary_step1.txt.ReadCount.lincRNA.train.DE.output"
            linc_FPKM = folder_path + "/summary_step1.txt.FPKM.lincRNA.train"
            PC_padj = "False"
            PC_FPKM = "False"
            Circ_padj = "False"
            Circ_FPKM = "False"
            microarray_padj = "False"
            microarray_FPKM = "False"
        elif gene_type == "microarray":
            linc_padj = "False"
            linc_FPKM = "False"
            PC_padj = "False"
            PC_FPKM = "False"
            Circ_padj = "False"
            Circ_FPKM = "False"
            microarray_padj = folder_path + "/summary_step1.txt.ReadCount.microarray.train.DE.output"
            microarray_FPKM = folder_path + "/summary_step1.txt.microarray.train"
        elif gene_type == "CircRNA":
            print("CircRNA relative data not defined")
            break

    elif High_Corr_Remove_type == "Original":
        linc_padj = "False"
        linc_FPKM = "False"
        PC_padj = "False"
        PC_FPKM = "False"
        Circ_padj = "False"
        Circ_FPKM = "False"
        microarray_padj = "False"
        microarray_FPKM = "False"

    if gene_type == "protein_coding":
        info_table = output_directory + '/Filtered_protein_coding_Gene_info.csv'
    elif gene_type == "lincRNA":
        info_table = output_directory + '/Filtered_lincRNA_Gene_info.csv'
    elif gene_type == "microarray":
        info_table = output_directory + '/Filtered_microarray_Gene_info.csv'

    Model_Selection_Methods = df.iloc[i].at['Model_Selection_Method']
    Model_Selection_Methods_split_list = Model_Selection_Methods.split()
    Model_Selection_Methods_Metric = Model_Selection_Methods_split_list[0]
    Model_Selection_Methods_Evaluation_Method = Model_Selection_Methods_split_list[1]

    if Model_Selection_Methods_Metric.upper() == "F1":
        Model_Selection_Methods_Metric_trans = "F1"
    elif Model_Selection_Methods_Metric.upper() == "AUC":
        Model_Selection_Methods_Metric_trans = "AUC"
    elif Model_Selection_Methods_Metric.upper() == "ACCURACY":
        Model_Selection_Methods_Metric_trans = "ACC"
    elif Model_Selection_Methods_Metric.upper() == "BALANCED":
        Model_Selection_Methods_Metric_trans = "BACC"
        Model_Selection_Methods_Evaluation_Method = Model_Selection_Methods_split_list[2]
    else:
        print("Model_Selection_Methods_Metric not found")
        error = error + 1

    if Model_Selection_Methods_Evaluation_Method.upper() == "10CV":
        Model_Selection_Methods_Evaluation_Method_trans = "cv"
    elif Model_Selection_Methods_Evaluation_Method.upper() == "100HOLDOUT":
        Model_Selection_Methods_Evaluation_Method_trans = "headout"
    elif Model_Selection_Methods_Evaluation_Method.upper() == "100BOOTSTRAP":
        Model_Selection_Methods_Evaluation_Method_trans = "bootstrap"
    else:
        print("Model_Selection_Methods_Evaluation_Method not found")
        error = error + 1

    High_Corr_Remove = str(df.iloc[i].at['High_Corr_Remove'])

    curve_point = str(df.iloc[i].at['curve_point'])

    save_directory = folder_path + "/" + str(ID)

    if os.path.isdir(save_directory):
        shutil.rmtree(save_directory)
    os.mkdir(save_directory)

    bash_string_0 = '''#! /bin/bash

    #BSUB -L /bin/bash
    #BSUB -q short
    #BSUB -n 20 -W 8:00
    #BSUB -o myjob.out
    #BSUB -e myjob.err
    #BSUB -R span[hosts=1]
    #BSUB -R rusage[mem=5000]
    #BSUB -u Peng.Zhou@umassmed.edu

    source activate /home/peng.zhou-umw/shared_env/share_env_20230216_ML_pipeline_R

    '''

    bash_string_1 = 'python {pipeline_path} --input {Data_address} --info_table {info_table} --model {model_list} \
    --metric {Model_Selection_Methods_Metric} --feature_selection True --testset {Testset_address} \
    --val_method_selection {Model_Selection_Methods_Evaluation_Method} --val_method_final cv headout bootstrap \
    --corrthreshold {High_Corr_Remove} --curve_point {curve_point} --pipeline_type {pipeline_type} --linc_padj {linc_padj} --PC_padj {PC_padj} \
    --linc_FPKM {linc_FPKM} --PC_FPKM {PC_FPKM} --Circ_padj {Circ_padj} --Circ_FPKM {Circ_FPKM} --pseudo_foldchange {pseudo_foldchange} --pvalue_type {pvalue_type} --pvalue_threshold {pvalue_threshold} --foldchange_threshold {foldchange_threshold} \
    --PC_maxmin_remove_threshold {PC_maxmin_remove_threshold} --linc_maxmin_remove_threshold {linc_maxmin_remove_threshold} --Circ_maxmin_remove_threshold {Circ_maxmin_remove_threshold} --final_evaluation False --microarray_padj {microarray_padj} --microarray_FPKM {microarray_FPKM} --microarray_maxmin_remove_threshold {microarray_maxmin_remove_threshold}'


    # bash_string_1 = 'python {pipeline_path} --input {Data_address} --model svm knn logReg decisionTree randomForest bayes xgboost stacking \
    # --metric {Model_Selection_Methods_Metric} --feature_selection True --testset {Testset_address} \
    # --val_method_selection {Model_Selection_Methods_Evaluation_Method} --val_method_final cv headout bootstrap \
    # --corrthreshold {High_Corr_Remove} --curve_point {curve_point} --pipeline_type {pipeline_type} --linc_padj {linc_padj} --PC_padj {PC_padj} \
    # --linc_FPKM {linc_FPKM} --PC_FPKM {PC_FPKM} --Circ_padj {Circ_padj} --Circ_FPKM {Circ_FPKM} --pseudo_foldchange {pseudo_foldchange} --pvalue_type {pvalue_type} --pvalue_threshold {pvalue_threshold} --foldchange_threshold {foldchange_threshold} \
    # --PC_maxmin_remove_threshold {PC_maxmin_remove_threshold} --linc_maxmin_remove_threshold {linc_maxmin_remove_threshold} --Circ_maxmin_remove_threshold {Circ_maxmin_remove_threshold} --final_evaluation False'
    #
    #


    bash_string_1_R = bash_string_1.format(Data_address=Data_address, pipeline_path=pipeline_path, \
                                           info_table=info_table, \
                                           Model_Selection_Methods_Metric=Model_Selection_Methods_Metric_trans, \
                                           Testset_address=Testset_address,
                                           Model_Selection_Methods_Evaluation_Method=Model_Selection_Methods_Evaluation_Method_trans, \
                                           High_Corr_Remove=High_Corr_Remove, curve_point=curve_point,
                                           pipeline_type=pipeline_type,
                                           linc_padj=linc_padj, PC_padj=PC_padj, linc_FPKM=linc_FPKM, PC_FPKM=PC_FPKM,
                                           Circ_padj=Circ_padj,
                                           Circ_FPKM=Circ_FPKM, pseudo_foldchange=pseudo_foldchange,
                                           pvalue_type=pvalue_type_list[current_partition], pvalue_threshold=pvalue_threshold_list[current_partition],
                                           foldchange_threshold=foldchange_threshold_list[current_partition],
                                           PC_maxmin_remove_threshold=PC_maxmin_remove_threshold,
                                           linc_maxmin_remove_threshold=linc_maxmin_remove_threshold,
                                           Circ_maxmin_remove_threshold=Circ_maxmin_remove_threshold,
                                           microarray_padj=microarray_padj, microarray_FPKM=microarray_FPKM, microarray_maxmin_remove_threshold=microarray_maxmin_remove_threshold, model_list= ' '.join(str(e) for e in model_list))

    bash_string = bash_string_0 + bash_string_1_R

    # output = bash_string.replace("^M",'abc')
    # print(output)
    with open(save_directory + '/pipeline_submit.sh', 'w') as rsh:
        rsh.write(bash_string)

    os.chdir(save_directory)
    # os.system('dos2unix ./pipeline_submit.sh')
    system_output = os.popen('bsub < ./pipeline_submit.sh').read()
    # system_output = "Job <1040418> is submitted to queue <large>."
    p1 = re.compile(r'[<](.*?)[>]', re.S)
    JobID = re.findall(p1, system_output)[0]
    df.loc[i:i, 'JobID'] = JobID
    os.chdir(ML_path)

    del Data_address
    del Model_Selection_Methods_Metric_trans
    del Testset_address
    del Model_Selection_Methods_Evaluation_Method_trans
    del High_Corr_Remove
# print("parent partition " + str(current_partition) + " finished")
print("There are " + str(error) + " errors be found.")
df.to_csv(ML_path + '/record.csv', index=False)


