import pandas as pd
import math
import os
import sys
import pandas as pd
from collections import Counter
from itertools import combinations
import numpy as np
import subprocess
from pdf2image import convert_from_path
import ast
import threading
import time
import multiprocessing as mp



# report_folder_path = "/nl/umw_chan_zhou/Billy/Biomarker_Detection_Paper/TCGA_Breast_Cancer/output/ML/Report_1234"
# report_path = "/nl/umw_chan_zhou/Billy/Biomarker_Detection_Paper/TCGA_Breast_Cancer/output/ML/Report_1234/Job_status.csv"
# Target_Metric = "Balanced_Accuracy"  ## "Accuracy", "Balanced_Accuracy", "ROC_AUC", "Precision", "Recall", "F1"
# Target_range = 0.9
## if it in[0, 1], it is the percentage of the best metric
## if it in(1, 100], it is the specific number of the best metric

# ML_path = "/nl/umw_chan_zhou/Billy/Biomarker_Detection_Paper/TCGA_Breast_Cancer/output/ML"
# DE_type = "wilcox"
# CV_partition = 10




report_folder_path = sys.argv[1]
report_path = sys.argv[2]
Target_Metric = sys.argv[3]
Target_range = float(sys.argv[4])
directory = sys.argv[5]
DE_type = sys.argv[6]
CV_partition = int(sys.argv[7])
summary_type = sys.argv[8]
model_list = ast.literal_eval(sys.argv[9])
gene_type_list = ast.literal_eval(sys.argv[10])
DE_type = sys.argv[6]
intersection_check = int(sys.argv[11])
src_path = sys.argv[12]

# gene_type_list = ["protein_coding", "lincRNA"]

best_model_accepted_treshold = 0.95


CV_performance_evaluation = True
result_organization = True


output_directory = directory + "/output"
ML_path = output_directory + "/ML"


Heatmap_info_dict = {}

Heatmap_info_dict["protein_coding"] = {}
Heatmap_info_dict["protein_coding"]["padj"] = output_directory + "/summary_step1.txt.ReadCount.ProteinCoding.DE.output"
Heatmap_info_dict["protein_coding"]["FPKM"] = output_directory + "/summary_step1.txt.FPKM.ProteinCoding"
Heatmap_info_dict["lincRNA"] = {}
Heatmap_info_dict["lincRNA"]["padj"] = output_directory + "/summary_step1.txt.ReadCount.lincRNA.DE.output"
Heatmap_info_dict["lincRNA"]["FPKM"] = output_directory + "/summary_step1.txt.FPKM.lincRNA"
Heatmap_info_dict["microarray"] = {}
Heatmap_info_dict["microarray"]["padj"] = output_directory + "/summary_step1.txt.microarray.DE.output"
Heatmap_info_dict["microarray"]["FPKM"] = output_directory + "/summary_step1.txt.microarray"
Heatmap_info_dict["combine"] = {}
Heatmap_info_dict["combine"]["padj"] = output_directory + "/summary_step1.txt.ReadCount.Combination.DE.output"
Heatmap_info_dict["combine"]["FPKM"] = output_directory + "/summary_step1.txt.FPKM.Combination"

# PC_padj = PC_padj_path
# PC_FPKM = PC_FPKM_path
# linc_padj = linc_padj_path
# linc_FPKM = linc_FPKM_path
# Circ_padj = Circ_padj_path
# Circ_FPKM = Circ_FPKM_path




#model_list = ["svm", "knn", "logReg", "decisionTree", "randomForest",  "bayes", "xgboost"] #["svm", "knn", "logReg", "decisionTree", "randomForest", "bayes", "xgboost", "stacking"]
record_metric_list = ['feature_num', 'effective_feature_num', 'Accuracy', 'Balanced_Accuracy', 'ROC_AUC', 'Precision', 'Recall', 'Specificity', 'F1', 'TN', 'FP', 'FN', 'TP']


report_df_ori = pd.read_csv(report_path)


check_name_list = report_df_ori["Dataset"].tolist()
for name in check_name_list:
    if name.split("_")[0] == "protein":
        protein_dataset = name
    elif name.split("_")[0] == "lincRNA":
        lincRNA_dataset = name
    elif name.split("_")[0] == "microarray":
        microarray_dataset = name
    elif name.split("_")[0] == "combine":
        combine_dataset = name


report_df_ori = report_df_ori[(report_df_ori["Status"] == "Finished")]
report_df_ori = report_df_ori[(report_df_ori["curve_point"] == "knee")]

report_df_ori.to_csv(report_path+".tem", index = False)
report_df_ori = pd.read_csv(report_path+".tem")

# CV performance evaluation
if CV_performance_evaluation:

    CV_performance_calculation_dict = {}

    for gene_type in gene_type_list:
        if gene_type == "protein_coding":
            report_df = report_df_ori[(report_df_ori["Dataset"] == protein_dataset)]
        elif gene_type == "lincRNA":
            report_df = report_df_ori[(report_df_ori["Dataset"] == lincRNA_dataset)]
        elif gene_type == "microarray":
            report_df = report_df_ori[(report_df_ori["Dataset"] == microarray_dataset)]
        elif gene_type == "combine":
            report_df = report_df_ori[(report_df_ori["Dataset"] == combine_dataset)]

        if len(report_df) == 0:
            continue

        CV_performance_calculation_dict[gene_type] = {}

        tem_column_name_list = report_df.columns.values.tolist()

        for model in model_list:
            if model + '_Accuracy' in tem_column_name_list:
                CV_performance_calculation_dict[gene_type][model] = {}

                if summary_type == "general":
                    CV_performance_calculation_dict[gene_type][model]["ID"] = []
                    for column_format_name in record_metric_list:
                        CV_performance_calculation_dict[gene_type][model][column_format_name] = []
                elif summary_type == "fresh":
                    CV_performance_calculation_dict[gene_type][model]['effective_feature_num'] = []
                    CV_performance_calculation_dict[gene_type][model][Target_Metric] = []

        for i in range(CV_partition):
            current_CV_partition_path = ML_path + "/" + gene_type + "/CV_partition_" + str(i)

            report_df_current_CV_partition = report_df[(report_df["Folder_path"] == current_CV_partition_path)]

            testset_path = report_df_current_CV_partition["TestDataset_path"].tolist()[0]
            testset_size = len(pd.read_csv(testset_path, index_col=0, sep="\t"))

            for model in model_list:
                if model + '_Accuracy' in tem_column_name_list:
                    target_metric_list_bestmodel_metric = report_df_current_CV_partition[model + "_" + "Accuracy"].tolist()
                    best_metric = float(max(target_metric_list_bestmodel_metric))
                    best_model_accepted_metric = best_model_accepted_treshold * best_metric


                    accepted_min_correct_case = int(best_model_accepted_metric * testset_size)

                    best_model_accepted_metric = (accepted_min_correct_case / testset_size) - 0.0001



                    acceptable_model_df = report_df_current_CV_partition[(report_df_current_CV_partition[model + "_" + "Accuracy"] >= best_model_accepted_metric)]
                    sorted_acceptable_model_df = acceptable_model_df.sort_values([model + "_feature_num", model + "_" + Target_Metric], ascending=[True, False])



                    if summary_type == "general":
                        CV_performance_calculation_dict[gene_type][model]["ID"].append(str(sorted_acceptable_model_df.iloc[0].at['ID']))
                        for column_format_name in record_metric_list:
                            CV_performance_calculation_dict[gene_type][model][column_format_name].append(sorted_acceptable_model_df.iloc[0].at[model + "_" + column_format_name])
                    elif summary_type == "fresh":
                        CV_performance_calculation_dict[gene_type][model]['effective_feature_num'].append(sorted_acceptable_model_df.iloc[0].at[model + "_" + 'effective_feature_num'])
                        CV_performance_calculation_dict[gene_type][model][Target_Metric].append(sorted_acceptable_model_df.iloc[0].at[model + "_" + Target_Metric])



    for gene_type in gene_type_list:
        CV_performance_output_dict = {}
        partition_list = []
        for i in range(CV_partition):
            partition_list.append("CV_partition_" + str(i))
        CV_performance_output_dict["category"] = partition_list
        CV_performance_mean_list = []
        CV_performance_std_list = []
        for model in model_list:
            if summary_type == "general":
                CV_performance_output_dict[model + "_ID"] = CV_performance_calculation_dict[gene_type][model]["ID"]
                CV_performance_mean_list.append("/")
                CV_performance_std_list.append("/")
                for column_format_name in record_metric_list:
                    CV_performance_output_dict[model+ "_" + column_format_name] = CV_performance_calculation_dict[gene_type][model][column_format_name]
                    mean_value = np.array(CV_performance_calculation_dict[gene_type][model][column_format_name]).mean()
                    std_value = np.array(CV_performance_calculation_dict[gene_type][model][column_format_name]).std()
                    if column_format_name == 'effective_feature_num' or column_format_name == 'feature_num':
                        CV_performance_mean_list.append(round(mean_value, 2))
                        CV_performance_std_list.append(round(std_value, 2))
                    else:
                        CV_performance_mean_list.append(round(mean_value, 4))
                        CV_performance_std_list.append(round(std_value, 2))
            elif summary_type == "fresh":
                CV_performance_output_dict[model + "_" + 'effective_feature_num'] = CV_performance_calculation_dict[gene_type][model]['effective_feature_num']
                mean_value = np.array(CV_performance_calculation_dict[gene_type][model]['effective_feature_num']).mean()
                std_value = np.array(CV_performance_calculation_dict[gene_type][model]['effective_feature_num']).std()
                CV_performance_mean_list.append(round(mean_value, 2))
                CV_performance_std_list.append(round(std_value, 2))

                CV_performance_output_dict[model + "_" + Target_Metric] = CV_performance_calculation_dict[gene_type][model][Target_Metric]
                mean_value = np.array(CV_performance_calculation_dict[gene_type][model][Target_Metric]).mean()
                std_value = np.array(CV_performance_calculation_dict[gene_type][model][Target_Metric]).std()
                CV_performance_mean_list.append(round(mean_value, 4))
                CV_performance_std_list.append(round(std_value, 4))






        # for key, value in CV_performance_output_dict.items():
        #     print(key, len(value))
        # print(CV_performance_output_dict)
        output_df = pd.DataFrame(CV_performance_output_dict).set_index("category", inplace=False).T

        # print(output_df)
        # print(CV_performance_std_list)
        output_df.insert(loc=0, column='std', value=CV_performance_std_list)
        output_df.insert(loc=0, column='mean', value=CV_performance_mean_list)
        output_path = report_folder_path + "/" + gene_type + "_CV_performance.csv"
        output_df.index.name = "name"
        output_df.to_csv(output_path, index = True)






# filter the identity models

if result_organization:
    for gene_type in gene_type_list:
        if gene_type == "protein_coding":
            report_df = report_df_ori[(report_df_ori["Dataset"] == protein_dataset)]
        elif gene_type == "lincRNA":
            report_df = report_df_ori[(report_df_ori["Dataset"] == lincRNA_dataset)]
        elif gene_type == "microarray":
            report_df = report_df_ori[(report_df_ori["Dataset"] == microarray_dataset)]
        elif gene_type == "combine":
            report_df = report_df_ori[(report_df_ori["Dataset"] == combine_dataset)]

        if len(report_df) == 0:
            continue

        record_dict = {}
        record_dict["name"] = []
        record_dict["ID"] = []
        record_dict["partition"] = []
        record_dict["model"] = []
        record_dict["feature_num"] = []
        record_dict["actual_feature_num"] = []
        record_dict["Accuracy"] = []
        record_dict["Balanced_Accuracy"] = []
        record_dict["ROC_AUC"] = []
        record_dict["Precision"] = []
        record_dict["Recall"] = []
        record_dict['Specificity'] = []
        record_dict["F1"] = []
        record_dict['TN'] = []
        record_dict['FP'] = []
        record_dict['FN'] = []
        record_dict['TP'] = []
        record_dict["features"] = []



        for part_num in range(CV_partition):
            current_CV_partition_path = ML_path + "/" + gene_type + "/CV_partition_" + str(part_num)

            report_df_current_CV_partition = report_df[(report_df["Folder_path"] == current_CV_partition_path)]

            if len(report_df_current_CV_partition) == 0:
                continue

            tem_column_name_list = report_df_current_CV_partition.columns.values.tolist()

            target_metric_list_bestmodel_all = []

            for model in model_list:
                if model + '_Accuracy' in tem_column_name_list:
                    target_metric_list_bestmodel_tem = report_df_current_CV_partition[model + "_" + "Accuracy"].tolist()
                    target_metric_list_bestmodel_all.extend(target_metric_list_bestmodel_tem)

            best_metric = float(max(target_metric_list_bestmodel_all))

            if Target_range >= 0 and Target_range <= 1:
                accepted_threshold = best_metric * Target_range
            elif Target_range > 1 and Target_range <= 100:
                accepted_threshold = Target_range / 100

            testset_path = report_df_current_CV_partition["TestDataset_path"].tolist()[0]
            testset_size = len(pd.read_csv(testset_path,index_col=0, sep="\t"))
            accepted_min_correct_case = int(accepted_threshold * testset_size)

            accepted_threshold = (accepted_min_correct_case/testset_size) - 0.0001

            print(current_CV_partition_path)
            print(accepted_threshold)




            for model in model_list:
                if model + '_Accuracy' in tem_column_name_list:
                    # target_metric_list_bestmodel = report_df_current_CV_partition[model + "_" + Target_Metric].tolist()
                    # best_metric = float(max(target_metric_list_bestmodel))
                    #
                    # if Target_range >= 0 and Target_range <= 1:
                    #     accepted_threshold = best_metric * Target_range
                    # elif Target_range > 1 and Target_range <= 100:
                    #     accepted_threshold = Target_range / 100

                    acceptable_model_df = report_df_current_CV_partition[(report_df_current_CV_partition[model + "_" + "Accuracy"] >= accepted_threshold)]

                    tem_dict = {}
                    tem_dict["features"] = []
                    # file_name = model + "_imp_features.txt"
                    tem_dict["ID"] = acceptable_model_df["ID"].tolist()
                    tem_dict["feature_num"] = acceptable_model_df[model + "_" + "feature_num"].tolist()
                    tem_dict["effective_feature_num"] = acceptable_model_df[model + "_" + "effective_feature_num"].tolist()
                    tem_dict["Accuracy"] = acceptable_model_df[model + "_" + "Accuracy"].tolist()
                    tem_dict["Balanced_Accuracy"] = acceptable_model_df[model + "_" + "Balanced_Accuracy"].tolist()
                    tem_dict["ROC_AUC"] = acceptable_model_df[model + "_" + "ROC_AUC"].tolist()
                    tem_dict["Precision"] = acceptable_model_df[model + "_" + "Precision"].tolist()
                    tem_dict["Recall"] = acceptable_model_df[model + "_" + "Recall"].tolist()
                    tem_dict['Specificity'] = acceptable_model_df[model + "_" + 'Specificity'].tolist()
                    tem_dict["F1"] = acceptable_model_df[model + "_" + "F1"].tolist()
                    tem_dict['TN'] = acceptable_model_df[model + "_" + 'TN'].tolist()
                    tem_dict['FP'] = acceptable_model_df[model + "_" + 'FP'].tolist()
                    tem_dict['FN'] = acceptable_model_df[model + "_" + 'FN'].tolist()
                    tem_dict['TP'] = acceptable_model_df[model + "_" + 'TP'].tolist()



                    tem_dict["model"] = []
                    tem_dict["partition"] =[]
                    for _ in range(len(acceptable_model_df[model + "_" + "F1"].tolist())):
                         tem_dict["model"].append(model)
                         tem_dict["partition"].append("CV_partition_" + str(part_num))

                    for i in range(len(tem_dict["Accuracy"])):
                        # if tem_dict[Target_Metric][i] != "/":
                        #     if float(tem_dict[Target_Metric][i]) >= accepted_threshold:
                        ID = tem_dict["ID"][i]
                        partition_info = tem_dict["partition"][i]
                        model_info = tem_dict["model"][i]
                        feature_list_path = ML_path + "/" + gene_type + "/" + partition_info + "/" + str(ID) + "/" + model_info + "_imp_features.txt"
                        feature_list_df = pd.read_csv(feature_list_path, sep="\t")
                        if tem_dict["model"][i] == "stacking" or tem_dict["model"][i] == "randomForest":
                            feature_list = \
                                feature_list_df[(feature_list_df["feature_type"] == "M") & (
                                            feature_list_df["feature_score"] != 0)][
                                    "feature_name"].tolist()
                        else:
                            feature_list = feature_list_df[feature_list_df["feature_type"] == "M"][
                                "feature_name"].tolist()
                        feature_list.sort()
                        tem_model = tem_dict["model"][i]
                        tem_feature_num = tem_dict["feature_num"][i]
                        tem_Accuracy = tem_dict["Accuracy"][i]
                        tem_Balanced_Accuracy = tem_dict["Balanced_Accuracy"][i]
                        tem_ROC_AUC = tem_dict["ROC_AUC"][i]
                        tem_Precision = tem_dict["Precision"][i]
                        tem_Recall = tem_dict["Recall"][i]
                        tem_Specificity = tem_dict['Specificity'][i]
                        tem_F1 = tem_dict["F1"][i]
                        tem_TN = tem_dict["TN"][i]
                        tem_FP = tem_dict["FP"][i]
                        tem_FN = tem_dict["FN"][i]
                        tem_TP = tem_dict["TP"][i]
                        tem_partition = tem_dict["partition"][i]

                        tem_combination = (tem_model, tem_feature_num, feature_list, tem_partition)

                        tem_combination_pool = list(
                            zip(record_dict["model"], record_dict["feature_num"], record_dict["features"], record_dict["partition"]))

                        if tem_combination not in tem_combination_pool:
                            if len(feature_list) != 0:
                                record_dict["actual_feature_num"].append(len(feature_list))
                                record_dict["name"].append(str(ID) + "_"  + tem_model + "_" + str(tem_feature_num) + "_" + tem_partition)
                                record_dict["ID"].append(ID)
                                record_dict["model"].append(tem_model)
                                record_dict["feature_num"].append(tem_feature_num)
                                record_dict["Accuracy"].append(tem_Accuracy)
                                record_dict["Balanced_Accuracy"].append(tem_Balanced_Accuracy)
                                record_dict["ROC_AUC"].append(tem_ROC_AUC)
                                record_dict["Precision"].append(tem_Precision)
                                record_dict["Recall"].append(tem_Recall)
                                record_dict["Specificity"].append(tem_Specificity)
                                record_dict["F1"].append(tem_F1)
                                record_dict["TN"].append(tem_TN)
                                record_dict["FP"].append(tem_FP)
                                record_dict["FN"].append(tem_FN)
                                record_dict["TP"].append(tem_TP)
                                record_dict["features"].append(feature_list)
                                record_dict["partition"].append(tem_partition)



                        else:
                            tem_index = tem_combination_pool.index(tem_combination)
                            if record_dict[Target_Metric][tem_index] < tem_dict[Target_Metric][i]:
                                # record_dict["actual_feature_num"].append(len(feature_list))
                                record_dict["name"][tem_index] = str(ID) + "_"  + tem_model + "_" + str(tem_feature_num) + "_" + tem_partition
                                record_dict["ID"][tem_index] = ID
                                # record_dict["model"].append(tem_model)
                                # record_dict["feature_num"].append(tem_feature_num)
                                record_dict["Accuracy"][tem_index] = tem_Accuracy
                                record_dict["Balanced_Accuracy"][tem_index] = tem_Balanced_Accuracy
                                record_dict["ROC_AUC"][tem_index] = tem_ROC_AUC
                                record_dict["Precision"][tem_index] = tem_Precision
                                record_dict["Recall"][tem_index] = tem_Recall
                                record_dict["Specificity"][tem_index] = tem_Specificity
                                record_dict["F1"][tem_index] = tem_F1
                                record_dict["TN"][tem_index] = tem_TN
                                record_dict["FP"][tem_index] = tem_FP
                                record_dict["FN"][tem_index] = tem_FN
                                record_dict["TP"][tem_index] = tem_TP

                                # record_dict["features"].append(feature_list)



        # df_output_1 = pd.DataFrame({"name": record_dict["name"], "partition": record_dict["partition"],
        #                             "model": record_dict["model"], "effective_feature_num": record_dict["actual_feature_num"],
        #                             "Accuracy": record_dict["Accuracy"], "Balanced_Accuracy": record_dict["Balanced_Accuracy"],
        #                             "ROC_AUC":record_dict["ROC_AUC"], "Precision": record_dict["Precision"], "Recall": record_dict["Recall"],
        #                             "F1": record_dict["F1"], "genes":record_dict["features"]})

        df_output_1 = pd.DataFrame({"partition": record_dict["partition"],
                                    "model": record_dict["model"], "effective_feature_num": record_dict["actual_feature_num"],
                                    "Accuracy": record_dict["Accuracy"], "Balanced_Accuracy": record_dict["Balanced_Accuracy"],
                                    "ROC_AUC":record_dict["ROC_AUC"], "Precision": record_dict["Precision"], "Recall": record_dict["Recall"],
                                    "Specificity": record_dict["Specificity"], "F1": record_dict["F1"],
                                    # "TN": record_dict["TN"],"FP": record_dict["FP"],
                                    # "FN": record_dict["FN"],"TP": record_dict["TP"],
                                    "genes":record_dict["features"]})


        # df_output_1 = df_output_1.sort_values([Target_Metric, 'effective_feature_num'], ascending=[False, True])

        df_output_1 = df_output_1.sort_values(["partition",'effective_feature_num', Target_Metric], ascending=[True, True, False])



        tem_str = report_folder_path + "/" + gene_type + "_selected_models.csv"

        df_output_1.to_csv(tem_str, index = False)


        record_dict_2 = {}

        record_dict_2["model"] = []
        record_dict_2["effective_feature_num"] = []
        record_dict_2["feature_num"] = []
        record_dict_2["Accuracy"] = []
        record_dict_2["Balanced_Accuracy"] = []
        record_dict_2["ROC_AUC"] = []
        record_dict_2["Precision"] = []
        record_dict_2["Recall"] = []
        record_dict_2["Specificity"] = []
        record_dict_2["F1"] = []
        record_dict_2["TN"] = []
        record_dict_2["FP"] = []
        record_dict_2["FN"] = []
        record_dict_2["TP"] = []
        record_dict_2["name_list"] = []
        record_dict_2["features"] = []

        for i in range(len(record_dict["name"])):
            tem_feature_list = record_dict["features"][i]
            tem_name = record_dict["name"][i]
            tem_model = record_dict["model"][i]
            tem_feature_num = record_dict["feature_num"][i]
            tem_effective_feature_num = record_dict["actual_feature_num"][i]
            tem_Accuracy = record_dict["Accuracy"][i]
            tem_Balanced_Accuracy = record_dict["Balanced_Accuracy"][i]
            tem_ROC_AUC = record_dict["ROC_AUC"][i]
            tem_Precision = record_dict["Precision"][i]
            tem_Recall = record_dict["Recall"][i]
            tem_Specificity = record_dict["Specificity"][i]
            tem_F1 = record_dict["F1"][i]
            tem_TN = record_dict["TN"][i]
            tem_FP = record_dict["FP"][i]
            tem_FN = record_dict["FN"][i]
            tem_TP = record_dict["TP"][i]


            tem_combination = (tem_model, tem_feature_num, tem_feature_list)

            tem_combination_pool = list(zip(record_dict_2["model"], record_dict_2["feature_num"], record_dict_2["features"]))

            if tem_combination not in tem_combination_pool:
                if len(tem_feature_list) != 0:
                    record_dict_2["model"].append(tem_model)
                    record_dict_2["effective_feature_num"].append(tem_effective_feature_num)
                    record_dict_2["feature_num"].append(tem_feature_num)
                    record_dict_2["Accuracy"].append(tem_Accuracy)
                    record_dict_2["Balanced_Accuracy"].append(tem_Balanced_Accuracy)
                    record_dict_2["ROC_AUC"].append(tem_ROC_AUC)
                    record_dict_2["Precision"].append(tem_Precision)
                    record_dict_2["Recall"].append(tem_Recall)
                    record_dict_2["Specificity"].append(tem_Specificity)
                    record_dict_2["F1"].append(tem_F1)
                    record_dict_2["TN"].append(tem_TN)
                    record_dict_2["FP"].append(tem_FP)
                    record_dict_2["FN"].append(tem_FN)
                    record_dict_2["TP"].append(tem_TP)
                    record_dict_2["name_list"].append([tem_name])
                    record_dict_2["features"].append(tem_feature_list)





            else:
                tem_index = tem_combination_pool.index(tem_combination)
                # record_dict_2["model"].append(tem_model)
                # record_dict_2["effective_feature_num"].append(tem_effective_feature_num)
                # record_dict_2["feature_num"].append(tem_feature_num)
                record_dict_2["Accuracy"][tem_index] = record_dict_2["Accuracy"][tem_index] + tem_Accuracy
                record_dict_2["Balanced_Accuracy"][tem_index] = record_dict_2["Balanced_Accuracy"][tem_index] + tem_Balanced_Accuracy
                record_dict_2["ROC_AUC"][tem_index] = record_dict_2["ROC_AUC"][tem_index] + tem_ROC_AUC
                record_dict_2["Precision"][tem_index] = record_dict_2["Precision"][tem_index] + tem_Precision
                record_dict_2["Recall"][tem_index] = record_dict_2["Recall"][tem_index] + tem_Recall
                record_dict_2["Specificity"][tem_index] = record_dict_2["Specificity"][tem_index] + tem_Specificity
                record_dict_2["F1"][tem_index] = record_dict_2["F1"][tem_index] + tem_F1
                record_dict_2["TN"][tem_index] = record_dict_2["TN"][tem_index] + tem_TN
                record_dict_2["FP"][tem_index] = record_dict_2["FP"][tem_index] + tem_FP
                record_dict_2["FN"][tem_index] = record_dict_2["FN"][tem_index] + tem_FN
                record_dict_2["TP"][tem_index] = record_dict_2["TP"][tem_index] + tem_TP



                record_dict_2["name_list"][tem_index].append(tem_name)
                # record_dict_2["features"].append(tem_feature_list)



        record_dict_2["Accuracy_mean"] = []
        record_dict_2["Balanced_Accuracy_mean"] = []
        record_dict_2["ROC_AUC_mean"] = []
        record_dict_2["Precision_mean"] = []
        record_dict_2["Recall_mean"] = []
        record_dict_2["Specificity_mean"] = []
        record_dict_2["F1_mean"] = []
        record_dict_2["TN_mean"] = []
        record_dict_2["FP_mean"] = []
        record_dict_2["FN_mean"] = []
        record_dict_2["TP_mean"] = []
        record_dict_2["repeat_time"] = []



        for i in range(len(record_dict_2["model"])):
            name_list_length = len(record_dict_2["name_list"][i])
            record_dict_2["repeat_time"].append(name_list_length)
            record_dict_2["Accuracy_mean"].append(record_dict_2["Accuracy"][i]/name_list_length)
            record_dict_2["Balanced_Accuracy_mean"].append(record_dict_2["Balanced_Accuracy"][i]/name_list_length)
            record_dict_2["ROC_AUC_mean"].append(record_dict_2["ROC_AUC"][i]/name_list_length)
            record_dict_2["Precision_mean"].append(record_dict_2["Precision"][i]/name_list_length)
            record_dict_2["Recall_mean"].append(record_dict_2["Recall"][i]/name_list_length)
            record_dict_2["Specificity_mean"].append(record_dict_2["Specificity"][i] / name_list_length)
            record_dict_2["F1_mean"].append(record_dict_2["F1"][i]/name_list_length)
            record_dict_2["TN_mean"].append(record_dict_2["TN"][i] / name_list_length)
            record_dict_2["FP_mean"].append(record_dict_2["FP"][i] / name_list_length)
            record_dict_2["FN_mean"].append(record_dict_2["FN"][i] / name_list_length)
            record_dict_2["TP_mean"].append(record_dict_2["TP"][i] / name_list_length)


        df_output_2 = pd.DataFrame({"model":record_dict_2["model"], "effective_feature_num": record_dict_2["effective_feature_num"],
                                    "repeat_time_in_different_partitions": record_dict_2["repeat_time"],
                                    "Accuracy_mean":record_dict_2["Accuracy_mean"], "Balanced_Accuracy_mean":record_dict_2["Balanced_Accuracy_mean"],
                                    "ROC_AUC_mean":record_dict_2["ROC_AUC_mean"], "Precision_mean": record_dict_2["Precision_mean"],
                                    "Recall_mean":record_dict_2["Recall_mean"], "Specificity_mean": record_dict_2["Specificity_mean"],
                                    "F1_mean":record_dict_2["F1_mean"],
                                    # "TN_mean": record_dict_2["TN_mean"], "FP_mean":record_dict_2["FP_mean"],
                                    # "FN_mean": record_dict_2["FN_mean"], "TP_mean":record_dict_2["TP_mean"],
                                    "genes":record_dict_2["features"], "model_names":record_dict_2["name_list"]})



        df_output_2 = df_output_2.sort_values(["repeat_time_in_different_partitions", Target_Metric+"_mean", "effective_feature_num", "model"], ascending=[False, False, True, True])

        tem_str = report_folder_path + "/" + gene_type + "_selected_models_stat_between_partitions.csv"

        df_output_2.to_csv(tem_str, index = False)


        tem_calculate_list = []
        for tem_list in record_dict["features"]:
            tem_calculate_list = tem_calculate_list + tem_list
        stat_calculate_dict = Counter(tem_calculate_list)
        sorted_stat_calculate_dict = stat_calculate_dict.most_common()
        # print(sorted_stat_calculate_dict)

        if summary_type == "general":
            shared_gene_name_list = []
            frequency_list = []
            name_list = []
            model_record_list = []
            gene_num_list = []
            accuracy_list = []
            balance_accuracy_list = []
            roc_auc_list = []
            precision_list = []
            recall_list = []
            specificity_list = []
            f1_list = []
            tn_list = []
            fp_list = []
            fn_list = []
            tp_list = []



            for tup in sorted_stat_calculate_dict:
                if tup[1] > 1:
                    for i in range(len(record_dict["features"])):
                        tem_features = record_dict["features"][i]
                        if tup[0] in tem_features:
                            shared_gene_name_list.append(tup[0])
                            frequency_list.append(tup[1])
                            name_list.append(record_dict["name"][i])
                            model_record_list.append(record_dict["model"][i])
                            gene_num_list.append(record_dict["actual_feature_num"][i])
                            accuracy_list.append(record_dict["Accuracy"][i])
                            balance_accuracy_list.append(record_dict["Balanced_Accuracy"][i])
                            roc_auc_list.append(record_dict["ROC_AUC"][i])
                            precision_list.append(record_dict["Precision"][i])
                            recall_list.append(record_dict["Recall"][i])
                            specificity_list.append(record_dict["Specificity"][i])
                            f1_list.append(record_dict["F1"][i])
                            tn_list.append(record_dict["TN"][i])
                            fp_list.append(record_dict["FP"][i])
                            fn_list.append(record_dict["FN"][i])
                            tp_list.append(record_dict["TP"][i])

            df_output_3 = pd.DataFrame({"shared_gene": shared_gene_name_list, "frequency": frequency_list,
                                        "existed_in": name_list, "model": model_record_list, "effective_feature_num": gene_num_list,
                                        "Accuracy": accuracy_list, "Balanced_Accuracy": balance_accuracy_list,
                                        "ROC_AUC":roc_auc_list, "Precision": precision_list, "Recall": recall_list,
                                        "Specificity":specificity_list,"F1": f1_list,
                                        # "TN": tn_list, "FP": fp_list,
                                        # "FN": fn_list, "TP": tp_list,

                                        })

            # df_output_2 = df_output_2.applymap(tuple)
            # print(df_output_2)
            df_output_3_revised = df_output_3.sort_values(["frequency", "shared_gene", Target_Metric, 'effective_feature_num'], ascending=[False,True, False, True])

            tem_str = report_folder_path + "/" + gene_type + "_high_frequency_gene.csv"

            df_output_3_revised.to_csv(tem_str, index = False)

        # elif summary_type == "fresh":
        #     shared_gene_name_list = []
        #     frequency_list = []
        #     model_record_list = []
        #     gene_num_list = []
        #     accuracy_list = []
        #     balance_accuracy_list = []
        #     roc_auc_list = []
        #     precision_list = []
        #     recall_list = []
        #     f1_list = []
        #     model_type_num_list = []
        #
        #
        #     for tup in sorted_stat_calculate_dict:
        #         if tup[1] > 1:
        #             shared_gene_name_list.append(tup[0])
        #             frequency_list.append(tup[1])
        #             tem_model_list = []
        #
        #             gene_num_summary = 0
        #             accuracy_summary = 0
        #             balance_accuracy_summary = 0
        #             roc_auc_summary = 0
        #             precision_summary = 0
        #             recall_summary = 0
        #             f1_summary = 0
        #
        #             for i in range(len(record_dict["features"])):
        #                 tem_features = record_dict["features"][i]
        #                 if tup[0] in tem_features:
        #                     if record_dict["model"][i] not in tem_model_list:
        #                         tem_model_list.append(record_dict["model"][i])
        #                     gene_num_summary = gene_num_summary + record_dict["actual_feature_num"][i]
        #                     accuracy_summary = accuracy_summary + record_dict["Accuracy"][i]
        #                     balance_accuracy_summary = balance_accuracy_summary + record_dict["Balanced_Accuracy"][i]
        #                     roc_auc_summary = roc_auc_summary + record_dict["ROC_AUC"][i]
        #                     precision_summary = precision_summary + record_dict["Precision"][i]
        #                     recall_summary = recall_summary + record_dict["Recall"][i]
        #                     f1_summary = f1_summary + record_dict["F1"][i]
        #
        #             model_record_list.append(tem_model_list)
        #             model_type_num_list.append(len(tem_model_list))
        #             gene_num_list.append(gene_num_summary/int(tup[1]))
        #             accuracy_list.append(accuracy_summary/int(tup[1]))
        #             balance_accuracy_list.append(balance_accuracy_summary/int(tup[1]))
        #             roc_auc_list.append(roc_auc_summary/int(tup[1]))
        #             precision_list.append(precision_summary/int(tup[1]))
        #             recall_list.append(recall_summary/int(tup[1]))
        #             f1_list.append(f1_summary/int(tup[1]))
        #
        #
        #     df_output_3 = pd.DataFrame({"shared_gene": shared_gene_name_list, "frequency": frequency_list,
        #                                 "model": model_record_list, "effective_feature_num_mean": gene_num_list,
        #                                 "Accuracy_mean": accuracy_list, "Balanced_Accuracy_mean": balance_accuracy_list,
        #                                 "ROC_AUC_mean":roc_auc_list, "Precision_mean": precision_list, "Recall_mean": recall_list,
        #                                 "F1_mean": f1_list, "model_type_num":model_type_num_list})
        #
        #
        #     df_output_3 = df_output_3.sort_values(["frequency", "shared_gene", Target_Metric+"_mean", "effective_feature_num_mean", "model_type_num"], ascending=[False,True, False, True, False])
        #
        #     df_output_3_revised = df_output_3.drop(columns=["model_type_num"])

        elif summary_type == "fresh":
            shared_gene_name_list = []
            frequency_list = []
            model_record_list = []
            gene_num_list = []
            accuracy_list = []
            balance_accuracy_list = []
            roc_auc_list = []
            precision_list = []
            recall_list = []
            specificity_list = []
            f1_list = []
            tn_list = []
            fp_list = []
            fn_list = []
            tp_list = []
            model_type_num_list = []
            partition_list = []
            partition_num_list = []
            gene_num_average_list = []
            gene_num_20_list = []
            gene_num_10_list = []
            gene_num_20_mean_list = []
            gene_num_10_mean_list = []
            gene_num_20_partition_list = []
            gene_num_10_partition_list = []
            gene_num_20_model_list = []
            gene_num_10_model_list = []



            for tup in sorted_stat_calculate_dict:
                if tup[1] > 1:
                    shared_gene_name_list.append(tup[0])
                    frequency_list.append(tup[1])
                    tem_model_list = []
                    tem_partition_list = []
                    tem_gene_num_list = []
                    tem_partition_20_list = []
                    tem_partition_10_list = []
                    tem_model_20_list = []
                    tem_model_10_list = []

                    gene_num_summary = 0
                    accuracy_summary = 0
                    balance_accuracy_summary = 0
                    roc_auc_summary = 0
                    precision_summary = 0
                    recall_summary = 0
                    specificity_summary = 0
                    f1_summary = 0
                    tn_list_summary = 0
                    fp_list_summary = 0
                    fn_list_summary = 0
                    tp_list_summary = 0







                    tem_flag = 0

                    tem_dict_high_frequency = {}


                    for i in range(len(record_dict["features"])):
                        tem_features = record_dict["features"][i]
                        if tup[0] in tem_features:

                            if int(record_dict["partition"][i].split("_")[-1]) not in tem_partition_list:
                                tem_partition_list.append(int(record_dict["partition"][i].split("_")[-1]))

                            tem_gene_num_list.append(int(record_dict["actual_feature_num"][i]))

                            if int(record_dict["actual_feature_num"][i])<= 20:
                                if int(record_dict["partition"][i].split("_")[-1]) not in tem_partition_20_list:
                                    tem_partition_20_list.append(int(record_dict["partition"][i].split("_")[-1]))
                                if record_dict["model"][i] not in tem_model_20_list:
                                    tem_model_20_list.append(record_dict["model"][i])

                            if int(record_dict["actual_feature_num"][i])<= 10:
                                if int(record_dict["partition"][i].split("_")[-1]) not in tem_partition_10_list:
                                    tem_partition_10_list.append(int(record_dict["partition"][i].split("_")[-1]))
                                if record_dict["model"][i] not in tem_model_10_list:
                                    tem_model_10_list.append(record_dict["model"][i])


                            if record_dict["model"][i] not in tem_model_list:
                                tem_model_list.append(record_dict["model"][i])

                            if tem_flag == 0:
                                tem_dict_high_frequency["actual_feature_num"] = record_dict["actual_feature_num"][i]
                                tem_dict_high_frequency["Accuracy"] = record_dict["Accuracy"][i]
                                tem_dict_high_frequency["Balanced_Accuracy"] = record_dict["Balanced_Accuracy"][i]
                                tem_dict_high_frequency["ROC_AUC"] = record_dict["ROC_AUC"][i]
                                tem_dict_high_frequency["Precision"] = record_dict["Precision"][i]
                                tem_dict_high_frequency["Recall"] = record_dict["Recall"][i]
                                tem_dict_high_frequency["Specificity"] = record_dict["Specificity"][i]
                                tem_dict_high_frequency["F1"] = record_dict["F1"][i]
                                tem_dict_high_frequency["TN"] = record_dict["TN"][i]
                                tem_dict_high_frequency["FP"] = record_dict["FP"][i]
                                tem_dict_high_frequency["FN"] = record_dict["FN"][i]
                                tem_dict_high_frequency["TP"] = record_dict["TP"][i]
                                tem_flag = 1

                            else:
                                if record_dict["actual_feature_num"][i] <= tem_dict_high_frequency["actual_feature_num"]:
                                    if record_dict["actual_feature_num"][i] < tem_dict_high_frequency["actual_feature_num"]:
                                        tem_dict_high_frequency["actual_feature_num"] = record_dict["actual_feature_num"][i]
                                        tem_dict_high_frequency["Accuracy"] = record_dict["Accuracy"][i]
                                        tem_dict_high_frequency["Balanced_Accuracy"] = record_dict["Balanced_Accuracy"][i]
                                        tem_dict_high_frequency["ROC_AUC"] = record_dict["ROC_AUC"][i]
                                        tem_dict_high_frequency["Precision"] = record_dict["Precision"][i]
                                        tem_dict_high_frequency["Recall"] = record_dict["Recall"][i]
                                        tem_dict_high_frequency["Specificity"] = record_dict["Specificity"][i]
                                        tem_dict_high_frequency["F1"] = record_dict["F1"][i]
                                        tem_dict_high_frequency["TN"] = record_dict["TN"][i]
                                        tem_dict_high_frequency["FP"] = record_dict["FP"][i]
                                        tem_dict_high_frequency["FN"] = record_dict["FN"][i]
                                        tem_dict_high_frequency["TP"] = record_dict["TP"][i]
                                    elif record_dict["actual_feature_num"][i] == tem_dict_high_frequency["actual_feature_num"]:
                                        if record_dict[Target_Metric][i] > tem_dict_high_frequency[Target_Metric]:
                                            tem_dict_high_frequency["actual_feature_num"] = record_dict["actual_feature_num"][i]
                                            tem_dict_high_frequency["Accuracy"] = record_dict["Accuracy"][i]
                                            tem_dict_high_frequency["Balanced_Accuracy"] = record_dict["Balanced_Accuracy"][i]
                                            tem_dict_high_frequency["ROC_AUC"] = record_dict["ROC_AUC"][i]
                                            tem_dict_high_frequency["Precision"] = record_dict["Precision"][i]
                                            tem_dict_high_frequency["Recall"] = record_dict["Recall"][i]
                                            tem_dict_high_frequency["Specificity"] = record_dict["Specificity"][i]
                                            tem_dict_high_frequency["F1"] = record_dict["F1"][i]
                                            tem_dict_high_frequency["TN"] = record_dict["TN"][i]
                                            tem_dict_high_frequency["FP"] = record_dict["FP"][i]
                                            tem_dict_high_frequency["FN"] = record_dict["FN"][i]
                                            tem_dict_high_frequency["TP"] = record_dict["TP"][i]



                    model_record_list.append(tem_model_list)
                    model_type_num_list.append(len(tem_model_list))
                    gene_num_list.append(tem_dict_high_frequency["actual_feature_num"])
                    accuracy_list.append(tem_dict_high_frequency["Accuracy"])
                    balance_accuracy_list.append(tem_dict_high_frequency["Balanced_Accuracy"])
                    roc_auc_list.append(tem_dict_high_frequency["ROC_AUC"])
                    precision_list.append(tem_dict_high_frequency["Precision"])
                    recall_list.append(tem_dict_high_frequency["Recall"])
                    specificity_list.append(tem_dict_high_frequency["Specificity"])
                    f1_list.append(tem_dict_high_frequency["F1"])
                    tn_list.append(tem_dict_high_frequency["TN"])
                    fp_list.append(tem_dict_high_frequency["FP"])
                    fn_list.append(tem_dict_high_frequency["FN"])
                    tp_list.append(tem_dict_high_frequency["TP"])

                    partition_list.append(tem_partition_list)
                    partition_num_list.append(len(tem_partition_list))

                    gene_num_average_list.append(sum(tem_gene_num_list)/len(tem_gene_num_list))

                    tem_less_than_20_list = [t for t in tem_gene_num_list if t <= 20]
                    tem_less_than_10_list = [t for t in tem_gene_num_list if t <= 10]

                    gene_num_20_list.append(len(tem_less_than_20_list ))
                    gene_num_10_list.append(len(tem_less_than_10_list ))

                    if len(tem_less_than_20_list)!= 0:
                        gene_num_20_mean_list.append(sum(tem_less_than_20_list)/len(tem_less_than_20_list))
                    else:
                        gene_num_20_mean_list.append(0)
                    if len(tem_less_than_10_list) != 0:
                        gene_num_10_mean_list.append(sum(tem_less_than_10_list)/len(tem_less_than_10_list))
                    else:
                        gene_num_10_mean_list.append(0)

                    gene_num_20_partition_list.append(len(tem_partition_20_list))
                    gene_num_10_partition_list.append(len(tem_partition_10_list))

                    gene_num_20_model_list.append(len(tem_model_20_list))
                    gene_num_10_model_list.append(len(tem_model_10_list))


            df_output_3 = pd.DataFrame({"shared_gene": shared_gene_name_list, "frequency": frequency_list,
                                        "model": model_record_list, "effective_feature_num_min": gene_num_list,
                                        "Accuracy": accuracy_list, "Balanced_Accuracy": balance_accuracy_list,
                                        "ROC_AUC": roc_auc_list, "Precision": precision_list,
                                        "Recall": recall_list, "Specificity": specificity_list,
                                        "F1": f1_list,
                                        "TN": tn_list, "FP": fp_list,
                                        "FN": fn_list, "TP": tp_list,
                                        "model_type_num": model_type_num_list})




            tem_Accuracy_record_list = []
            tem_Balanced_Accuracy_record_list = []
            tem_ROC_AUC_record_list = []
            tem_Precision_record_list = []
            tem_Recall_record_list = []
            tem_Specificity_record_list = []
            tem_F1_record_list = []
            tem_TN_record_list = []
            tem_FP_record_list = []
            tem_FN_record_list = []
            tem_TP_record_list = []
            tem_repeat_record_list = []
            tem_model_name_record_list = []


            for j in range(len(shared_gene_name_list)):
                high_freq_gene = shared_gene_name_list[j]
                effect_gene_num = gene_num_list[j]
                tem_similar_number_list = []
                tem_Accuracy_record = 0
                tem_Balanced_Accuracy_record = 0
                tem_ROC_AUC_record = 0
                tem_Precision_record = 0
                tem_Recall_record = 0
                tem_Specificity_record = 0
                tem_F1_record = 0
                tem_TN_record = 0
                tem_FP_record = 0
                tem_FN_record = 0
                tem_TP_record = 0

                for k in range(len(record_dict["features"])):
                    tem_features = record_dict["features"][k]
                    if high_freq_gene in tem_features:
                        if record_dict["actual_feature_num"][k] == effect_gene_num:
                            tem_similar_number_list.append(record_dict["name"][k])
                            tem_Accuracy_record = record_dict["Accuracy"][k] + tem_Accuracy_record
                            tem_Balanced_Accuracy_record = record_dict["Balanced_Accuracy"][k] + tem_Balanced_Accuracy_record
                            tem_ROC_AUC_record = record_dict["ROC_AUC"][k] + tem_ROC_AUC_record
                            tem_Precision_record = record_dict["Precision"][k] + tem_Precision_record
                            tem_Recall_record = record_dict["Recall"][k] + tem_Recall_record
                            tem_Specificity_record = record_dict["Specificity"][k] + tem_Specificity_record
                            tem_F1_record = record_dict["F1"][k] + tem_F1_record
                            tem_TN_record = record_dict["TN"][k] + tem_TN_record
                            tem_FP_record = record_dict["FP"][k] + tem_FP_record
                            tem_FN_record = record_dict["FN"][k] + tem_FN_record
                            tem_TP_record = record_dict["TP"][k] + tem_TP_record

                repeat_time = len(tem_similar_number_list)
                tem_repeat_record_list.append(repeat_time)
                tem_model_name_record_list.append(tem_similar_number_list)

                tem_Accuracy_record_list.append(tem_Accuracy_record/repeat_time)
                tem_Balanced_Accuracy_record_list.append(tem_Balanced_Accuracy_record/repeat_time)
                tem_ROC_AUC_record_list.append(tem_ROC_AUC_record/repeat_time)
                tem_Precision_record_list.append(tem_Precision_record/repeat_time)
                tem_Recall_record_list.append(tem_Recall_record/repeat_time)
                tem_Specificity_record_list.append(tem_Specificity_record / repeat_time)
                tem_F1_record_list.append(tem_F1_record/repeat_time)
                tem_TN_record_list.append(tem_TN_record / repeat_time)
                tem_FP_record_list.append(tem_FP_record / repeat_time)
                tem_FN_record_list.append(tem_FN_record / repeat_time)
                tem_TP_record_list.append(tem_TP_record / repeat_time)

            df_output_3_organized = pd.DataFrame({"shared_gene": shared_gene_name_list, "frequency": frequency_list,
                                        "model": model_record_list, "model_type_num": model_type_num_list,
                                        "covered_partitions":partition_list, "covered_partitions_num":partition_num_list,
                                        "avg_effective_feature_num":gene_num_average_list,
                                        "num_effective_feature_<=20":gene_num_20_list,
                                        "partition_num_effective_feature_<=20": gene_num_20_partition_list,
                                        "model_type_effective_feature_<=20": gene_num_20_model_list,
                                        "avg_effective_feature_<=20":gene_num_20_mean_list,
                                        "num_effective_feature_<=10": gene_num_10_list,
                                        "partition_num_effective_feature_<=10": gene_num_10_partition_list,
                                        "model_type_effective_feature_<=10": gene_num_10_model_list,
                                        "avg_effective_feature_<=10": gene_num_10_mean_list,
                                        "effective_feature_num_min": gene_num_list, "min_gene_model_repeat": tem_repeat_record_list,
                                        "Accuracy_mean": tem_Accuracy_record_list, "Balanced_Accuracy_mean": tem_Balanced_Accuracy_record_list,
                                        "ROC_AUC_mean": tem_ROC_AUC_record_list, "Precision_mean": tem_Precision_record_list,
                                        "Recall_mean": tem_Recall_record_list,
                                        "Specificity_mean": tem_Specificity_record_list,
                                        "F1_mean": tem_F1_record_list,
                                        # "TN_mean": tem_TN_record_list,"FP_mean": tem_FP_record_list,
                                        # "FN_mean": tem_FN_record_list,"TP_mean": tem_TP_record_list,
                                        "names":tem_model_name_record_list})


            # df_output_3_organized = df_output_3_organized.sort_values(["frequency",  "effective_feature_num_min", Target_Metric, "model_type_num", "shared_gene"],ascending=[False, True, False, False, True])

            df_output_3_organized = df_output_3_organized.sort_values(
                ["effective_feature_num_min", "min_gene_model_repeat", "frequency", Target_Metric+"_mean", "model_type_num", "shared_gene"],
                ascending=[True, False, False, False, False, True])






            # df_output_3_revised = df_output_3.drop(columns=["model_type_num"])







            tem_str = report_folder_path + "/" + gene_type + "_high_frequency_gene.csv"

            df_output_3_organized.to_csv(tem_str, index = False)


            ###heatmap generation

        def add_information_heatmap(padj_path, FPKM_path, feature_imp_df, gene_type):

            feature_imp_df =feature_imp_df[["shared_gene", "frequency"]].set_index("shared_gene")

            padj_df = pd.read_csv(padj_path, sep="\t", index_col=0)
            padj_df = pd.merge(feature_imp_df, padj_df, left_index=True, right_index=True, how='left')

            FPKM_df = pd.read_csv(FPKM_path, sep="\t", index_col=0)

            FPKM_df = pd.merge(feature_imp_df, FPKM_df, left_index=True, right_index=True, how='left').drop(
                ["frequency"], axis=1)

            FPKM_df_transposition = FPKM_df.T

            Label_list = []

            for i in range(len(FPKM_df_transposition)):
                if str(FPKM_df_transposition.index.values[i]).endswith('_C'):
                    Label_list.append(0)
                elif str(FPKM_df_transposition.index.values[i]).endswith('_T'):
                    Label_list.append(1)
                else:
                    print("Label not match")

            FPKM_df_transposition['label'] = Label_list
            FPKM_df_transposition_0 = FPKM_df_transposition[FPKM_df_transposition['label'].isin([0])]
            FPKM_df_transposition_1 = FPKM_df_transposition[FPKM_df_transposition['label'].isin([1])]

            C_number = len(FPKM_df_transposition_0)
            T_number = len(FPKM_df_transposition_1)

            FPKM_df_transposition_0_transform = FPKM_df_transposition_0.drop(['label'], axis=1).T
            FPKM_df_transposition_1_transform = FPKM_df_transposition_1.drop(['label'], axis=1).T

            FPKM_df = pd.merge(FPKM_df_transposition_0_transform, FPKM_df_transposition_1_transform,
                               left_index=True,
                               right_index=True)

            return padj_df, FPKM_df, C_number, T_number

        # PC_padj = PC_padj_path
        # PC_FPKM = PC_FPKM_path
        # linc_padj = linc_padj_path
        # linc_FPKM = linc_FPKM_path
        # Circ_padj = Circ_padj_path
        # Circ_FPKM = Circ_FPKM_path

        # Heatmap_info_dict["protein_coding"] = {}
        # Heatmap_info_dict["protein_coding"]["padj"] = output_directory + "/" + DE_type + "_DE_output/summary_step1.txt.ReadCount.ProteinCoding.DE.output"
        # Heatmap_info_dict["protein_coding"]["FPKM"] = output_directory + "/summary_step1.txt.FPKM.ProteinCoding"
        # Heatmap_info_dict["lincRNA"]["padj"] = output_directory + "/" + DE_type + "_DE_output/summary_step1.txt.ReadCount.lincRNA.DE.output"
        # Heatmap_info_dict["lincRNA"]["FPKM"] = output_directory + "/summary_step1.txt.FPKM.lincRNA"




        # feature_imp_path = ML_path + "/" + ID + "/best_model_imp_features.txt"
        # feature_imp_df = pd.read_csv(feature_imp_path, sep="\t", index_col=0)
        # linc_flag = 0
        # PC_flag = 0
        # Circ_flag = 0
        # add_info_flag = 0

        padj_df, FPKM_df, C_number, T_number = add_information_heatmap(Heatmap_info_dict[gene_type]["padj"], Heatmap_info_dict[gene_type]["FPKM"],
                                                                                 df_output_3_organized, gene_type)

        Heatmap_file_folder = report_folder_path + "/Heatmap_files"

        if not os.path.isdir(Heatmap_file_folder):
            os.mkdir(Heatmap_file_folder)

        total_padj_df = padj_df
        total_FPKM_df = FPKM_df
        total_padj_df.index.name = None
        padj_save_path = Heatmap_file_folder + "/" + gene_type + "_biomarker_padj_df.txt"
        FPKM_save_path = Heatmap_file_folder + "/" + gene_type + "_biomarker_FPKM_df.txt"

        if os.path.isfile(padj_save_path):
            tem_str = "rm " + padj_save_path
            os.system(tem_str)
        if os.path.isfile(FPKM_save_path):
            tem_str = "rm " + FPKM_save_path
            os.system(tem_str)
        total_padj_df.to_csv(padj_save_path, sep="\t", na_rep='NA')
        total_FPKM_df.to_csv(FPKM_save_path, sep="\t", na_rep='NA')
        os.chdir(Heatmap_file_folder)
        R_script_path = src_path + "/Heatmap_High_Frequency_Gene_General.R"
        # heatmap_generation = subprocess.Popen(
        #     ['Rscript', R_script_path, str(C_number), gene_type + "_biomarker_padj_df.txt",
        #      gene_type + "_biomarker_FPKM_df.txt", gene_type + "_heatmap.pdf"
        #      #   , pvalue_type, pvalue_threshold,linc_min_FPKM_expression
        #      ])
        # heatmap_generation.communicate()
        os.chdir(report_folder_path)
        # result_heatmap_path = Heatmap_file_folder + "/" + gene_type + "_heatmap.pdf"


        # if os.path.isfile(result_heatmap_path):
        #     # tem_str = "ID:" + ID + " Balanced_Accuracy:" + str(
        #     #     report_df.iloc[0].at['Balanced_Accuracy']) + " Best Model:" + \
        #     #           str(report_df.iloc[0].at['BestModel']) + " feature number：" + str(
        #     #     report_df.iloc[0].at['feature_num'])
        #     # output_report_doc.add_paragraph(tem_str)
        #     images = convert_from_path(result_heatmap_path)
        #     tem_pic_path = Heatmap_file_folder + "/" + gene_type + "_heatmap.jpg"
        #     images[0].save(tem_pic_path, 'JPEG')
        #     # pic_no = pic_no + 1
        #     # output_report_doc.add_picture(tem_pic_path, width=Inches(6.42), height=Inches(3.21))
        #     # output_report_doc = gene_name_translation(gene_name_trans_file, ID, ML_path, output_report_doc)
        # else:
        #     print(gene_type + "'s heatmap cannot be generated")




############ new intersection check function
        if intersection_check == int(1):
            if not os.path.isdir(report_folder_path + "/tem"):
                os.mkdir(report_folder_path + "/tem")

            t = time.time()
            print("2 intersection begin")

            intersection_model_list = []
            intersection_model_num_list = []
            intersection_gene_num_list = []
            intersection_gene_list = []

            #############################
            # all_shared_genes_list = record_dict["features"].copy()
            # all_model_combination_list = record_dict["name"].copy()

            all_shared_genes_list = record_dict_2["features"].copy()
            all_model_combination_list = []
            for i in range(len(record_dict_2["name_list"])):
                all_model_combination_list.append("repeat"+str(record_dict_2["repeat_time"][i])+"_"+str(record_dict_2["name_list"][i][0]))




            for i in range(len(all_shared_genes_list)):
                for j in range(i):
                    tem_list = []
                    tem_list.append(all_shared_genes_list[i])
                    tem_list.append(all_shared_genes_list[j])
                    share_gene = set.intersection(*map(set, tem_list))
                    num = len(share_gene)
                    if num > 1:
                        tem_name_list = []
                        tem_name_list.append(all_model_combination_list[i])
                        tem_name_list.append(all_model_combination_list[j])
                        tem_name_list.sort()
                        tem_share_gene_list = list(share_gene)
                        tem_share_gene_list.sort()
                        intersection_gene_num_list.append(num)
                        intersection_model_list.append(tem_name_list)
                        intersection_model_num_list.append(2)
                        intersection_gene_list.append(tem_share_gene_list)

            tem_df_2 = pd.DataFrame({"shared_gene_num": intersection_gene_num_list,
                                     "model_num": intersection_model_num_list,
                                     "shared_genes": intersection_gene_list,
                                     "model_combination": intersection_model_list}).sort_values(['model_num', 'shared_gene_num'], ascending=[False, False])

            tem_str = report_folder_path + "/tem/" + str(2) + "_" + gene_type + "_gene_intersection_relationship.csv"
            tem_df_2.to_csv(tem_str, index=False)

            last_intersect_model_num = 2

            # intersect_model_num = 2

            print("2 intersection finished")
            print(len(intersection_gene_list))
            print("time: " + str(time.time() - t))

            t = time.time()
            print("2 intersection dictionary generation")
            tem_record_dict = {}

            for i in range(len(all_model_combination_list)):
                if all_model_combination_list[i] not in tem_record_dict:
                    tem_record_dict[all_model_combination_list[i]] = []
                for j in range(len(intersection_model_list)):
                    if all_model_combination_list[i] in intersection_model_list[j]:
                        for model_name in intersection_model_list[j]:
                            if model_name != all_model_combination_list[i]:
                                tem_record_dict[all_model_combination_list[i]].append(model_name)


            print("2 intersection dictionary generation finished")
            print("time: " + str(time.time() - t))


            t = time.time()
            print("3 intersection begin")


            shared_genes_list_2 = intersection_gene_list.copy()
            model_combination_list_2 = intersection_model_list.copy()



            intersection_model_list = []
            intersection_model_num_list = []
            intersection_gene_num_list = []
            intersection_gene_list = []

            for i in range(len(model_combination_list_2)):
                model_combination = model_combination_list_2[i]
                ori_shared_gene = shared_genes_list_2[i]
                model_A = model_combination[0]
                model_B = model_combination[1]
                tem_list = []
                tem_list.append(tem_record_dict[model_A])
                tem_list.append(tem_record_dict[model_B])
                candidate_model_list = list(set.intersection(*map(set, tem_list)))
                for j in range(len(candidate_model_list)):
                    candidate_model = candidate_model_list[j]
                    if candidate_model not in model_combination:
                        tem_name_list = model_combination.copy()
                        tem_name_list.append(candidate_model)
                        tem_name_list.sort()
                        if tem_name_list not in intersection_model_list:
                            tem_list = []
                            tem_list.append(ori_shared_gene)
                            index = all_model_combination_list.index(candidate_model)
                            tem_list.append(all_shared_genes_list[index])
                            share_gene = set.intersection(*map(set, tem_list))
                            num = len(share_gene)
                            if num > 1:
                                tem_share_gene_list = list(share_gene)
                                tem_share_gene_list.sort()
                                intersection_gene_num_list.append(num)
                                intersection_model_list.append(tem_name_list)
                                intersection_model_num_list.append(3)
                                intersection_gene_list.append(tem_share_gene_list)


            tem_df_3 = pd.DataFrame({"shared_gene_num": intersection_gene_num_list,
                                     "model_num": intersection_model_num_list,
                                     "shared_genes": intersection_gene_list,
                                     "model_combination": intersection_model_list}).sort_values(['model_num', 'shared_gene_num'], ascending=[False, False])

            tem_str = report_folder_path + "/tem/" + str(3) + "_" + gene_type + "_gene_intersection_relationship.csv"
            tem_df_3.to_csv(tem_str, index=False)

            if len(tem_df_3) != 0:

                last_intersect_model_num = 3

                # intersect_model_num = 2

                print("3 intersection finished")
                print(len(intersection_gene_list))
                print("time: " + str(time.time() - t))

                t = time.time()
                print("3 intersection dictionary generation")

                tem_record_dict = {}

                for i in range(len(model_combination_list_2)):
                    tem_name = ""
                    for name in model_combination_list_2[i]:
                        tem_name = tem_name + name + ","
                    if tem_name not in tem_record_dict:
                        tem_record_dict[tem_name] = []
                    for j in range(len(intersection_model_list)):
                        flag = 0
                        for k in range(len(model_combination_list_2[i])):
                            if model_combination_list_2[i][k] in intersection_model_list[j]:
                                flag = flag + 1
                        if flag == len(model_combination_list_2[i]):
                            for model_name in intersection_model_list[j]:
                                if model_name not in model_combination_list_2[i]:
                                    tem_record_dict[tem_name].append(model_name)
                    if len(tem_record_dict[tem_name]) == 0:
                        tem_record_dict.pop(tem_name, None)

                print("3 intersection dictionary generation finished")
                print("time: " + str(time.time() - t))


                def func(current_intersection_gene_list, current_intersection_model_list,tem_record_dict):
                    output_pool = {}
                    output_pool["tem_intersection_model_name_list"] = []
                    output_pool["tem_intersection_model_list"] = []
                    output_pool["tem_intersection_gene_list"] = []
                    output_pool["tem_intersection_model_add_list"] = []

                    model_number_current_pursue = len(ast.literal_eval(current_intersection_model_list[0])) + 1

                    for i in range(len(current_intersection_model_list)):
                        current_intersection_model = ast.literal_eval(current_intersection_model_list[i])
                        index_list = [*range(model_number_current_pursue - 1)]
                        combination_list = list(combinations(index_list, model_number_current_pursue - 2))
                        record_combination_model_info = []
                        flag = 0
                        for combination in combination_list:
                            tem_name_list = []
                            for tem_index in combination:
                                tem_name_list.append(current_intersection_model[tem_index])
                            tem_name_list.sort()
                            tem_name = ""
                            for name in tem_name_list:
                                tem_name = tem_name + name + ","
                            curent_combination_model_info = tem_record_dict[tem_name]
                            if len(record_combination_model_info) == 0:
                                record_combination_model_info = curent_combination_model_info.copy()
                            else:
                                tem_list = []
                                tem_list.append(record_combination_model_info)
                                tem_list.append(curent_combination_model_info)
                                share_model = list(set.intersection(*map(set, tem_list)))
                                if len(share_model) == 0:
                                    flag = 1
                                    break
                                else:
                                    record_combination_model_info = share_model.copy()

                        if flag == 0:
                            for j in range(len(record_combination_model_info)):
                                if record_combination_model_info[j] not in ast.literal_eval(current_intersection_model_list[i]):
                                    tem_name_list = ast.literal_eval(current_intersection_model_list[i]).copy()
                                    tem_name_list.append(record_combination_model_info[j])
                                    tem_name_list.sort()
                                    str_tem_name_list = ' '.join(str(e) for e in tem_name_list)
                                    if str_tem_name_list not in output_pool["tem_intersection_model_name_list"]:
                                        output_pool["tem_intersection_model_list"].append(tem_name_list)
                                        output_pool["tem_intersection_model_name_list"].append(str_tem_name_list)
                                        output_pool["tem_intersection_gene_list"].append(ast.literal_eval(current_intersection_gene_list[i]))
                                        output_pool["tem_intersection_model_add_list"].append(record_combination_model_info[j])

                    return output_pool

                def func_2(candidate_tem_intersection_model_list, candidate_tem_intersection_gene_list, candidate_tem_intersection_model_add_list, all_shared_genes_list, all_model_combination_list):
                    output_pool = {}
                    output_pool["tem_intersection_model_list"] = []
                    output_pool["tem_intersection_model_num_list"] = []
                    output_pool["tem_intersection_gene_num_list"] = []
                    output_pool["tem_intersection_gene_list"] = []

                    for i in range(len(candidate_tem_intersection_model_list)):
                        tem_list = []
                        tem_list.append(candidate_tem_intersection_gene_list[i])
                        tem_index = all_model_combination_list.index(candidate_tem_intersection_model_add_list[i])
                        tem_list.append(all_shared_genes_list[tem_index])
                        share_gene = set.intersection(*map(set, tem_list))
                        num = len(share_gene)
                        if num > 1:
                            tem_share_gene_list = list(share_gene)
                            tem_share_gene_list.sort()
                            output_pool["tem_intersection_gene_num_list"].append(num)
                            output_pool["tem_intersection_model_list"].append(candidate_tem_intersection_model_list[i])
                            output_pool["tem_intersection_model_num_list"].append(len(candidate_tem_intersection_model_list[i]))
                            output_pool["tem_intersection_gene_list"].append(tem_share_gene_list)


                    return output_pool

                def func_3(current_intersection_model_list, tem_intersection_model_list):
                    tem_record_dict = {}

                    for i in range(len(current_intersection_model_list)):
                        tem_name = ""
                        transform_current_intersection_model = ast.literal_eval(current_intersection_model_list[i])
                        for name in transform_current_intersection_model:
                            tem_name = tem_name + name + ","
                        if tem_name not in tem_record_dict:
                            tem_record_dict[tem_name] = []
                        for j in range(len(tem_intersection_model_list)):
                            flag = 0
                            for k in range(len(transform_current_intersection_model)):
                                if transform_current_intersection_model[k] in tem_intersection_model_list[j]:
                                    flag = flag + 1
                            if flag == len(transform_current_intersection_model):
                                for model_name in tem_intersection_model_list[j]:
                                    if model_name not in transform_current_intersection_model:
                                        tem_record_dict[tem_name].append(model_name)
                        if len(tem_record_dict[tem_name]) == 0:
                            tem_record_dict.pop(tem_name, None)

                    return tem_record_dict

                # all_shared_genes_list = record_dict["features"].copy()
                # all_model_combination_list = record_dict["name"].copy()


                for intersect_model_num in range(4, len(all_shared_genes_list)+1):
                    # intersect_model_num = intersect_model_num + 1
                    # if intersect_model_num == 3:

                    t = time.time()
                    print(str(intersect_model_num) + " intersection begin")


                    tem_df = pd.read_csv(report_folder_path + "/tem/" + str(intersect_model_num-1) + "_" + gene_type + "_gene_intersection_relationship.csv")

                    current_intersection_model_list = tem_df["model_combination"].tolist()
                    # current_intersection_model_num_list = intersection_model_num_list
                    # current_intersection_gene_num_list = intersection_gene_num_list
                    current_intersection_gene_list = tem_df["shared_genes"].tolist()

                    num_cores = int(mp.cpu_count())


                    if len(current_intersection_model_list) >= num_cores:
                        separation_length = int(len(current_intersection_model_list) / num_cores)
                        last_end_index = 0
                        tem_list_dict = {}
                        for i in range(num_cores-1):
                            tem_list_dict["model_list_" + str(i)] = current_intersection_model_list[last_end_index:(last_end_index + separation_length)]
                            tem_list_dict["gene_list_" + str(i)] = current_intersection_gene_list[last_end_index:(last_end_index + separation_length)]
                            last_end_index = last_end_index + separation_length

                        tem_list_dict["model_list_" + str(i+1)] = current_intersection_model_list[last_end_index:]
                        tem_list_dict["gene_list_" + str(i+1)] = current_intersection_gene_list[last_end_index:]

                        pool = mp.Pool(num_cores)
                        process_pool = {}
                        output_pool_1 = {}
                        output_pool_1["tem_intersection_model_name_list"] = []
                        output_pool_1["tem_intersection_model_list"] = []
                        output_pool_1["tem_intersection_gene_list"] = []
                        output_pool_1["tem_intersection_model_add_list"] = []
                        for i in range(num_cores):
                            process_pool["process_"+str(i)] = pool.apply_async(func, args=(tem_list_dict["gene_list_" + str(i)], tem_list_dict["model_list_" + str(i)], tem_record_dict))

                        pool.close()
                        pool.join()

                        for key in process_pool.keys():
                            tem_collect_dict = process_pool[key].get()
                            output_pool_1["tem_intersection_model_name_list"].extend(tem_collect_dict["tem_intersection_model_name_list"])
                            output_pool_1["tem_intersection_model_list"].extend(tem_collect_dict["tem_intersection_model_list"])
                            output_pool_1["tem_intersection_gene_list"].extend(tem_collect_dict["tem_intersection_gene_list"])
                            output_pool_1["tem_intersection_model_add_list"].extend(tem_collect_dict["tem_intersection_model_add_list"])

                        tem_df_candidate_combination = pd.DataFrame(output_pool_1)
                        tem_df_candidate_combination.drop_duplicates(subset=["tem_intersection_model_name_list"], keep='first', ignore_index=True, inplace = True)

                        candidate_tem_intersection_model_list = tem_df_candidate_combination["tem_intersection_model_list"].tolist()
                        candidate_tem_intersection_gene_list = tem_df_candidate_combination["tem_intersection_gene_list"].tolist()
                        candidate_tem_intersection_model_add_list = tem_df_candidate_combination["tem_intersection_model_add_list"].tolist()

                        if len(candidate_tem_intersection_model_list) >= num_cores:
                            separation_length = int(len(candidate_tem_intersection_model_list) / num_cores)
                            last_end_index = 0
                            tem_list_dict = {}
                            for i in range(num_cores - 1):
                                tem_list_dict["model_list_" + str(i)] = candidate_tem_intersection_model_list[last_end_index:(last_end_index + separation_length)]
                                tem_list_dict["gene_list_" + str(i)] = candidate_tem_intersection_gene_list[last_end_index:(last_end_index + separation_length)]
                                tem_list_dict["model_add_list_" + str(i)] = candidate_tem_intersection_model_add_list[last_end_index:(last_end_index + separation_length)]
                                last_end_index = last_end_index + separation_length

                            tem_list_dict["model_list_" + str(i + 1)] = candidate_tem_intersection_model_list[last_end_index:]
                            tem_list_dict["gene_list_" + str(i + 1)] = candidate_tem_intersection_gene_list[last_end_index:]
                            tem_list_dict["model_add_list_" + str(i + 1)] = candidate_tem_intersection_model_add_list[last_end_index:]

                            pool = mp.Pool(num_cores)
                            process_pool = {}
                            output_pool = {}
                            for i in range(num_cores):
                                process_pool["process_" + str(i)] = pool.apply_async(func_2, args=(tem_list_dict["model_list_" + str(i)],
                                                                                                   tem_list_dict["gene_list_" + str(i)],
                                                                                                   tem_list_dict["model_add_list_" + str(i)],
                                                                                                   all_shared_genes_list, all_model_combination_list))

                            pool.close()
                            pool.join()

                            tem_intersection_model_list = []
                            tem_intersection_model_num_list = []
                            tem_intersection_gene_num_list = []
                            tem_intersection_gene_list = []

                            for key in process_pool.keys():
                                tem_collect_dict = process_pool[key].get()
                                tem_intersection_model_list.extend(tem_collect_dict["tem_intersection_model_list"])
                                tem_intersection_model_num_list.extend(tem_collect_dict["tem_intersection_model_num_list"])
                                tem_intersection_gene_num_list.extend(tem_collect_dict["tem_intersection_gene_num_list"])
                                tem_intersection_gene_list.extend(tem_collect_dict["tem_intersection_gene_list"])

                        else:
                            tem_intersection_model_list = []
                            tem_intersection_model_num_list = []
                            tem_intersection_gene_num_list = []
                            tem_intersection_gene_list = []

                            for i in range(len(candidate_tem_intersection_model_list)):
                                tem_list = []
                                tem_list.append(candidate_tem_intersection_gene_list[i])
                                tem_index = all_model_combination_list.index(candidate_tem_intersection_model_add_list[i])
                                tem_list.append(all_shared_genes_list[tem_index])
                                share_gene = set.intersection(*map(set, tem_list))
                                num = len(share_gene)
                                if num > 1:
                                    tem_share_gene_list = list(share_gene)
                                    tem_share_gene_list.sort()
                                    tem_intersection_gene_num_list.append(num)
                                    tem_intersection_model_list.append(candidate_tem_intersection_model_list[i])
                                    tem_intersection_model_num_list.append(len(candidate_tem_intersection_model_list[i]))
                                    tem_intersection_gene_list.append(tem_share_gene_list)



                        # for j in range(num_cores):
                        #     for k in range(len(output_pool[str(j)]["tem_intersection_model_list"])):
                        #         if output_pool[str(j)]["tem_intersection_model_list"][k] not in tem_intersection_model_list:
                        #             tem_intersection_model_list.append(output_pool[str(j)]["tem_intersection_model_list"][k])
                        #             tem_intersection_model_num_list.append(output_pool[str(j)]["tem_intersection_model_num_list"][k])
                        #             tem_intersection_gene_num_list.append(output_pool[str(j)]["tem_intersection_gene_num_list"][k])
                        #             tem_intersection_gene_list.append(output_pool[str(j)]["tem_intersection_gene_list"][k])

                    else:
                        tem_intersection_model_list = []
                        tem_intersection_model_num_list = []
                        tem_intersection_gene_num_list = []
                        tem_intersection_gene_list = []
                        for i in range(len(current_intersection_model_list)):
                            for j in range(len(all_shared_genes_list)):
                                if all_model_combination_list[j] not in ast.literal_eval(current_intersection_model_list[i]):
                                    tem_name_list = ast.literal_eval(current_intersection_model_list[i]).copy()
                                    tem_name_list.append(all_model_combination_list[j])
                                    tem_name_list.sort()
                                    if tem_name_list not in tem_intersection_model_list:
                                        tem_list = []
                                        tem_list.append(ast.literal_eval(current_intersection_gene_list[i]))
                                        tem_list.append(all_shared_genes_list[j])
                                        share_gene = set.intersection(*map(set, tem_list))
                                        num = len(share_gene)
                                        if num > 1:
                                            tem_share_gene_list = list(share_gene)
                                            tem_share_gene_list.sort()
                                            tem_intersection_gene_num_list.append(num)
                                            tem_intersection_model_list.append(tem_name_list)
                                            tem_intersection_model_num_list.append(intersect_model_num)
                                            tem_intersection_gene_list.append(tem_share_gene_list)




                    print(str(intersect_model_num) + " intersection finished")
                    print(len(tem_intersection_gene_num_list))
                    print("time: " + str(time.time() - t))



                    if len(tem_intersection_gene_num_list) == 0:
                        break

                    tem_df_2 = pd.DataFrame({"shared_gene_num": tem_intersection_gene_num_list, "model_num": tem_intersection_model_num_list,
                                             "shared_genes": tem_intersection_gene_list, "model_combination": tem_intersection_model_list}).sort_values(['model_num', 'shared_gene_num'], ascending=[False, False])

                    tem_str = report_folder_path + "/tem/" + str(intersect_model_num) + "_" + gene_type + "_gene_intersection_relationship.csv"
                    tem_df_2.to_csv(tem_str, index=False)
                    last_intersect_model_num = intersect_model_num


                    t = time.time()
                    print(str(intersect_model_num) + " intersection dictionary generation")

                    num_cores = int(mp.cpu_count())

                    if len(current_intersection_model_list) >= num_cores:
                        separation_length = int(len(current_intersection_model_list) / num_cores)
                        last_end_index = 0
                        tem_list_dict = {}
                        for i in range(num_cores-1):
                            tem_list_dict["model_list_" + str(i)] = current_intersection_model_list[last_end_index:(last_end_index + separation_length)]
                            last_end_index = last_end_index + separation_length

                        tem_list_dict["model_list_" + str(i+1)] = current_intersection_model_list[last_end_index:]

                        pool = mp.Pool(num_cores)
                        process_pool = {}
                        tem_record_dict = {}
                        for i in range(num_cores):
                            process_pool["process_"+str(i)] = pool.apply_async(func_3, args=(tem_list_dict["model_list_" + str(i)], tem_intersection_model_list))

                        pool.close()
                        pool.join()

                        for key in process_pool.keys():
                            tem_record_dict.update(process_pool[key].get())


                    else:
                        tem_record_dict = {}

                        for i in range(len(current_intersection_model_list)):
                            tem_name = ""
                            transform_current_intersection_model = ast.literal_eval(current_intersection_model_list[i])
                            for name in transform_current_intersection_model:
                                tem_name = tem_name + name + ","
                            if tem_name not in tem_record_dict:
                                tem_record_dict[tem_name] = []
                            for j in range(len(tem_intersection_model_list)):
                                flag = 0
                                for k in range(len(transform_current_intersection_model)):
                                    if transform_current_intersection_model[k] in tem_intersection_model_list[j]:
                                        flag = flag + 1
                                if flag == len(transform_current_intersection_model):
                                    for model_name in tem_intersection_model_list[j]:
                                        if model_name not in transform_current_intersection_model:
                                            tem_record_dict[tem_name].append(model_name)
                            if len(tem_record_dict[tem_name]) == 0:
                                tem_record_dict.pop(tem_name, None)

                    print(str(intersect_model_num) + " intersection dictionary generation finished")
                    print("time: " + str(time.time() - t))

                    tem_intersection_model_list = []
                    tem_intersection_model_num_list = []
                    tem_intersection_gene_num_list = []
                    tem_intersection_gene_list = []





                print("redundancy check")

                df_combine_list = []
                for i in range(2, last_intersect_model_num+1):
                    df_combine_list.append(pd.read_csv(report_folder_path + "/tem/" + str(i) + "_" + gene_type + "_gene_intersection_relationship.csv"))

                combine_df = pd.concat(df_combine_list, ignore_index=True)
                combine_df_sorted = combine_df.sort_values(['model_num', 'shared_gene_num'], ascending=[False, False])

                # intersection_model_list = []
                # intersection_model_num_list = []
                # intersection_gene_num_list = []
                # intersection_gene_list = []
                #
                # current_intersection_model_list = tem_intersection_model_list.copy()
                # current_intersection_model_num_list = tem_intersection_model_num_list.copy()
                # current_intersection_gene_num_list = tem_intersection_gene_num_list.copy()
                # current_intersection_gene_list = tem_intersection_gene_list.copy()
                #
                # tem_intersection_model_list = []
                # tem_intersection_model_num_list = []
                # tem_intersection_gene_num_list = []
                # tem_intersection_gene_list = []


                model_combination_list = combine_df_sorted["model_combination"].tolist()
                model_num_list = combine_df_sorted["model_num"].tolist()
                shared_genes_list = combine_df_sorted["shared_genes"].tolist()
                shared_gene_num_list = combine_df_sorted["shared_gene_num"].tolist()

                # tem_df_2 = " "


                intersection_model_list = [ast.literal_eval(model_combination_list[0])]
                intersection_model_num_list = [model_num_list[0]]
                intersection_gene_num_list = [shared_gene_num_list[0]]
                intersection_gene_list = [ast.literal_eval(shared_genes_list[0])]

                for i in range(1, len(model_combination_list)):
                    flag = 0
                    for j in range(len(intersection_model_list)):
                        if set(ast.literal_eval(model_combination_list[i])).issubset(intersection_model_list[j]) and (ast.literal_eval(shared_genes_list[i]) == intersection_gene_list[j]):
                            flag = 1
                            break
                    if flag == 0:
                        intersection_model_list.append(ast.literal_eval(model_combination_list[i]))
                        intersection_model_num_list.append(model_num_list[i])
                        intersection_gene_list.append(ast.literal_eval(shared_genes_list[i]))
                        intersection_gene_num_list.append(shared_gene_num_list[i])

                print(str(last_intersect_model_num) + " intersection finished")

                df_output_4 = pd.DataFrame({"shared_gene_num": intersection_gene_list, "model_num": intersection_model_num_list,
                                            "shared_genes": intersection_gene_num_list,
                                            "model_combination": intersection_model_list})



                tem_str = report_folder_path + "/" + gene_type + "_gene_intersection_relationship.csv"

                df_output_4.to_csv(tem_str, index = False)




        # index_list = [*range(len(record_dict["features"]))]
        #
        # intersection_model_list = []
        # intersection_model_num_list = []
        # intersection_gene_num_list = []
        # intersection_gene_list = []
        #
        # for i in range(2, len(record_dict["features"]) + 1):
        #     combination_list = list(combinations(index_list, i))
        #     for combination in combination_list:
        #         tem_list = []
        #         tem_name_list = []
        #         for index in combination:
        #             tem_list.append(record_dict["features"][index])
        #             tem_name_list.append(record_dict["name"][index])
        #         share_gene = set.intersection(*map(set, tem_list))
        #         num = len(share_gene)
        #         if num > 0:
        #             intersection_gene_num_list.append(num)
        #             intersection_model_list.append(tem_name_list)
        #             intersection_model_num_list.append(i)
        #             intersection_gene_list.append(share_gene)
        #
        # tem_df_1 = pd.DataFrame({"shared_gene_num": intersection_gene_num_list, "model_num": intersection_model_num_list, "shared_genes": intersection_gene_list,
        #      "model_combination": intersection_model_list})
        #
        # tem_df_2 = tem_df_1.sort_values(['model_num', 'shared_gene_num'], ascending=[False, False])
        #
        # model_combination_list = tem_df_2["model_combination"].tolist()
        # model_num_list = tem_df_2["model_num"].tolist()
        # shared_genes_list = tem_df_2["shared_genes"].tolist()
        # shared_gene_num_list = tem_df_2["shared_gene_num"].tolist()
        #
        # model_combination_list_update = [model_combination_list[0]]
        # model_num_list_update = [model_num_list[0]]
        # shared_genes_list_update = [shared_genes_list[0]]
        # shared_gene_num_list_update = [shared_gene_num_list[0]]
        #
        # for i in range(1, len(model_combination_list)):
        #     flag = 0
        #     for j in range(len(model_combination_list_update)):
        #         if set(model_combination_list[i]).issubset(model_combination_list_update[j]) and (shared_genes_list[i] == shared_genes_list_update[j]):
        #             flag = 1
        #             break
        #     if flag == 0:
        #         model_combination_list_update.append(model_combination_list[i])
        #         model_num_list_update.append(model_num_list[i])
        #         shared_genes_list_update.append(shared_genes_list[i])
        #         shared_gene_num_list_update.append(shared_gene_num_list[i])
        #
        # df_output_4 = pd.DataFrame({"shared_gene_num": shared_gene_num_list_update, "model_num": model_num_list_update,
        #                             "shared_genes": shared_genes_list_update,
        #                             "model_combination": model_combination_list_update,  })
        #
        #
        #
        # tem_str = report_folder_path + "/" + gene_type + "_gene_intersection_relationship.csv"
        #
        # df_output_4.to_csv(tem_str, index = False)




