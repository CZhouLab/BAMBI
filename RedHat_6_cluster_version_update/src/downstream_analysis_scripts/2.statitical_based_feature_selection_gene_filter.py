import os
import re
import pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold
import shutil
import subprocess
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
from statistics import median
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import sys
import argparse

############################################################
### parameter need to set by user:

# biomarker_target_gene_type = "protein_coding" ## "protein_coding" or "lincRNA"
#
# target_pvalue_type = "padj" # "pvalue" or "padj"
# target_pvalue_threshold = 0.05
# target_foldchange_threshold = 1.0
# target_maxmin_remove_threshold = 1.0
# target_overlap_area_threshold = 0.01


############################################################

def statitical_based_feature_selection_gene_filter(directory, biomarker_target_gene_type="protein_coding", target_pvalue_type = "padj",
                                                       target_pvalue_threshold = 0.05, target_foldchange_threshold = 1.0,
                                                   target_maxmin_remove_threshold = 1.0, target_overlap_area_threshold = 0.01, HPC_parallel = False):

    folder_path = directory

    gene_type_list = [biomarker_target_gene_type] # protein_coding lincRNA microarray

    outlier_remove = True
    #==================================

    pvalue_type_list = [target_pvalue_type]
    pvalue_threshold_list = [target_pvalue_threshold]
    foldchange_threshold_list = [target_foldchange_threshold]
    maxmin_remove_threshold_list = [target_maxmin_remove_threshold]
    overlap_area_threshold_list = [target_overlap_area_threshold]

    output_dict = {}

    # output_dict["partitions"] = []
    # if outlier_remove:
    #     pvalue_type_list = pvalue_type_list * 2
    #     pvalue_threshold_list = pvalue_threshold_list * 2
    #     foldchange_threshold_list = foldchange_threshold_list * 2
    #     overlap_area_threshold_list = overlap_area_threshold_list * 2
    #     maxmin_remove_threshold_list = maxmin_remove_threshold_list * 2


    output_dict["pvalue_type"] = pvalue_type_list
    output_dict["pvalue_threshold"] = pvalue_threshold_list
    output_dict["foldchange_threshold"] = foldchange_threshold_list
    output_dict["maxmin_remove_threshold"] = maxmin_remove_threshold_list
    output_dict["overlap_area_threshold"] = overlap_area_threshold_list

    output_dict["num_Original"] = []
    output_dict["num_after_DE_pvalue_filter"] = []
    output_dict["num_after_foldchange_filter"] = []
    output_dict["num_after_maximum_value_filter"] = []
    output_dict["num_after_overlap_area_filter"] = []

    folder_list = []
    outlier_remove_list = []

    # general_output_folder = folder_path + "/" + name
    #
    # if os.path.isdir(general_output_folder):
    #     shutil.rmtree(general_output_folder)
    # os.mkdir(general_output_folder)
    #
    # os.system("cp -R " + folder_path + "/output " + general_output_folder + "/output")
    #
    # folder_list.append(general_output_folder)
    # outlier_remove_list.append(0)
    #
    # if outlier_remove:
    #     outlier_remove_output_folder = folder_path + "/" + name + "_outlier_remove"
    #
    #     if os.path.isdir(outlier_remove_output_folder):
    #         shutil.rmtree(outlier_remove_output_folder)
    #     os.mkdir(outlier_remove_output_folder)
    #
    #     os.system("cp -R " + folder_path + "/output " + outlier_remove_output_folder + "/output")
    #
    #     folder_list.append(outlier_remove_output_folder)
    #     outlier_remove_list.append(1)

    folder_list = [folder_path]
    if outlier_remove:
        outlier_remove_list = [1]
    else:
        outlier_remove_list = [0]


    # for i in range(len(outlier_remove_list)):
    i = 0
    if outlier_remove_list[i] == 1:
        outlier_info_folder_path = folder_list[i] + "/outlier_remove_info"

        if os.path.isdir(outlier_info_folder_path):
            shutil.rmtree(outlier_info_folder_path)
        os.mkdir(outlier_info_folder_path)

        outlier_info_pic_folder_path = outlier_info_folder_path + "/pic"
        os.mkdir(outlier_info_pic_folder_path)

        outlier_record_dict = {}
        # outlier_record_dict["partition"] = []
        outlier_record_dict["target_gene"] = []
        outlier_record_dict["C_upper_bound"] = []
        outlier_record_dict["T_upper_bound"] = []

        if "microarray" in gene_type_list:
            outlier_record_dict["C_lower_bound"] = []
            outlier_record_dict["T_lower_bound"] = []
            outlier_record_dict["outlier_remove_upper_threshold"] = []
            outlier_record_dict["outlier_remove_lower_threshold"] = []
        else:
            outlier_record_dict["outlier_remove_threshold"] = []





    output_directory = folder_list[i] + "/output"

    Filter_num_record_dict = {}
    Filter_num_record_dict["step"] = []
    Filter_num_record_dict["Gene_number"] = []

    for gene_type in gene_type_list:
        if gene_type == "protein_coding":
            output_name = "ProteinCoding"
            FPKM_name = "summary_step1.txt.FPKM.ProteinCoding"
            ReadCount_name = "summary_step1.txt.ReadCount.ProteinCoding"
            # FPKM_standard_file_path = standard_file_path + ".FPKM.ProteinCoding"
        elif gene_type == "lincRNA":
            output_name = "lincRNA"
            FPKM_name = "summary_step1.txt.FPKM.lincRNA"
            ReadCount_name = "summary_step1.txt.ReadCount.lincRNA"
            # FPKM_standard_file_path = standard_file_path + ".txt.FPKM.lincRNA"
        elif gene_type == "microarray":
            output_name = "microarray"
            FPKM_name = "summary_step1.txt.microarray"
            ReadCount_name = "summary_step1.txt.microarray"
            # FPKM_standard_file_path = standard_file_path + ".microarray"

        ## filter by DE p-value

        info_df = pd.read_csv(output_directory + '/' + gene_type + '_Gene_info.csv', sep='\t')

        original_gene_size = len(pd.read_csv(output_directory + '/' + FPKM_name, sep="\t", index_col=0))

        Filter_num_record_dict["step"].append(gene_type + " Original")
        Filter_num_record_dict["Gene_number"].append(original_gene_size)
        output_dict["num_Original"].append(original_gene_size)


        DE_filtered_info_df = info_df[(info_df[pvalue_type_list[i]] <= float(pvalue_threshold_list[i]))]

        Filter_num_record_dict["step"].append(gene_type + " after DE pvalue filter")
        Filter_num_record_dict["Gene_number"].append(len(DE_filtered_info_df))
        output_dict["num_after_DE_pvalue_filter"].append(len(DE_filtered_info_df))

        ## filter by foldchange

        FC_filtered_info_df = DE_filtered_info_df[
            (DE_filtered_info_df['abs(log2(fold-change))'] >= float(foldchange_threshold_list[i]))]

        Filter_num_record_dict["step"].append(gene_type + " after foldchange filter")
        Filter_num_record_dict["Gene_number"].append(len(FC_filtered_info_df))
        output_dict["num_after_foldchange_filter"].append(len(FC_filtered_info_df))

        ## filter by maximum value

        max_val_filtered_info_df = FC_filtered_info_df[
            FC_filtered_info_df["max_val"] >= maxmin_remove_threshold_list[i]]

        Filter_num_record_dict["step"].append(gene_type + " after maximum value filter")
        Filter_num_record_dict["Gene_number"].append(len(max_val_filtered_info_df))
        output_dict["num_after_maximum_value_filter"].append(len(max_val_filtered_info_df))

        ## filter by overlap area

        overlap_filtered_info_df = max_val_filtered_info_df[
            max_val_filtered_info_df["overlap_area"] <= overlap_area_threshold_list[i]]

        Filter_num_record_dict["step"].append(gene_type + " after overlap area filter")
        Filter_num_record_dict["Gene_number"].append(len(overlap_filtered_info_df))

        output_dict["num_after_overlap_area_filter"].append(len(overlap_filtered_info_df))

        if os.path.exists(output_directory + '/Filtered_' + gene_type + '_Gene_info.csv'):
            os.remove(output_directory + '/Filtered_' + gene_type + '_Gene_info.csv')

        overlap_filtered_info_df.to_csv(output_directory + '/Filtered_' + gene_type + '_Gene_info.csv', sep="\t",
                                        index=False)

        filtered_Gene_list = overlap_filtered_info_df["Gene"].tolist()




        # standard_FPKM_df_T = pd.read_csv(output_directory + "/" + FPKM_name, sep="\t", index_col=0).T
        #
        # standard_FPKM_df_T_transform = (standard_FPKM_df_T) / (standard_FPKM_df_T.max(axis=0))
        #
        # standard_FPKM_df_transform = standard_FPKM_df_T_transform.T
        #
        # FPKM_df = standard_FPKM_df_transform.reset_index()



        # FPKM_df = pd.read_csv(output_directory + "/" + FPKM_name, sep="\t")

        def IQR_remove_calculation(target_gene, expression_df):
            expression_list = expression_df[target_gene].tolist()
            Q1 = np.percentile(expression_list, 25, interpolation='midpoint')
            Q3 = np.percentile(expression_list, 75, interpolation='midpoint')
            IQR = Q3 - Q1

            upper_bound = Q3 + 1.5 * IQR
            lower_bound = Q1 - 1.5 * IQR

            return upper_bound, lower_bound


        if outlier_remove_list[i] == 1:
            pic_flag = 0
            ori_FPKM_df = pd.read_csv(output_directory + "/" + FPKM_name, sep="\t", index_col=0).T[filtered_Gene_list].copy()

            label_list = []
            for tem in ori_FPKM_df.index.tolist():
                if tem.split("_")[-1] == "C":
                    label_list.append(int(0))
                elif tem.split("_")[-1] == "T":
                    label_list.append(int(1))
            ori_FPKM_df["y"] = label_list

            ori_FPKM_C_sub_df = ori_FPKM_df[ori_FPKM_df["y"] == 0].copy()
            ori_FPKM_T_sub_df = ori_FPKM_df[ori_FPKM_df["y"] == 1].copy()

            FPKM_df = pd.read_csv(output_directory + "/" + FPKM_name, sep="\t", index_col=0).T[filtered_Gene_list].copy()
            # FPKM_df_test = pd.read_csv(output_directory + "/" + FPKM_name + ".test", sep="\t", index_col=0).T[filtered_Gene_list].copy()


            for target_gene in filtered_Gene_list:
                tem_C_upper_bound, tem_C_lower_bound = IQR_remove_calculation(target_gene=target_gene, expression_df=ori_FPKM_C_sub_df)
                tem_T_upper_bound,  tem_T_lower_bound = IQR_remove_calculation(target_gene=target_gene, expression_df=ori_FPKM_T_sub_df)

                # outlier_record_dict["partition"].append("CV_partition_" + str(partition))
                outlier_record_dict["target_gene"].append(target_gene)
                if gene_type == "microarray":
                    outlier_record_dict["C_upper_bound"].append(tem_C_upper_bound)
                    outlier_record_dict["T_upper_bound"].append(tem_T_upper_bound)
                    outlier_record_dict["C_lower_bound"].append(tem_C_lower_bound)
                    outlier_record_dict["T_lower_bound"].append(tem_T_lower_bound)
                    outlier_record_dict["outlier_remove_upper_threshold"].append(max([tem_C_upper_bound, tem_T_upper_bound, tem_C_lower_bound, tem_T_lower_bound]))
                    outlier_record_dict["outlier_remove_lower_threshold"].append(min([tem_C_upper_bound, tem_T_upper_bound, tem_C_lower_bound, tem_T_lower_bound]))


                else:
                    outlier_record_dict["C_upper_bound"].append(tem_C_upper_bound)
                    outlier_record_dict["T_upper_bound"].append(tem_T_upper_bound)
                    outlier_record_dict["outlier_remove_threshold"].append(max([tem_C_upper_bound, tem_T_upper_bound]))




                Control_num = 1
                Patient_num = 1
                Color_list = []

                median_C = median(ori_FPKM_C_sub_df[target_gene].tolist())
                median_T = median(ori_FPKM_T_sub_df[target_gene].tolist())

                if median_C >= median_T:
                    ascending_str = False
                else:
                    ascending_str = True


                update_df_C = ori_FPKM_C_sub_df.sort_values([target_gene], ascending=[ascending_str])
                update_df_T = ori_FPKM_T_sub_df.sort_values([target_gene], ascending=[ascending_str])
                y_list = update_df_C["y"].tolist() + update_df_T["y"].tolist()
                expression_list = update_df_C[target_gene].tolist() + update_df_T[target_gene].tolist()


                tem_label_list = []
                for label in y_list:
                    if label == 0:
                        tem_label_list.append("Control_" + str(Control_num))
                        Control_num += 1
                        Color_list.append("#1f77b4")
                    elif label == 1:
                        tem_label_list.append("Patient_" + str(Patient_num))
                        Patient_num += 1
                        Color_list.append("Orange")
                fig = plt.figure(figsize=(14, 7), dpi=100)
                plt.bar(tem_label_list, expression_list, color=Color_list)
                plt.xticks(color="white")



                if tem_T_upper_bound >= tem_C_upper_bound:
                    outlier_remove_threshold = tem_T_upper_bound
                    custom_lines = [Line2D([0], [0], color="#1F77B4", lw=4, label='Control'),
                                    Line2D([0], [0], color="orange", lw=4, label='Patient'),
                                    Line2D([0], [0], color="purple", lw=4, label='higher_group_boundary_T'),
                                    Line2D([0], [0], color="green", lw=4, label='lower_group_boundary_C')]
                    plt.axhline(y=tem_T_upper_bound, color='purple', linestyle='--')
                    plt.axhline(y=tem_C_upper_bound, color='green', linestyle='--')

                    if gene_type == "microarray":
                        plt.axhline(y=tem_T_lower_bound, color='purple', linestyle='--')
                        plt.axhline(y=tem_C_lower_bound, color='green', linestyle='--')


                else:
                    outlier_remove_threshold = tem_C_upper_bound
                    custom_lines = [Line2D([0], [0], color="#1F77B4", lw=4, label='Control'),
                                    Line2D([0], [0], color="orange", lw=4, label='Patient'),
                                    Line2D([0], [0], color="purple", lw=4, label='higher_group_boundary_C'),
                                    Line2D([0], [0], color="green", lw=4, label='lower_group_boundary_T')]

                    plt.axhline(y=tem_C_upper_bound, color='purple', linestyle='--')
                    plt.axhline(y=tem_T_upper_bound, color='green', linestyle='--')
                    if gene_type == "microarray":
                        plt.axhline(y=tem_C_lower_bound, color='purple', linestyle='--')
                        plt.axhline(y=tem_T_lower_bound, color='green', linestyle='--')



                plt.legend(handles=custom_lines)
                plt.title("FPKM distribution: " + target_gene, fontsize=20)
                plt.savefig(outlier_info_pic_folder_path + "/" + str(pic_flag) + ".jpg", bbox_inches='tight')
                pic_flag = pic_flag + 1

                if gene_type == "microarray":
                    microarray_outlier_remove_upper_threshold = max([tem_C_upper_bound, tem_T_upper_bound, tem_C_lower_bound, tem_T_lower_bound])
                    microarray_outlier_remove_lower_threshold = min([tem_C_upper_bound, tem_T_upper_bound, tem_C_lower_bound, tem_T_lower_bound])

                # XFPKM_df = pd.read_csv(output_directory + "/" + FPKM_name, sep="\t", index_col=0).T[filtered_Gene_list].copy()
                # XFPKM_df_test = pd.read_csv(output_directory + "/" + FPKM_name + ".test", sep="\t", index_col=0).T[filtered_Gene_list].copy()

                if gene_type == "microarray":
                    tem_update_target_gene_list = []
                    for tem_target_gene_expression in FPKM_df[target_gene].tolist():
                        if tem_target_gene_expression > microarray_outlier_remove_upper_threshold:
                            tem_update_target_gene_list.append(microarray_outlier_remove_upper_threshold)
                        elif tem_target_gene_expression < microarray_outlier_remove_lower_threshold:
                            tem_update_target_gene_list.append(microarray_outlier_remove_lower_threshold)
                        else:
                            tem_update_target_gene_list.append(tem_target_gene_expression)
                    FPKM_df[target_gene] = tem_update_target_gene_list

                    # tem_update_target_gene_list = []
                    # for tem_target_gene_expression in FPKM_df_test[target_gene].tolist():
                    #     if tem_target_gene_expression > microarray_outlier_remove_upper_threshold:
                    #         tem_update_target_gene_list.append(microarray_outlier_remove_upper_threshold)
                    #     elif tem_target_gene_expression < microarray_outlier_remove_lower_threshold:
                    #         tem_update_target_gene_list.append(microarray_outlier_remove_lower_threshold)
                    #     else:
                    #         tem_update_target_gene_list.append(tem_target_gene_expression)
                    # FPKM_df_test[target_gene] = tem_update_target_gene_list


                else:
                    tem_update_target_gene_list = []
                    for tem_target_gene_expression in FPKM_df[target_gene].tolist():
                        if tem_target_gene_expression > outlier_remove_threshold:
                            tem_update_target_gene_list.append(outlier_remove_threshold)
                        else:
                            tem_update_target_gene_list.append(tem_target_gene_expression)
                    FPKM_df[target_gene] = tem_update_target_gene_list


                    # tem_update_target_gene_list = []
                    # for tem_target_gene_expression in FPKM_df_test[target_gene].tolist():
                    #     if tem_target_gene_expression > outlier_remove_threshold:
                    #         tem_update_target_gene_list.append(outlier_remove_threshold)
                    #     else:
                    #         tem_update_target_gene_list.append(tem_target_gene_expression)
                    # FPKM_df_test[target_gene] = tem_update_target_gene_list



        else:
            FPKM_df = pd.read_csv(output_directory + "/" + FPKM_name, sep="\t", index_col=0).T[filtered_Gene_list].copy()
            # FPKM_df_test = pd.read_csv(output_directory + "/" + FPKM_name + ".test", sep="\t", index_col=0).T[filtered_Gene_list].copy()

        train_sample_list = FPKM_df.index.tolist()
        # test_sample_list = FPKM_df_test.index.tolist()

        # overall_FPKM_df = pd.concat([FPKM_df, FPKM_df_test])

        overall_FPKM_df = FPKM_df.copy()


        if gene_type == "microarray":
            target_max = 1
            target_min = -1
            X_std = (overall_FPKM_df - overall_FPKM_df.min(axis=0)) / (overall_FPKM_df.max(axis=0) - overall_FPKM_df.min(axis=0))
            overall_FPKM_df_transform = X_std * (target_max - target_min) + target_min
        else:
            overall_FPKM_df_transform = (overall_FPKM_df) / (overall_FPKM_df.max(axis=0))

        for tem_num in range(3):
            tem_gene = filtered_Gene_list[tem_num]
            tem_expression_list = overall_FPKM_df_transform[tem_gene].tolist()
            print(tem_gene + " max:" + str(max(tem_expression_list)) + " min:" + str(min(tem_expression_list)))

        Filtered_FPKM_df = overall_FPKM_df_transform.T[train_sample_list].copy()

        # Filtered_FPKM_df = FPKM_df[FPKM_df['Gene'].isin(filtered_Gene_list)]

        if os.path.exists(output_directory + "/Filtered_" + FPKM_name):
            os.remove(output_directory + "/Filtered_" + FPKM_name)

        Filtered_FPKM_df.to_csv(output_directory + "/Filtered_" + FPKM_name, sep="\t")

        # Filtered_FPKM_df_test = overall_FPKM_df_transform.T[test_sample_list].copy()
        # # Filtered_FPKM_df_test = FPKM_df_test[FPKM_df_test['Gene'].isin(filtered_Gene_list)]
        #
        # if os.path.exists(output_directory + "/Filtered_" + FPKM_name + ".test"):
        #     os.remove(output_directory + "/Filtered_" + FPKM_name + ".test")
        #
        # Filtered_FPKM_df_test.to_csv(output_directory + "/Filtered_" + FPKM_name + ".test", sep="\t")



        if gene_type != "microarray":
            ReadCount_df = pd.read_csv(output_directory + "/" + ReadCount_name, sep="\t")

            Filtered_ReadCount_df = ReadCount_df[ReadCount_df['Gene'].isin(filtered_Gene_list)]

            if os.path.exists(output_directory + "/Filtered_" + ReadCount_name):
                os.remove(output_directory + "/Filtered_" + ReadCount_name)

            Filtered_ReadCount_df.to_csv(output_directory + "/Filtered_" + ReadCount_name, sep="\t", index=False)





        def reformat_file_generation(output_directory, FPKM_name, gene_type, testset):
            df_ori_input = pd.read_csv(output_directory + "/Filtered_" + FPKM_name, sep="\t", index_col=0).T

            df_index_name = df_ori_input.index.tolist()

            label_list = []

            for tem in df_index_name:
                if tem.split("_")[-1] == "C":
                    label_list.append(int(0))
                elif tem.split("_")[-1] == "T":
                    label_list.append(int(1))

            # y_dict = {"index": df_index_name,
            #           "y": label_list}
            #
            # df_y = pd.DataFrame(y_dict)
            #
            # df_y.set_index(["index"], inplace=True)

            # revised_df = pd.concat([df_ori_input, df_y], axis=1)

            df_ori_input["y"] = label_list
            if testset:
                reformat_output_path = output_directory + '/Filtered_' + gene_type + '_input.test.txt'
            else:
                reformat_output_path = output_directory + '/Filtered_' + gene_type + '_input.all.txt'

            if os.path.exists(reformat_output_path):
                os.remove(reformat_output_path)

            df_ori_input.to_csv(reformat_output_path, header=True, sep='\t')


        reformat_file_generation(output_directory=output_directory, FPKM_name=FPKM_name, gene_type=gene_type, testset=False)
        # reformat_file_generation(output_directory=output_directory, FPKM_name=FPKM_name + ".test", gene_type=gene_type, testset=True)


    Filter_num_record_output_path = output_directory + '/Filtered_info_record.csv'

    if os.path.exists(Filter_num_record_output_path):
        os.remove(Filter_num_record_output_path)

    pd.DataFrame(Filter_num_record_dict).to_csv(Filter_num_record_output_path, index=False)

    if outlier_remove_list[i] == 1:
        pd.DataFrame(outlier_record_dict).to_csv(outlier_info_folder_path + "/outlier_record.csv", index=False)

    if os.path.exists(folder_path + "/summary_gene_filter_info.csv"):
        os.remove(folder_path + "/summary_gene_filter_info.csv")
    pd.DataFrame(output_dict).to_csv(folder_path + "/summary_gene_filter_info.csv", index=False)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--directory", default=None, type=str, required=True)
    parser.add_argument("--biomarker_target_gene_type", choices=["protein_coding", "lincRNA", "microarray"], default="protein_coding", required=True)
    parser.add_argument("--target_pvalue_type", default=None, type=str, required=True)
    parser.add_argument("--target_pvalue_threshold", default=None, type=float, required=True)
    parser.add_argument("--target_foldchange_threshold", default=None, type=float, required=True)
    parser.add_argument("--target_maxmin_remove_threshold", default=None, type=float, required=True)
    parser.add_argument("--target_overlap_area_threshold", default=None, type=float, required=True)

    args = parser.parse_args()

    statitical_based_feature_selection_gene_filter(directory=args.directory, biomarker_target_gene_type=args.biomarker_target_gene_type,
                                                   target_pvalue_type=args.target_pvalue_type,
                                                   target_pvalue_threshold=args.target_pvalue_threshold, target_foldchange_threshold=args.target_foldchange_threshold,
                                                   target_maxmin_remove_threshold=args.target_maxmin_remove_threshold,
                                                   target_overlap_area_threshold=args.target_overlap_area_threshold)
    # stab_selection_main(folder_path=args.folder_path, file_name=args.train_file_name, testfile_name=args.test_file_name, threshold=args.threshold)



if __name__ == '__main__':
    if not sys.warnoptions:
        warnings.simplefilter("ignore")
        os.environ["PYTHONWARNINGS"] = "ignore"  # Also affect subprocesses
    main()
