import os
import re
import pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold
import shutil
import subprocess
import time
from sklearn.preprocessing import MinMaxScaler
############################################################
### parameter need to set by user:

biomarker_target_gene_type = "protein_coding" ## "protein_coding" or "lincRNA"

############################################################

### Internal Code from here, further edit not required

submit_script_string = '''
#! /bin/bash

#BSUB -L /bin/bash
#BSUB -q short
#BSUB -n 20 -W 8:00
#BSUB -o myjob.out
#BSUB -e myjob.err
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=5000]

'''

standard_file_path = "./output/summary_step1.txt"
datatype = "BMGD_processed"
folder_path = "."

gene_type_list = [biomarker_target_gene_type]

DE_type = "wilcox"

microarray_logtransform = str(0) ## if data has been log transformed: 0, if not or not sure: 1

preprocessing_python_script_path = "./src/main_CV_20220725_add_microarray_sliverman_combined_DE_revision.py"

for gene_type in gene_type_list:
    if gene_type == "protein_coding":
        FPKM_standard_file_path = standard_file_path + ".FPKM.ProteinCoding"
        ReadCount_standard_file_path = standard_file_path + ".ReadCount.ProteinCoding"
        ReadCount_overall_file_path = standard_file_path + ".ReadCount"
    elif gene_type == "lincRNA":
        FPKM_standard_file_path = standard_file_path + ".FPKM.lincRNA"
        ReadCount_standard_file_path = standard_file_path + ".ReadCount.lincRNA"
        ReadCount_overall_file_path = standard_file_path + ".ReadCount"
    elif gene_type == "microarray":
        FPKM_standard_file_path = standard_file_path + ".microarray"
        ReadCount_standard_file_path = standard_file_path + ".microarray"

    if datatype == "BMGD_processed":
        standard_FPKM_df = pd.read_csv(FPKM_standard_file_path, sep="\t", index_col=0)
        standard_ReadCount_df = pd.read_csv(ReadCount_standard_file_path, sep="\t", index_col=0)
        sample_name_list = standard_FPKM_df.columns.tolist()


    partition_folder_path = os.getcwd()

    tem_output_directory = partition_folder_path + "/output"
    os.mkdir(tem_output_directory)

    os.chdir(tem_output_directory)

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


    bash_string_1 = 'python {preprocessing_python_script_path} --gene_type_list {gene_type_list} --only_filter_num {only_filter_num} --directory {directory} --DE_type {DE_type} --filter_requirment_table_path {filter_requirment_table_path}'

    bash_string_1_R = bash_string_1.format(preprocessing_python_script_path=preprocessing_python_script_path,
                                           gene_type_list=' '.join(str(e) for e in [gene_type]),
                                           only_filter_num=0, directory=partition_folder_path,
                                           DE_type=DE_type, filter_requirment_table_path="./")

    bash_string = submit_script_string + bash_string_1_R

    with open(partition_folder_path + "/submit.sh", 'w') as rsh:
        rsh.write(bash_string)

    os.chdir(partition_folder_path)
    system_output = os.popen('bsub < ./submit.sh').read()

