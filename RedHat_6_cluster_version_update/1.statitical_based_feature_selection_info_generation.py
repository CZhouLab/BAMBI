import os
import re
import pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold
import shutil
import subprocess
import time
from sklearn.preprocessing import MinMaxScaler
import argparse
############################################################
### parameter need to set by user:

# biomarker_target_gene_type = "protein_coding" ## "protein_coding" or "lincRNA"

############################################################
def statitical_based_feature_selection_info_generation(biomarker_target_gene_type="protein_coding", RNASeq_FPKM_table_path=None,
                                                       RNASeq_ReadCount_table_path=None, microarray_table_path=None):

    # submit_script_string = '''
    # #! /bin/bash
    #
    # #BSUB -L /bin/bash
    # #BSUB -q short
    # #BSUB -n 20 -W 8:00
    # #BSUB -o myjob.out
    # #BSUB -e myjob.err
    # #BSUB -R span[hosts=1]
    # #BSUB -R rusage[mem=5000]
    #
    # '''



    # datatype = "BMGD_processed" ##  "BMGD_processed" "original_microarray"
    ### original_microarray requires standard file
    ### BMGD_processed requies XXX/output/summary_step1.txt"

    directory = os.getcwd()
    if RNASeq_FPKM_table_path is None and RNASeq_ReadCount_table_path is None and microarray_table_path is None:
        datatype = "BMGD_processed"
        standard_file_path = directory + "/output/summary_step1.txt"
    elif RNASeq_FPKM_table_path is None and RNASeq_ReadCount_table_path is None and microarray_table_path is not None:
        datatype = "original_microarray"
    elif RNASeq_FPKM_table_path is not None and RNASeq_ReadCount_table_path is not None and microarray_table_path is None:
        datatype = "original_RNA_Seq"






    gene_type_list = [biomarker_target_gene_type] # protein_coding lincRNA

    if biomarker_target_gene_type in ["protein_coding", "lincRNA"]:
        DE_type = "wilcox" ## "wilcox", "DEseq2", "limma"
    elif biomarker_target_gene_type == "microarray":
        DE_type = "limma"

    microarray_logtransform = str(0) ## if data has been log transformed: 0, if not or not sure: 1

    output_directory = directory + "/output"
    src_path = directory + "/src"

    preprocessing_python_script_path = src_path + "/main_CV_20220725_add_microarray_sliverman_combined_DE_revision.py"

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    ML_path = output_directory + "/ML"
    if not os.path.isdir(ML_path):
        os.mkdir(ML_path)



    #==================================




    for gene_type in gene_type_list:
        if datatype == "BMGD_processed":
            if gene_type == "protein_coding":
                FPKM_standard_file_path = standard_file_path + ".FPKM.ProteinCoding"
                ReadCount_standard_file_path = standard_file_path + ".ReadCount.ProteinCoding"
                # ReadCount_overall_file_path = standard_file_path + ".ReadCount"
            elif gene_type == "lincRNA":
                FPKM_standard_file_path = standard_file_path + ".FPKM.lincRNA"
                ReadCount_standard_file_path = standard_file_path + ".ReadCount.lincRNA"
                # ReadCount_overall_file_path = standard_file_path + ".ReadCount"
            # elif gene_type == "microarray":
            #     FPKM_standard_file_path = standard_file_path + ".microarray"
            #     ReadCount_standard_file_path = standard_file_path + ".microarray"
            standard_FPKM_df = pd.read_csv(FPKM_standard_file_path, sep="\t", index_col=0)
            standard_ReadCount_df = pd.read_csv(ReadCount_standard_file_path, sep="\t", index_col=0)
            # sample_name_list = standard_FPKM_df.columns.tolist()

        if datatype == "original_RNA_Seq":
            standard_FPKM_df = pd.read_csv(RNASeq_FPKM_table_path, sep="\t", index_col=0)
            standard_ReadCount_df = pd.read_csv(RNASeq_ReadCount_table_path, sep="\t", index_col=0)

        elif datatype == "original_microarray":
            run_Pipeline_step_2_2_6_para_1 = "./src/limmaDE/limma_Preprocessing.R"

            run_Pipeline_step_2_2_6_para_2 = microarray_table_path
            run_Pipeline_step_2_2_6_para_3 = microarray_logtransform
            run_Pipeline_step_2_2_6_para_4 = output_directory + "/summary_step1.txt.microarray"
            run_Pipeline_step_2_2_6 = subprocess.Popen(['Rscript', run_Pipeline_step_2_2_6_para_1,
                                                        run_Pipeline_step_2_2_6_para_2, run_Pipeline_step_2_2_6_para_3,
                                                        run_Pipeline_step_2_2_6_para_4])
            run_Pipeline_step_2_2_6.communicate()



            standard_FPKM_df = pd.read_csv(output_directory + "/summary_step1.txt.microarray", sep="\t", index_col=0)
            standard_ReadCount_df = pd.read_csv(output_directory + "/summary_step1.txt.microarray", sep="\t", index_col=0)
            # sample_name_list = standard_FPKM_df.columns.tolist()


        partition_folder_path = os.getcwd()

        tem_output_directory = partition_folder_path + "/output"
        os.mkdir(tem_output_directory)

        os.chdir(tem_output_directory)
        if datatype == "original_RNA_Seq":
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


            FPKM_output_train_path = tem_output_directory + "/" + FPKM_name

            standard_FPKM_df.to_csv(FPKM_output_train_path, sep="\t")

            ReadCount_output_train_path = tem_output_directory + "/" + ReadCount_name

            standard_ReadCount_df.to_csv(ReadCount_output_train_path, sep="\t")


        # if gene_type in ["lincRNA", "protein_coding"]:
        #     shutil.copy2(ReadCount_overall_file_path, tem_output_directory + "/summary_step1.txt.ReadCount")



        # bash_string_1 = 'python {preprocessing_python_script_path} --gene_type_list {gene_type_list} --only_filter_num {only_filter_num} --directory {directory} --DE_type {DE_type} --filter_requirment_table_path {filter_requirment_table_path}'
        #
        # bash_string_1_R = bash_string_1.format(preprocessing_python_script_path=preprocessing_python_script_path,
        #                                        gene_type_list=' '.join(str(e) for e in [gene_type]),
        #                                        only_filter_num=0, directory=partition_folder_path,
        #                                        DE_type=DE_type, filter_requirment_table_path="./")
        #
        # bash_string = submit_script_string + bash_string_1_R
        #
        # with open(partition_folder_path + "/submit.sh", 'w') as rsh:
        #     rsh.write(bash_string)
        #
        # os.chdir(partition_folder_path)
        # system_output = os.popen('bsub < ./submit.sh').read()

        # bash_string_1 = 'python {preprocessing_python_script_path} --gene_type_list {gene_type_list} --only_filter_num {only_filter_num} --directory {directory} --DE_type {DE_type} --filter_requirment_table_path {filter_requirment_table_path}'
        #
        # bash_string_1_R = bash_string_1.format(preprocessing_python_script_path=preprocessing_python_script_path,
        #                                        gene_type_list=' '.join(str(e) for e in [gene_type]),
        #                                        only_filter_num=0, directory=partition_folder_path,
        #                                        DE_type=DE_type, filter_requirment_table_path="./")

        # bash_string = submit_script_string + bash_string_1_R
        #
        # with open(partition_folder_path + "/submit.sh", 'w') as rsh:
        #     rsh.write(bash_string)
        #
        # os.chdir(partition_folder_path)
        # system_output = os.popen('bsub < ./submit.sh').read()

        run_Pipeline_tem = subprocess.Popen(['python', preprocessing_python_script_path,
                                             '--gene_type_list', ' '.join(str(e) for e in [gene_type]),
                                             '--only_filter_num', str(0),
                                             '--directory', partition_folder_path,
                                             '--DE_type', DE_type,
                                             '--filter_requirment_table_path', "./"])

        run_Pipeline_tem.communicate()


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--biomarker_target_gene_type", choices=["protein_coding", "lincRNA", "microarray"], default="protein_coding", required=True)
    parser.add_argument("--RNASeq_FPKM_table_path", default=None, type=str, required=False)
    parser.add_argument("--RNASeq_ReadCount_table_path", default=None, type=str, required=False)
    parser.add_argument("--microarray_table_path", default=None, type=str, required=False)

    args = parser.parse_args()

    statitical_based_feature_selection_info_generation(biomarker_target_gene_type=args.biomarker_target_gene_type,
                                                       RNASeq_FPKM_table_path=args.RNASeq_FPKM_table_path,
                                                       RNASeq_ReadCount_table_path=args.RNASeq_ReadCount_table_path,
                                                       microarray_table_path=args.microarray_table_path
                                                       )
    # stab_selection_main(folder_path=args.folder_path, file_name=args.train_file_name, testfile_name=args.test_file_name, threshold=args.threshold)



if __name__ == '__main__':
    if not sys.warnoptions:
        warnings.simplefilter("ignore")
        os.environ["PYTHONWARNINGS"] = "ignore"  # Also affect subprocesses
    main()
