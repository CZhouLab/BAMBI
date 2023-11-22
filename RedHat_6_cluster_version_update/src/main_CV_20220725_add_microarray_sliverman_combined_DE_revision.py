
import os
import re
import pandas as pd
import subprocess
import os.path
import time
from KDEpy import FFTKDE
from KDEpy import NaiveKDE
import sys
import math
import random
import numpy as np
import argparse
from sklearn.base import clone
import warnings
############################################################





def job_submission(sample_name, strandness, directory, Known_splice_site, GTF_file, ExonLength_file, sequence_type, R1_input="NA", R2_input="NA", unpaired_input = "NA"):

    output_directory = directory + "/output"

    # bash_string_0 = '''#! /bin/bash
    #
    #     #BSUB -L /bin/bash
    #     #BSUB -q large
    #     #BSUB -n 20 -W 96:00
    #     #BSUB -o myjob.out
    #     #BSUB -e myjob.err
    #     #BSUB -R span[hosts=1]
    #     #BSUB -R rusage[mem=5000]
    #     #BSUB -u Peng.Zhou@umassmed.edu
    #
    #     module load hisat2/2.0.5
    #     module load sambamba/0.6.7
    #     module load python/2.7.9_packages/HTSeq/0.6.1
    #     module load python/2.7.9_packages/pysam/0.8.4
    #
    #     '''
    bash_string_0 = '''#! /bin/bash

        #BSUB -L /bin/bash
        #BSUB -q large
        #BSUB -n 20 -W 96:00
        #BSUB -o myjob.out
        #BSUB -e myjob.err
        #BSUB -R span[hosts=1]
        #BSUB -u Peng.Zhou@umassmed.edu
        
        source activate /home/peng.zhou-umw/shared_env/share_env_20230317_ML_pipeline_raw_sequencing

        '''
    bash_string_1 = 'python {pyscript_path} --sample_name {sample_name} --R1_input {R1_input} --R2_input {R2_input} --unpaired_input {unpaired_input} --strandness {strandness} --directory {directory} --Known_splice_site {Known_splice_site} --GTF_file {GTF_file} --ExonLength_file {ExonLength_file} --sequence_type {sequence_type}'

    pyscript_path= directory + "/cluster_job_submission_FPKM_revision.py"

    bash_string_1_R = bash_string_1.format(pyscript_path=pyscript_path, sample_name=sample_name, R1_input=R1_input, R2_input=R2_input, unpaired_input=unpaired_input, strandness=strandness, directory=directory,
                                           Known_splice_site=Known_splice_site, GTF_file=GTF_file, ExonLength_file=ExonLength_file, sequence_type=sequence_type)

    bash_string = bash_string_0 + bash_string_1_R

    save_directory = output_directory + "/" + str(sample_name)

    if not os.path.isdir(save_directory):
        os.mkdir(save_directory)

    with open(save_directory + '/job_submission.lsf', 'w') as rsh:
        rsh.write(bash_string)

    os.chdir(save_directory)
    os.system('dos2unix ./job_submission.lsf')
    system_output = os.popen('bsub < ./job_submission.lsf').read()
    # system_output = "Job <1040418> is submitted to queue <large>."
    p1 = re.compile(r'[<](.*?)[>]', re.S)
    JobID = re.findall(p1, system_output)[0]
    # df.loc[i:i, 'JobID'] = JobID
    os.chdir(output_directory)

    return JobID








def preprocessing(gene_type_list, only_filter_num, directory, DE_type, filter_requirment_table_path):


    ### option relative
    from_PreProcessing = False
    from_PreProcessing_submission = False
    from_FPKM_ReadCount_generation = False
    from_gene_info_generation = True
    from_DE = True
    from_info_table_generation = True
    from_filter = False

    if only_filter_num == 0:
        only_filter = False
    elif only_filter_num == 1:
        only_filter = True
    else:
        print("only filter received a wrong input")


    if "microarray" in gene_type_list:
        from_PreProcessing = False
        from_PreProcessing_submission = False


    if only_filter:
        from_PreProcessing = False
        from_PreProcessing_submission = False
        from_FPKM_ReadCount_generation = False
        from_gene_info_generation = False
        from_DE = False
        from_info_table_generation = False
        from_filter = True

    # gene_type_list = ["protein_coding", "lincRNA"] # ["protein_coding", "lincRNA", "microarray"]

    ## path relative
    # directory = "/home/peng.zhou-umw/project/Biomarker_Detection/Biomarker_Detection_Paper/NewCluster_test/GSE54456_wilcox_New_FPKM"

    output_directory = directory + "/output"
    src_path = directory + "/src"

    ## filter condition
    # DE_type = "wilcox" ## "wilcox", "DEseq2", "limma"

    # # pvalue_type = "padj" ##"pvalue"
    # pvalue_type_dict = {}
    # pvalue_type_dict["protein_coding"] = "padj"
    # pvalue_type_dict["lincRNA"] = "padj"
    # pvalue_type_dict["microarray"] = "padj"
    # pvalue_type_dict["Circ"] = "padj"
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
    #
    # # overlap_area_threshold = 0.1
    # overlap_area_threshold_dict = {}
    # overlap_area_threshold_dict["protein_coding"] = 0.1
    # overlap_area_threshold_dict["lincRNA"] = 0.1
    # overlap_area_threshold_dict["microarray"] = 0
    # overlap_area_threshold_dict["Circ"] = 0
    #
    #
    # maxmin_remove_threshold_dict = {}
    # maxmin_remove_threshold_dict["protein_coding"] = 1.0
    # maxmin_remove_threshold_dict["lincRNA"] = 0.001
    # maxmin_remove_threshold_dict["microarray"] = 0
    # maxmin_remove_threshold_dict["Circ"] = 0







    ## doesn't matter for micro array
    # sequence_type = "Single"    # ["Single","Paired"]
    # annotation_file = "LncBook_Version2.0_all"  # ["LncBook_Version2.0_all","gencode_v22", "gencode_v29", "gencode_v37", any path to gtf]
    # inputCSV = directory + "/NewCluster_GSE54456_revised_input_info.csv"

    ##microarray related
    # microarray_file = "/home/peng.zhou-umw/project/Biomarker_Detection/Biomarker_Detection_Paper/Compare_Datasets_2023/Gentles2015/Gentles2015_imputed_Xstandardized_labeled_upsample_standard.csv"
    # microarray_logtransform = str(0) ## if data has been log transformed: 0, if not or not sure: 1
    # microarray_normalization = 1 ## if data has been normalized, or not need to normalize any more: 0, if not or not sure: 1




    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    ML_path = output_directory + "/ML"
    if not os.path.isdir(ML_path):
        os.mkdir(ML_path)

    summary_step1_path = output_directory + "/summary_step1.txt"

    # os.chdir(directory)
    #
    # os.system('chmod +x ./sambamba_view.sh')
    #
    # os.system('chmod +x ./sambamba_sort.sh')
    #
    # os.system('chmod +x ./sed_1.sh')
    #
    # os.system('chmod +x ./sed_2.sh')
    #
    # os.system('chmod +x ./sed_3.sh')

    os.chdir(output_directory)

    # if annotation_file == "gencode_v22":  ## ["gencode_v22", "gencode_v29", "gencode_v37", any path to gtf]
    #     Known_splice_site = '/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/gencode_v22/gencode.v22.annotation.splice_site.txt'
    #     GTF_file = '/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/gencode_v22/gencode.v22.annotation.gtf'
    #     ExonLength_file = '/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/gencode_v22/gencode.v22.annotation.ExonLength'
    # elif annotation_file == "gencode_v29":
    #     Known_splice_site = '/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/gencode_v29/gencode.v29.annotation.splice_site.txt'
    #     GTF_file = '/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/gencode_v29/gencode.v29.annotation.gtf'
    #     ExonLength_file = '/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/gencode_v29/gencode.v29.annotation.ExonLength'
    # elif annotation_file == "gencode_v37":
    #     Known_splice_site = '/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/gencode_v37/gencode.v37.annotation.splice_site.txt'
    #     GTF_file = '/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/gencode_v37/gencode.v37.annotation.gtf'
    #     ExonLength_file = '/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/gencode_v37/gencode.v37.annotation.ExonLength'
    # elif annotation_file == "LncBook_Version2.0_all":
    #     Known_splice_site = '/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/LncBook_Version2.0_all/LncBook_Version2.0_all.annotation.splice_site.txt'
    #     GTF_file = '/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/LncBook_Version2.0_all/LncBook_Version2.0_all.gtf'
    #     ExonLength_file = '/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/LncBook_Version2.0_all/LncBook_Version2.0_all.annotation.ExonLength'
    # else:
    #     create_Known_splice_site_file = subprocess.Popen(["/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/extract_splice_sites.py", annotation_file],stdout=subprocess.PIPE, )
    #     stdout_value = create_Known_splice_site_file.communicate()[0]
    #     tem_folder_path = src_path + "/annotation_file"
    #     if not os.path.isdir(tem_folder_path):
    #         os.mkdir(tem_folder_path)
    #     create_Known_splice_site_file_output = tem_folder_path + "/customized_annotation.splice_site.txt"
    #     logfile = open(create_Known_splice_site_file_output, 'w')
    #     logfile.write(stdout_value.decode('utf-8'))
    #     logfile.close()
    #     Known_splice_site = create_Known_splice_site_file_output
    #     GTF_file = annotation_file
    #     create_ExonLength_file_output = tem_folder_path + "/customized_annotation.ExonLength"
    #     create_ExonLength_file = subprocess.Popen(["Rscript", "/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/ExonLength_forGene_fromgtf.R", annotation_file, create_ExonLength_file_output])
    #     create_ExonLength_file.communicate()
    #     ExonLength_file = create_ExonLength_file_output
    #
    #
    #
    #
    #
    #
    #
    #
    # if from_PreProcessing:
    #
    #     input_df = pd.read_csv(inputCSV, sep=",")
    #
    #     sample_name_list = input_df["sample_name"].values.tolist()
    #     Label_list = input_df["Label"].values.tolist()
    #
    #     Strandness_list = input_df["Strandness"].values.tolist()
    #
    #     if sequence_type == "Paired":
    #         R1_input_list = input_df["R1_input"].values.tolist()
    #         R2_input_list = input_df["R2_input"].values.tolist()
    #     elif sequence_type == "Single":
    #         unpaired_input_list = input_df["unpaired_input"].values.tolist()
    #
    #     Jobid_list = []
    #
    #     if from_PreProcessing_submission:
    #         if sequence_type == "Paired":
    #             for i in range(len(sample_name_list)):
    #                 sample_name = sample_name_list[i]
    #                 R1_input = R1_input_list[i]
    #                 R2_input = R2_input_list[i]
    #                 strandness = Strandness_list[i]
    #                 JobID = job_submission(sample_name=sample_name, R1_input=R1_input, R2_input=R2_input, strandness=strandness,
    #                                        directory=directory, Known_splice_site=Known_splice_site, GTF_file=GTF_file,
    #                                        ExonLength_file=ExonLength_file, sequence_type=sequence_type)
    #                 Jobid_list.append(JobID)
    #         elif sequence_type == "Single":
    #             for i in range(len(sample_name_list)):
    #                 sample_name = sample_name_list[i]
    #                 unpaired_input = unpaired_input_list[i]
    #                 strandness = Strandness_list[i]
    #                 JobID = job_submission(sample_name=sample_name, unpaired_input=unpaired_input, strandness=strandness,
    #                                        directory=directory, Known_splice_site=Known_splice_site, GTF_file=GTF_file,
    #                                        ExonLength_file=ExonLength_file, sequence_type=sequence_type)
    #                 Jobid_list.append(JobID)
    #
    #
    #
    #
    #         JobID_df = pd.DataFrame({"sample" : sample_name_list, "JobID" : Jobid_list})
    #
    #         JobID_path = output_directory + "/JobID_information.txt"
    #
    #         JobID_df.to_csv(JobID_path, sep='\t', index=False)
    #
    #     if not from_PreProcessing_submission:
    #         JobID_path = output_directory + "/JobID_information.txt"
    #         JobID_df = pd.read_csv(JobID_path, sep='\t', index_col=None)
    #
    #     for sample_name in sample_name_list:
    #         file_path = output_directory + "/" + str(sample_name) + "/done.txt"
    #         print(str(sample_name) + " is checking")
    #         while not os.path.exists(file_path):
    #             time.sleep(30)
    #
    #     HTSEQ_count_path_list = []
    #     Gene_FPKM_path_list = []
    #
    #     err_number = 0
    #     for sample_name in sample_name_list:
    #         HTSEQ_count_path = output_directory + "/" + str(sample_name) + "/dump/" + str(sample_name) + "/HTSEQ_count.txt"
    #         Gene_FPKM_path = output_directory + "/" + str(sample_name) + "/dump/" + str(sample_name) + "/Gene_FPKM.txt"
    #         if not os.path.isfile(HTSEQ_count_path):
    #             print(HTSEQ_count_path + " not exist")
    #             err_number = err_number + 1
    #         if not os.path.isfile(Gene_FPKM_path):
    #             print(Gene_FPKM_path + " not exist")
    #             err_number = err_number + 1
    #         HTSEQ_count_path_list.append(HTSEQ_count_path)
    #         Gene_FPKM_path_list.append(Gene_FPKM_path)
    #
    #     if err_number == 0:
    #         print("all files are found")
    #
    #
    #     summary_step1_df = pd.DataFrame({"sample" : sample_name_list, "HTSEQ_count_path" : HTSEQ_count_path_list,
    #                                      "Gene_FPKM" : Gene_FPKM_path_list, "Label" : Label_list})
    #
    #     summary_step1_path = output_directory + "/summary_step1.txt"
    #
    #     summary_step1_df.to_csv(summary_step1_path, sep='\t', index=False, header=False)
    #
    #     # for sample_name in sample_name_list:
    #     #     BAM_path_all = output_directory + "/" + str(sample_name) + "/dump/" + str(sample_name) + "/*bam*"
    #     #     tem_str = "rm -rf " + BAM_path_all
    #     #     os.system(tem_str)
    #
    #     print("step_1 done")
    #
    # os.chdir(output_directory)
    #
    # if from_FPKM_ReadCount_generation:
    #     # step2_output_foldchange_threshold = foldchange_threshold
    #
    #     # step2_output_pvalue_threshold = pvalue_threshold
    #
    #     # if pvalue_type == "padj":
    #     #     step2_output_pvalue_column = str(7)
    #     # elif pvalue_type == "pvalue":
    #     #     step2_output_pvalue_column = str(6)
    #
    #
    #     if annotation_file == "gencode_v22":  ## ["gencode_v22", "gencode_v29", "gencode_v37", any path to gtf]
    #         Pipeline_step_2_1_para_3 = "/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/gencode_v22/gencode_v22_ENSG_Type.txt"
    #     elif annotation_file == "gencode_v29":
    #         Pipeline_step_2_1_para_3 = "/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/gencode_v29/gencode_v29_ENSG_Type.txt"
    #     elif annotation_file == "gencode_v37":
    #         Pipeline_step_2_1_para_3 = "/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/gencode_v37/gencode_v37_ENSG_Type.txt"
    #     elif annotation_file == "LncBook_Version2.0_all":
    #         Pipeline_step_2_1_para_3 = "/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/LncBook_Version2.0_all/LncBook_Version2.0_all_ENSG_Type.txt"
    #     else:
    #         Pipeline_step_2_1_1_tem = subprocess.Popen(['python', "/home/peng.zhou-umw/project/Biomarker_Detection/annotation_file/ENSG_Type.py", annotation_file], stdout=subprocess.PIPE,)
    #         stdout_value = Pipeline_step_2_1_1_tem.communicate()[0]
    #         tem_folder_path = src_path + "/annotation_file"
    #         if not os.path.isdir(tem_folder_path):
    #             os.mkdir(tem_folder_path)
    #         Pipeline_step_2_1_1_tem_output = tem_folder_path + "/customized_ENSG_Type.txt"
    #         logfile = open(Pipeline_step_2_1_1_tem_output, 'w')
    #         logfile.write(stdout_value.decode('utf-8'))
    #         logfile.close()
    #         Pipeline_step_2_1_para_3 = Pipeline_step_2_1_1_tem_output
    #
    #
    #     if "microarray" not in gene_type_list:
    #         Pipeline_step_2_1_para_1 = src_path + '/Expression_Matrix_RNAType_updated.py'
    #         Pipeline_step_2_1_para_2 = summary_step1_path
    #
    #         run_Pipeline_step_2_1 = subprocess.Popen(['python', Pipeline_step_2_1_para_1, Pipeline_step_2_1_para_2, Pipeline_step_2_1_para_3])
    #         run_Pipeline_step_2_1.communicate()
    #
    #     for gene_type in gene_type_list:
    #         if gene_type == "protein_coding":
    #             step2_output_path = summary_step1_path + ".FPKM.ProteinCoding.DEAoutput"
    #             step2_output_reformat_path = output_directory + "/ProteinCoding_output.txt"
    #         elif gene_type == "lincRNA":
    #             step2_output_path = summary_step1_path + ".FPKM.lincRNA.DEAoutput"
    #             step2_output_reformat_path = output_directory + "/lincRNA_output.txt"
    #
    #         def step_2_2(summary_step1_path, src_path, gene_type):
    #
    #             if gene_type == "protein_coding":
    #                 gene_type_para_a = """'NR==1 || $1=="protein_coding" {print}'"""
    #                 gene_type_para_b = "ProteinCoding"
    #             elif gene_type == "lincRNA":
    #                 gene_type_para_a = """'NR==1 || $1=="lincRNA" {print}'"""
    #                 gene_type_para_b = "lincRNA"
    #
    #             run_Pipeline_step_2_2_1_para = summary_step1_path + ".ReadCount"
    #
    #             tem_str = 'awk '+ gene_type_para_a + ' FS="\t" OFS="\t" ' + run_Pipeline_step_2_2_1_para
    #             run_Pipeline_step_2_2_1 = subprocess.Popen(tem_str, shell=True, stdout=subprocess.PIPE, )
    #
    #             run_Pipeline_step_2_2_2 = subprocess.Popen(['cut', "-f", "2-"], stdin=run_Pipeline_step_2_2_1.stdout, stdout=subprocess.PIPE,  )
    #             stdout_value = run_Pipeline_step_2_2_2.communicate()[0]
    #             run_Pipeline_step_2_2_2_output = summary_step1_path + ".ReadCount."+ gene_type_para_b
    #             logfile = open(run_Pipeline_step_2_2_2_output, 'w')
    #             logfile.write(stdout_value.decode('utf-8'))
    #             logfile.close()
    #
    #             run_Pipeline_step_2_2_3_para = summary_step1_path + ".FPKM"
    #
    #             tem_str = 'awk ' + gene_type_para_a + ' FS="\t" OFS="\t" ' + run_Pipeline_step_2_2_3_para
    #             run_Pipeline_step_2_2_3 = subprocess.Popen(tem_str, shell=True, stdout=subprocess.PIPE, )
    #
    #             run_Pipeline_step_2_2_4 = subprocess.Popen(['cut', "-f", "2-"], stdin=run_Pipeline_step_2_2_3.stdout, stdout=subprocess.PIPE,  )
    #             stdout_value = run_Pipeline_step_2_2_4.communicate()[0]
    #             run_Pipeline_step_2_2_4_output = summary_step1_path + ".FPKM."+ gene_type_para_b
    #             logfile = open(run_Pipeline_step_2_2_4_output, 'w')
    #             logfile.write(stdout_value.decode('utf-8'))
    #             logfile.close()
    #
    #             run_Pipeline_step_2_2_5_para_1 = src_path + "/sparse_check.py"
    #             run_Pipeline_step_2_2_5_para_2 = summary_step1_path + ".ReadCount." + gene_type_para_b
    #             run_Pipeline_step_2_2_5 = subprocess.Popen(['python', run_Pipeline_step_2_2_5_para_1, run_Pipeline_step_2_2_5_para_2])
    #             run_Pipeline_step_2_2_5.communicate()
    #
    #
    #         if gene_type != "microarray":
    #             # '($3 > 1 || $3 < -1) && $6<0.05 {print $1}' Popen awk do not need " ' "
    #             # awk_condition = "'($3 > " + step2_output_foldchange_threshold + ' || $3 < -' + step2_output_foldchange_threshold + ') && $' + step2_output_pvalue_column + '<' + step2_output_pvalue_threshold + " {print $1}'"
    #             # print(awk_condition)
    #             step_2_2(summary_step1_path, src_path, gene_type)
    #             print(gene_type + " step_2 done")
    #
    #         elif gene_type == "microarray":
    #             if microarray_normalization == 1:
    #                 run_Pipeline_step_2_2_6_para_1 = src_path + "/limma_Preprocessing.R"
    #             elif microarray_normalization == 0:
    #                 run_Pipeline_step_2_2_6_para_1 = src_path + "/limma_Preprocessing_xNormalization.R"
    #             run_Pipeline_step_2_2_6_para_2 = microarray_file
    #             run_Pipeline_step_2_2_6_para_3 = microarray_logtransform
    #             run_Pipeline_step_2_2_6_para_4 = output_directory + "/summary_step1.txt.microarray"
    #             #run_Pipeline_step_2_2_6_para_5 = output_directory + "/summary_step1.txt.microarray.label"
    #             run_Pipeline_step_2_2_6 = subprocess.Popen(['Rscript', run_Pipeline_step_2_2_6_para_1,
    #                                                         run_Pipeline_step_2_2_6_para_2, run_Pipeline_step_2_2_6_para_3,
    #                                                         run_Pipeline_step_2_2_6_para_4])
    #             run_Pipeline_step_2_2_6.communicate()

    ## filter information tables generation
    if from_gene_info_generation:
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

            if gene_type in ["protein_coding", "lincRNA"]:
                run_Pipeline_step_2_2_1_para = output_directory + "/summary_step1.txt.ReadCount"

                gene_type_para_a = """'NR==1 || $1=="protein_coding" || $1=="lincRNA" {print}'"""

                tem_str = 'awk '+ gene_type_para_a + ' FS="\t" OFS="\t" ' + run_Pipeline_step_2_2_1_para
                run_Pipeline_step_2_2_1 = subprocess.Popen(tem_str, shell=True, stdout=subprocess.PIPE, )

                run_Pipeline_step_2_2_2 = subprocess.Popen(['cut', "-f", "1,2"], stdin=run_Pipeline_step_2_2_1.stdout, stdout=subprocess.PIPE,  )
                stdout_value = run_Pipeline_step_2_2_2.communicate()[0]
                run_Pipeline_step_2_2_2_output = output_directory + "/summary_step1.txt.ReadCount.GeneTypeInfo"
                logfile = open(run_Pipeline_step_2_2_2_output, 'w')
                logfile.write(stdout_value.decode('utf-8'))
                logfile.close()

                run_Pipeline_step_2_2_1_para = output_directory + "/summary_step1.txt.ReadCount"

                gene_type_para_a = """'NR==1 || $1=="protein_coding" || $1=="lincRNA" {print}'"""

                tem_str = 'awk '+ gene_type_para_a + ' FS="\t" OFS="\t" ' + run_Pipeline_step_2_2_1_para
                run_Pipeline_step_2_2_1 = subprocess.Popen(tem_str, shell=True, stdout=subprocess.PIPE, )

                run_Pipeline_step_2_2_2 = subprocess.Popen(['cut', "-f", "2-"], stdin=run_Pipeline_step_2_2_1.stdout, stdout=subprocess.PIPE,  )
                stdout_value = run_Pipeline_step_2_2_2.communicate()[0]
                run_Pipeline_step_2_2_2_output = output_directory + "/summary_step1.txt.ReadCount.all"
                logfile = open(run_Pipeline_step_2_2_2_output, 'w')
                logfile.write(stdout_value.decode('utf-8'))
                logfile.close()


            FPKM_path = output_directory + "/" + FPKM_name

            ReadCount_path = output_directory + "/" + ReadCount_name

            if from_DE:
                if gene_type == "microarray":
                    microarray_two_power_output_path = ReadCount_path + ".two_power"
                    tem_microarray_df = pd.read_csv(ReadCount_path, sep="\t", index_col=0).T
                    tem_column_name_list = tem_microarray_df.columns.values.tolist()
                    for column_name in tem_column_name_list:
                        tem_microarray_df[column_name] = 2 ** tem_microarray_df[column_name]
                    tem_microarray_two_power_df = tem_microarray_df.T
                    tem_microarray_two_power_df.to_csv(microarray_two_power_output_path, sep="\t", index=True)
                else:
                    microarray_two_power_output_path = "./"


                if DE_type == "DEseq2":
                    run_Pipeline_step_2_2_6_para_1 = src_path + "/DEseq2.R"
                    if gene_type == "microarray":
                        run_Pipeline_step_2_2_6_para_2 = microarray_two_power_output_path
                    else:
                        run_Pipeline_step_2_2_6_para_2 = ReadCount_path
                    run_Pipeline_step_2_2_6 = subprocess.Popen(
                        ['Rscript', run_Pipeline_step_2_2_6_para_1, run_Pipeline_step_2_2_6_para_2])
                    run_Pipeline_step_2_2_6.communicate()
                    # tem_target_path = output_directory + "/" + ReadCount_name + ".DE_tem.output"
                    # tem_str = "mv " + run_Pipeline_step_2_2_6_para_2 + ".DESeq2.output" + " " + tem_target_path
                    # os.system(tem_str)

                elif DE_type == "wilcox":
                    run_Pipeline_step_2_2_6_para_1 = src_path + "/wilcox_test/wilcox_pipeline.py"
                    if gene_type in ["protein_coding", "lincRNA"]:
                        run_Pipeline_step_2_2_6_para_2 = output_directory + "/summary_step1.txt.ReadCount.all"
                    else:
                        run_Pipeline_step_2_2_6_para_2 = ReadCount_path
                    run_Pipeline_step_2_2_6_para_3 = src_path + "/wilcox_test/wilcox_test.R"
                    run_Pipeline_step_2_2_6_para_4 = src_path + "/wilcox_test/wilcox_test_microarray.R"
                    run_Pipeline_step_2_2_6_para_5 = microarray_two_power_output_path

                    if gene_type == "microarray":
                        run_Pipeline_step_2_2_6_para_6 = str(1)
                    else:
                        run_Pipeline_step_2_2_6_para_6 = str(0)


                    run_Pipeline_step_2_2_6 = subprocess.Popen(
                        ['python', run_Pipeline_step_2_2_6_para_1, run_Pipeline_step_2_2_6_para_2,
                         run_Pipeline_step_2_2_6_para_3, run_Pipeline_step_2_2_6_para_4,
                         run_Pipeline_step_2_2_6_para_5, run_Pipeline_step_2_2_6_para_6])
                    run_Pipeline_step_2_2_6.communicate()
                    # tem_ori_path = ReadCount_path + ".DE_tem.output"
                    # tem_target_path = output_directory + "/" + ReadCount_name + ".DE_tem.output"
                    # tem_str = "mv " + tem_ori_path + " " + tem_target_path
                    # os.system(tem_str)
                    #
                    # tem_ori_path = ReadCount_path + ".condition.csv"
                    # tem_target_path = output_directory + "/" + ReadCount_name + ".condition.csv"
                    # tem_str = "mv " + tem_ori_path + " " + tem_target_path
                    # os.system(tem_str)
                elif DE_type == "limma":
                    # microarray limma
                    run_Pipeline_step_2_2_6_para_1 = src_path + "/limmaDE/limmaDE_pipeline.py"
                    run_Pipeline_step_2_2_6_para_2 = ReadCount_path
                    run_Pipeline_step_2_2_6_para_3 = src_path + "/limmaDE/limmaDE.R"
                    run_Pipeline_step_2_2_6 = subprocess.Popen(
                        ['python', run_Pipeline_step_2_2_6_para_1, run_Pipeline_step_2_2_6_para_2,
                         run_Pipeline_step_2_2_6_para_3])
                    run_Pipeline_step_2_2_6.communicate()
                    # tem_ori_path = ReadCount_path + ".DE_tem.output"
                    # tem_target_path = DE_output_folder_path + "/" + ReadCount_name + ".DE_tem.output"
                    # tem_str = "mv " + tem_ori_path + " " + tem_target_path
                    # os.system(tem_str)
                    #
                    # tem_ori_path = ReadCount_path + ".condition.csv"
                    # tem_target_path = DE_output_folder_path + "/" + ReadCount_name + ".condition.csv"
                    # tem_str = "mv " + tem_ori_path + " " + tem_target_path
                    # os.system(tem_str)

                if ((gene_type in ["protein_coding", "lincRNA"]) and (DE_type == "wilcox")):
                    # os.system("cp " + output_directory + "/summary_step1.txt.ReadCount.all" + ".DE_tem.output " + output_directory + "/" + ReadCount_name + ".DE_tem.output")
                    tem_overall_df = pd.read_csv(output_directory + "/summary_step1.txt.ReadCount.all.DE_tem.output", sep="\t", index_col=0)
                    gene_type_info_df = pd.read_csv(output_directory + "/summary_step1.txt.ReadCount.GeneTypeInfo", sep="\t")

                    # required_gene_type = protein_coding # "lincRNA"

                    required_gene_list = gene_type_info_df[gene_type_info_df["GeneType"] == gene_type].copy()["Gene"].tolist()

                    tem_list = []
                    for tem_gene_name in tem_overall_df.index.tolist():
                        if tem_gene_name in required_gene_list:
                            tem_list.append(1)
                        else:
                            tem_list.append(0)

                    tem_overall_df["flag"] = tem_list

                    tem_output_df = tem_overall_df[tem_overall_df["flag"] == 1].copy()

                    tem_output_df.drop(["flag"], axis=1, inplace=True)

                    tem_output_df.to_csv(output_directory + "/" + ReadCount_name + ".DE_tem.output", sep="\t")



                run_Pipeline_step_2_2_6_5_para_0 = src_path + "/FPKM_log2foldchange_calculation.py"
                run_Pipeline_step_2_2_6_5_para_1 = output_directory + "/" + ReadCount_name + ".DE_tem.output"
                run_Pipeline_step_2_2_6_5_para_2 = FPKM_path
                run_Pipeline_step_2_2_6_5_para_3 = output_directory + "/" + ReadCount_name + ".DE.output"
                if gene_type == "microarray":
                    run_Pipeline_step_2_2_6_5_para_4 = str(1)
                else:
                    run_Pipeline_step_2_2_6_5_para_4 = str(0)

                run_Pipeline_step_2_2_6_5 = subprocess.Popen(
                    ['python', run_Pipeline_step_2_2_6_5_para_0, run_Pipeline_step_2_2_6_5_para_1,
                     run_Pipeline_step_2_2_6_5_para_2, run_Pipeline_step_2_2_6_5_para_3, run_Pipeline_step_2_2_6_5_para_4])
                run_Pipeline_step_2_2_6_5.communicate()

            if from_info_table_generation:

                reformat_ori = FPKM_path

                reformat_to = output_directory + "/" + output_name + "_output.all.txt"

                if os.path.exists(reformat_to):
                    os.remove(reformat_to)


                reformat_1 = subprocess.Popen(['cp', '-i', reformat_ori, reformat_to])

                reformat_1.communicate()

                df_ori_input = pd.read_csv(reformat_to, sep="\t", index_col=0)

                df_ori_input = df_ori_input.T

                df_index_name = df_ori_input.index.tolist()

                label_list = []

                for i in df_index_name:
                    if i.split("_")[-1] == "C":
                        label_list.append(int(0))
                    elif i.split("_")[-1] == "T":
                        label_list.append(int(1))

                y_dict = {"index": df_index_name,
                          "y": label_list}

                df_y = pd.DataFrame(y_dict)

                df_y.set_index(["index"], inplace=True)

                revised_df = pd.concat([df_ori_input, df_y], axis=1)

                reformat_output_path = output_directory + '/reformat_' + gene_type + '_input.all.txt'

                if os.path.exists(reformat_output_path):
                    os.remove(reformat_output_path)


                revised_df.to_csv(reformat_output_path, header=True, sep='\t')

                #calculate the bioinfo information

                one_record_list = []
                overlap_record = []

                X_column_list = df_ori_input.columns.tolist()

                padj_df = pd.read_csv(output_directory + "/" + ReadCount_name + ".DE.output", sep="\t", index_col=0)[
                    ["padj", "pvalue", "log2FoldChange", "log2FoldChange_mean"]]
                FPKM_df = pd.read_csv(FPKM_path, sep="\t", index_col=0)
                if gene_type == "microarray":
                    tem_FPKM_df_transposition = FPKM_df.T
                    column_name_list = tem_FPKM_df_transposition.columns.values.tolist()
                    for column_name in column_name_list:
                        tem_FPKM_df_transposition[column_name] = 2 ** tem_FPKM_df_transposition[column_name]
                    FPKM_df = tem_FPKM_df_transposition.T
                FPKM_df = FPKM_df.loc[~(FPKM_df == 0).all(axis=1)]
                FPKM_df_transposition = FPKM_df.T#[X_column_list]

                Label_list = []
                for i in range(len(FPKM_df_transposition)):
                    if str(FPKM_df_transposition.index.values[i]).endswith('_C'):
                        Label_list.append(0)
                    elif str(FPKM_df_transposition.index.values[i]).endswith('_T'):
                        Label_list.append(1)
                    else:
                        print("Label not match")

                pseudo_foldchange = float(1e-6)
                FPKM_df_transposition['label'] = Label_list
                FPKM_df_transposition_0 = FPKM_df_transposition[FPKM_df_transposition['label'].isin([0])]
                FPKM_df_transposition_1 = FPKM_df_transposition[FPKM_df_transposition['label'].isin([1])]

                FPKM_df_transposition_0_transform = FPKM_df_transposition_0.drop(['label'], axis=1).T
                FPKM_df_transposition_1_transform = FPKM_df_transposition_1.drop(['label'], axis=1).T

                FPKM_df_0_copy = FPKM_df_transposition_0_transform.copy()
                FPKM_df_1_copy = FPKM_df_transposition_1_transform.copy()

                FPKM_mean_C_series = FPKM_df_transposition_0_transform.mean(axis=1)
                FPKM_mean_T_series = FPKM_df_transposition_1_transform.mean(axis=1)

                FPKM_std_C_series = FPKM_df_transposition_0_transform.std(axis=1)
                FPKM_std_T_series = FPKM_df_transposition_1_transform.std(axis=1)


                # FPKM_df_transposition_0_transform["mean"] = FPKM_df_transposition_0_transform.mean(axis=1)
                # FPKM_df_transposition_1_transform["mean"] = FPKM_df_transposition_1_transform.mean(axis=1)
                #
                # FPKM_df_transposition_0_transform["std"] = FPKM_df_transposition_0_transform.std(axis=1)
                # FPKM_df_transposition_1_transform["std"] = FPKM_df_transposition_1_transform.std(axis=1)

                FPKM_df_transposition_0_transform["mean"] = FPKM_mean_C_series
                FPKM_df_transposition_1_transform["mean"] = FPKM_mean_T_series

                FPKM_df_transposition_0_transform["std"] = FPKM_std_C_series
                FPKM_df_transposition_1_transform["std"] = FPKM_std_T_series


                FPKM_df_C_mean = FPKM_df_transposition_0_transform[["mean"]].rename(columns={'mean': 'mean_Control'})
                FPKM_df_T_mean = FPKM_df_transposition_1_transform[["mean"]].rename(columns={'mean': 'mean_Treatment'})

                FPKM_df_C_std = FPKM_df_transposition_0_transform[["std"]].rename(columns={'std': 'std_Control'})
                FPKM_df_T_std = FPKM_df_transposition_1_transform[["std"]].rename(columns={'std': 'std_Treatment'})

                FPKM_df_with_mean = pd.merge(FPKM_df_C_mean, FPKM_df_T_mean, left_index=True, right_index=True)
                FPKM_df_with_std = pd.merge(FPKM_df_C_std, FPKM_df_T_std, left_index=True, right_index=True)

                FPKM_df_with_all = pd.merge(FPKM_df_with_mean, FPKM_df_with_std, left_index=True, right_index=True)

                # abs_log2_list = []
                # foldchange_sd_score_list = []
                # for i in range(len(FPKM_df_with_all)):
                #     mean_C = FPKM_df_with_all.iloc[i].at['mean_Control']
                #     mean_T = FPKM_df_with_all.iloc[i].at['mean_Treatment']
                #     if mean_C == 0 or mean_T == 0:
                #         mean_C = mean_C + pseudo_foldchange
                #         mean_T = mean_T + pseudo_foldchange
                #     tem_abs_log2_foldchange = abs(math.log(mean_T / mean_C, 2))
                #     abs_log2_list.append(tem_abs_log2_foldchange)
                #
                #     std_C = FPKM_df_with_all.iloc[i].at['std_Control']
                #     std_T = FPKM_df_with_all.iloc[i].at['std_Treatment']
                #     max_std = max(std_C, std_T)
                #     if max_std == 0:
                #         max_std = max_std + pseudo_foldchange
                #         tem_abs_log2_foldchange = tem_abs_log2_foldchange + pseudo_foldchange
                #
                #     tem_foldchange_sd_score = tem_abs_log2_foldchange / max_std
                #     foldchange_sd_score_list.append(tem_foldchange_sd_score)

                # FPKM_df_with_all["abs(log2(fold-change))"] = abs_log2_list
                # FPKM_df_with_all["fold-change_std_score"] = foldchange_sd_score_list

                abs_log2_list = []
                for tem_foldchange in padj_df["log2FoldChange"].tolist():
                    abs_log2_list.append(abs(tem_foldchange))
                padj_df["abs(log2(fold-change))"] = abs_log2_list

                abs_log2_mean_list = []
                for tem_foldchange in padj_df["log2FoldChange_mean"].tolist():
                    abs_log2_mean_list.append(abs(tem_foldchange))
                padj_df["abs(log2(fold-change_mean))"] = abs_log2_mean_list



                info_df = pd.merge(FPKM_df_with_all, padj_df, left_index=True, right_index=True, how='left')

                dict_FPKM_0 = FPKM_df_0_copy.T.to_dict(orient='list')
                dict_FPKM_1 = FPKM_df_1_copy.T.to_dict(orient='list')

                tem_overlap_dict = {}
                tem_overlap_dict["tem_index"] = []
                tem_overlap_dict["overlap_area"] = []
                random.seed(7)

                for gene_name in dict_FPKM_0:
                    tem_overlap_dict["tem_index"].append(gene_name)
                    # print(gene_name)
                    tem_C_list = dict_FPKM_0[gene_name]
                    tem_T_list = dict_FPKM_1[gene_name]

                    ## remove outlier
                    # Q1_C = np.percentile(tem_C_list, 25, interpolation='midpoint')
                    #
                    # Q3_C = np.percentile(tem_C_list, 75, interpolation='midpoint')
                    # IQR_C = Q3_C - Q1_C
                    #
                    # tem_C_list = [item for item in tem_C_list if item <= (Q3_C+(1.5*IQR_C)) and item >= (Q1_C-(1.5*IQR_C))]
                    #
                    # Q1_T = np.percentile(tem_T_list, 25, interpolation='midpoint')
                    #
                    # Q3_T = np.percentile(tem_T_list, 75, interpolation='midpoint')
                    # IQR_T = Q3_T - Q1_T
                    #
                    # tem_T_list = [item for item in tem_T_list if item <= (Q3_T + (1.5 * IQR_T)) and item >= (Q1_T - (1.5 * IQR_T))]



                    # only keep the range between distributions
                    xmax_min = max(min(tem_C_list), min(tem_T_list))
                    xmin_max = min(max(tem_C_list), max(tem_T_list))

                    # it is possible two distribution totally depart, so it need exchange
                    xmin = min(xmax_min, xmin_max)
                    xmax = max(xmax_min, xmin_max)

                    dx = 0.2 * (xmax - xmin)  # add a 20% margin, as the kde is wider than the data
                    xmin -= dx
                    xmax += dx
                    if xmin < 0:
                        xmin = 0

                    generated_sample_x = np.linspace(xmin, xmax, 1024)
                    diff_record = generated_sample_x[1] - generated_sample_x[0]
                    epsilon_info = sys.float_info.epsilon
                    flag_allZero_C = 0
                    flag_allZero_T = 0
                    try:
                        y_C_kernel = NaiveKDE(kernel='gaussian', bw='silverman').fit(np.array(tem_C_list))
                        y_C_sample = y_C_kernel.evaluate(generated_sample_x)
                    except:
                        flag_allZero_C = 1
                    try:
                        y_T_kernel = NaiveKDE(kernel='gaussian', bw='silverman').fit(np.array(tem_T_list))
                        y_T_sample = y_T_kernel.evaluate(generated_sample_x)
                    except:
                        flag_allZero_T = 1

                    if flag_allZero_C == 1 and flag_allZero_T == 1:
                        area_inters_x = 1.0
                    elif flag_allZero_C == 1 or flag_allZero_T == 1:
                        area_inters_x = 0.0
                    else:
                        inters_x = np.minimum(y_C_sample, y_T_sample)

                        # because the overlap area maybe very narraw which does not include enough sample to calculate the area
                        # so we need to resample for each non-zero area

                        for index_non_epsilon in range(len(inters_x)):
                            if inters_x[index_non_epsilon] != epsilon_info:
                                break

                        flag_index = 0
                        calculation_dict = {}
                        calculation_dict[flag_index] = {}
                        calculation_dict[flag_index]["x"] = []
                        calculation_dict[flag_index]["y"] = []
                        for i in range(index_non_epsilon, len(inters_x)):
                            if inters_x[i] != epsilon_info:
                                calculation_dict[flag_index]["x"].append(generated_sample_x[i])
                                calculation_dict[flag_index]["y"].append(inters_x[i])
                            else:
                                if len(calculation_dict[flag_index]["x"]) == 0:
                                    continue
                                else:
                                    flag_index = flag_index + 1
                                    calculation_dict[flag_index] = {}
                                    calculation_dict[flag_index]["x"] = []
                                    calculation_dict[flag_index]["y"] = []
                        sum_record = 0
                        for key in calculation_dict:
                            if len(calculation_dict[key]["x"]) != 0:
                                tem_min = min(calculation_dict[key]["x"])
                                tem_max = max(calculation_dict[key]["x"])
                                tem_generated_sample_x = np.linspace(tem_min, tem_max, 1024)
                                tem_y_C_sample = y_C_kernel.evaluate(tem_generated_sample_x)
                                tem_y_T_sample = y_T_kernel.evaluate(tem_generated_sample_x)
                                tem_inters_x = np.minimum(tem_y_C_sample, tem_y_T_sample)
                                tem_score = np.trapz(tem_inters_x, tem_generated_sample_x)

                                ## it is possible non-zero area only include few sample which not enough for area calculation
                                if tem_score > 1:
                                    tem_length_one_list = calculation_dict[key]["x"].copy()
                                    tem_length_one_extension_min = min(tem_length_one_list) #- diff_record
                                    if tem_length_one_extension_min < 0:
                                        tem_length_one_extension_min = 0
                                    tem_length_one_extension_max = max(tem_length_one_list) #+ diff_record
                                    tem_length_one_extension_generated_sample_x = np.linspace(tem_length_one_extension_min,
                                                                                              tem_length_one_extension_max, 1024 * 1024)
                                    tem_length_one_extension_y_C_sample = y_C_kernel.evaluate(
                                        tem_length_one_extension_generated_sample_x)
                                    tem_length_one_extension_y_T_sample = y_T_kernel.evaluate(
                                        tem_length_one_extension_generated_sample_x)
                                    tem_length_one_extension_inters_x = np.minimum(tem_length_one_extension_y_C_sample,
                                                                                   tem_length_one_extension_y_T_sample)
                                    sum_record = sum_record + np.trapz(tem_length_one_extension_inters_x,
                                                                       tem_length_one_extension_generated_sample_x)
                                else:
                                    sum_record = sum_record + tem_score
                        area_inters_x = sum_record
                    if area_inters_x > 1:
                        area_inters_x = 1.0
                    if area_inters_x == 1.0:
                        one_record_list.append(gene_name)
                        overlap_record.append(area_inters_x)
                    tem_overlap_dict["overlap_area"].append(area_inters_x)

                overlap_info_df = pd.DataFrame(tem_overlap_dict)
                overlap_info_df.set_index("tem_index", inplace=True)
                info_df = pd.merge(info_df, overlap_info_df, left_index=True, right_index=True, how='left')

                #mean difference
                difference_list = []
                C_tem_list = info_df["mean_Control"].tolist()
                T_tem_list = info_df["mean_Treatment"].tolist()

                for i in range(len(info_df["mean_Control"].tolist())):
                    difference_list.append(abs(C_tem_list[i] - T_tem_list[i]))
                info_df["mean_diff"] = difference_list

                #median difference
                median_record_dict = {}
                median_record_dict["gene_name"] = []
                median_record_dict["C_median"] = []
                median_record_dict["T_median"] = []
                median_record_dict["median_diff"] = []

                for gene_name in dict_FPKM_0:
                    median_record_dict["gene_name"].append(gene_name)
                    tem_C_median = np.median(np.array(dict_FPKM_0[gene_name]))
                    tem_T_median = np.median(np.array(dict_FPKM_1[gene_name]))
                    median_record_dict["C_median"].append(tem_C_median)
                    median_record_dict["T_median"].append(tem_T_median)
                    median_record_dict["median_diff"].append(abs(tem_C_median - tem_T_median))
                median_info_df = pd.DataFrame(median_record_dict)
                median_info_df.set_index("gene_name", inplace=True)
                info_df = pd.merge(info_df, median_info_df, left_index=True, right_index=True, how='left')

                overlap_list = info_df["overlap_area"].tolist()
                median_diff_list = info_df["median_diff"].tolist()

                update_exp_overlap_list = []
                update_score_list = []

                for i in range(len(overlap_list)):
                    overlap = overlap_list[i]
                    exp_overlap = math.exp(overlap * 10)
                    update_exp_overlap_list.append(exp_overlap)
                    update_score_list.append(median_diff_list[i] / exp_overlap)

                info_df["exp_overlap"] = update_exp_overlap_list
                info_df["score"] = update_score_list

                #maximum expression of genes
                max_val_record_dict = {}
                max_val_record_dict["gene_name"] = []
                max_val_record_dict["max_val"] = []


                for gene_name in dict_FPKM_0:
                    max_val_record_dict["gene_name"].append(gene_name)
                    tem_combine_list = dict_FPKM_0[gene_name] + dict_FPKM_1[gene_name]
                    max_val_record_dict["max_val"].append(max(list(map(abs,tem_combine_list))))

                max_val_info_df = pd.DataFrame(max_val_record_dict)
                max_val_info_df.set_index("gene_name", inplace=True)
                info_df = pd.merge(info_df, max_val_info_df, left_index=True, right_index=True, how='left')

                # info_df.to_csv(output_directory + '/' + gene_type+".csv", header=True, sep='\t')

                info_df = info_df[
                    ["pvalue", "padj", "abs(log2(fold-change))", "abs(log2(fold-change_mean))", "max_val", "mean_Control", "mean_Treatment",
                     "mean_diff", "std_Control", "std_Treatment", "C_median", "T_median", "median_diff", "overlap_area",
                     "exp_overlap", "score"]]

                info_table_output_path = output_directory + '/' + gene_type + '_Gene_info.csv'

                if os.path.exists(info_table_output_path):
                    os.remove(info_table_output_path)

                info_df.to_csv(info_table_output_path, header=True, sep='\t')


                info_table_for_review = {"Gene":info_df["Gene"].tolist(),
                                         "pvalue":info_df["pvalue"].tolist(),
                                         "padj":info_df["padj"].tolist(),
                                         "abs(log2(fold-change))":info_df["abs(log2(fold-change))"].tolist(),
                                         "max_val":info_df["max_val"].tolist(),
                                         "distribution_overlap_area":info_df["overlap_area"].tolist(),
                                         }

                current_time = time.strftime('%Y%m%d%H%M', time.localtime(time.time()))
                Copy_info_table_output_path = directory + '/' + gene_type + '_Gene_info_' + current_time + '.xlsx'

                if os.path.exists(Copy_info_table_output_path):
                    os.remove(Copy_info_table_output_path)


                info_table_for_review.to_excel(Copy_info_table_output_path, header=True, index=False)

                record_df = pd.DataFrame({"gene": one_record_list, "overlap": overlap_record})

                if os.path.exists(output_directory + '/' + gene_type + 'overlap_one_info.csv'):
                    os.remove(output_directory + '/' + gene_type + 'overlap_one_info.csv')

                record_df.to_csv(output_directory + '/' + gene_type + 'overlap_one_info.csv', header=True, sep='\t', index=False)

                info_table_output_path




    if from_filter:

        Filter_num_record_dict = {}
        Filter_num_record_dict["step"] = []
        Filter_num_record_dict["Gene_number"] = []


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

            ## filter by DE p-value

            info_df = pd.read_csv(output_directory + '/' + gene_type + '_Gene_info.csv', sep='\t')

            Filter_num_record_dict["step"].append(gene_type + " Original")
            Filter_num_record_dict["Gene_number"].append(len(info_df))

            DE_filtered_info_df = info_df[(info_df[pvalue_type_dict[gene_type]] <= float(pvalue_threshold_dict[gene_type]))]

            Filter_num_record_dict["step"].append(gene_type + " after DE pvalue filter")
            Filter_num_record_dict["Gene_number"].append(len(DE_filtered_info_df))

            ## filter by foldchange

            FC_filtered_info_df = DE_filtered_info_df[(DE_filtered_info_df['abs(log2(fold-change))'] >= float(foldchange_threshold_dict[gene_type]))]

            Filter_num_record_dict["step"].append(gene_type + " after foldchange filter")
            Filter_num_record_dict["Gene_number"].append(len(FC_filtered_info_df))


            ## filter by maximum value

            max_val_filtered_info_df = FC_filtered_info_df[FC_filtered_info_df["max_val"]>=maxmin_remove_threshold_dict[gene_type]]

            Filter_num_record_dict["step"].append(gene_type + " after maximum value filter")
            Filter_num_record_dict["Gene_number"].append(len(max_val_filtered_info_df))

            ## filter by overlap area

            overlap_filtered_info_df = max_val_filtered_info_df[max_val_filtered_info_df["overlap_area"]<=overlap_area_threshold_dict[gene_type]]

            Filter_num_record_dict["step"].append(gene_type + " after overlap area filter")
            Filter_num_record_dict["Gene_number"].append(len(overlap_filtered_info_df))

            if os.path.exists(output_directory + '/Filtered_' + gene_type + '_Gene_info.csv'):
                os.remove(output_directory + '/Filtered_' + gene_type + '_Gene_info.csv')

            overlap_filtered_info_df.to_csv(output_directory + '/Filtered_' + gene_type + '_Gene_info.csv', sep="\t", index=False)


            filtered_Gene_list = overlap_filtered_info_df["Gene"].tolist()

            FPKM_df = pd.read_csv(output_directory + "/" + FPKM_name, sep="\t")

            Filtered_FPKM_df = FPKM_df[FPKM_df['Gene'].isin(filtered_Gene_list)]

            if os.path.exists(output_directory + "/Filtered_" + FPKM_name):
                os.remove(output_directory + "/Filtered_" + FPKM_name)

            Filtered_FPKM_df.to_csv(output_directory + "/Filtered_" + FPKM_name, sep="\t", index=False)


            ReadCount_df = pd.read_csv(output_directory + "/" + ReadCount_name, sep="\t")

            Filtered_ReadCount_df = ReadCount_df[ReadCount_df['Gene'].isin(filtered_Gene_list)]

            if os.path.exists(output_directory + "/Filtered_" + ReadCount_name):
                os.remove(output_directory + "/Filtered_" + ReadCount_name)

            Filtered_ReadCount_df.to_csv(output_directory + "/Filtered_" + ReadCount_name, sep="\t", index=False)



            df_ori_input = pd.read_csv(output_directory + "/Filtered_" + FPKM_name, sep="\t", index_col=0).T

            df_index_name = df_ori_input.index.tolist()

            label_list = []

            for i in df_index_name:
                if i.split("_")[-1] == "C":
                    label_list.append(int(0))
                elif i.split("_")[-1] == "T":
                    label_list.append(int(1))

            # y_dict = {"index": df_index_name,
            #           "y": label_list}
            #
            # df_y = pd.DataFrame(y_dict)
            #
            # df_y.set_index(["index"], inplace=True)

            # revised_df = pd.concat([df_ori_input, df_y], axis=1)

            df_ori_input["y"] = label_list

            reformat_output_path = output_directory + '/Filtered_' + gene_type + '_input.all.txt'

            if os.path.exists(reformat_output_path):
                os.remove(reformat_output_path)

            df_ori_input.to_csv(reformat_output_path, header=True, sep='\t')

        Filter_num_record_output_path = output_directory + '/Filtered_info_record.csv'

        if os.path.exists(Filter_num_record_output_path):
            os.remove(Filter_num_record_output_path)

        pd.DataFrame(Filter_num_record_dict).to_csv(Filter_num_record_output_path, index=False)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--gene_type_list", nargs='+', default=["protein_coding", "lincRNA", "microarray"], required=True)
    parser.add_argument("--only_filter_num", default=None, type=int, required=True)
    parser.add_argument("--directory", default=None, type=str, required=True)
    parser.add_argument("--DE_type", default=None, type=str, required=True)
    parser.add_argument("--filter_requirment_table_path", default=None, type=str, required=True)



    args = parser.parse_args()

    preprocessing(gene_type_list=args.gene_type_list, only_filter_num=args.only_filter_num,
                          directory=args.directory, DE_type=args.DE_type,
                          filter_requirment_table_path=args.filter_requirment_table_path)



    # stab_selection_main(folder_path=args.folder_path, file_name=args.train_file_name, testfile_name=args.test_file_name, threshold=args.threshold)



if __name__ == '__main__':
    if not sys.warnoptions:
        warnings.simplefilter("ignore")
        os.environ["PYTHONWARNINGS"] = "ignore"  # Also affect subprocesses
    main()
