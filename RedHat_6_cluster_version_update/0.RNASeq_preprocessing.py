
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
############################################################

### parameter need to set by user:

# gene_type_list = ["protein_coding", "lincRNA"] # ["protein_coding", "lincRNA"]
# sequence_type = "Paired"    # ["Single","Paired"]
# annotation_file = "LncBook_Version2.0_all"  # ["LncBook_Version2.0_all","gencode_v22", "gencode_v29", "gencode_v37", any path to gtf]
# inputCSV = "./0.RNASeq_preprocessing_input_sample_Paired-End.csv"

############################################################

def RNASeq_preprocessing(biomarker_target_gene_type = "protein_coding", sequence_type="Paired", annotation_file="LncBook_Version2.0_all", inputCSV="./0.RNASeq_preprocessing_input_sample_Paired-End.csv"):
    ### option relative
    if biomarker_target_gene_type == "both":
        gene_type_list = ["protein_coding", "lincRNA"]
    else:
        gene_type_list = [biomarker_target_gene_type]

    from_PreProcessing = True
    from_PreProcessing_submission = True
    from_FPKM_ReadCount_generation = True

    ## path relative
    directory = os.getcwd()
    output_directory = directory + "/output"
    src_path = directory + "/src"

    annotation_file_folder_path = src_path + "/annotation_file"

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    ML_path = output_directory + "/ML"
    if not os.path.isdir(ML_path):
        os.mkdir(ML_path)

    summary_step1_path = output_directory + "/summary_step1.txt"

    os.chdir(directory)

    os.system('chmod +x ./sambamba_view.sh')

    os.system('chmod +x ./sambamba_sort.sh')

    os.system('chmod +x ./sed_1.sh')

    os.system('chmod +x ./sed_2.sh')

    os.system('chmod +x ./sed_3.sh')

    os.chdir(output_directory)

    if annotation_file == "gencode_v22":  ## ["gencode_v22", "gencode_v29", "gencode_v37", any path to gtf]
        Known_splice_site = annotation_file_folder_path + '/gencode_v22/gencode.v22.annotation.splice_site.txt'
        GTF_file = annotation_file_folder_path + '/gencode_v22/gencode.v22.annotation.gtf'
        ExonLength_file = annotation_file_folder_path + '/gencode_v22/gencode.v22.annotation.ExonLength'
    elif annotation_file == "gencode_v29":
        Known_splice_site = annotation_file_folder_path + '/gencode_v29/gencode.v29.annotation.splice_site.txt'
        GTF_file = annotation_file_folder_path + '/gencode_v29/gencode.v29.annotation.gtf'
        ExonLength_file = annotation_file_folder_path + '/gencode_v29/gencode.v29.annotation.ExonLength'
    elif annotation_file == "gencode_v37":
        Known_splice_site = annotation_file_folder_path + '/gencode_v37/gencode.v37.annotation.splice_site.txt'
        GTF_file = annotation_file_folder_path + '/gencode_v37/gencode.v37.annotation.gtf'
        ExonLength_file = annotation_file_folder_path + '/gencode_v37/gencode.v37.annotation.ExonLength'
    elif annotation_file == "LncBook_Version2.0_all":
        Known_splice_site = annotation_file_folder_path + '/LncBook_Version2.0_all/LncBook_Version2.0_all.annotation.splice_site.txt'
        GTF_file = annotation_file_folder_path + '/LncBook_Version2.0_all/LncBook_Version2.0_all.gtf'
        ExonLength_file = annotation_file_folder_path + '/LncBook_Version2.0_all/LncBook_Version2.0_all.annotation.ExonLength'
    else:
        create_Known_splice_site_file = subprocess.Popen(["./src/annotation_file/extract_splice_sites.py", annotation_file],stdout=subprocess.PIPE, )
        stdout_value = create_Known_splice_site_file.communicate()[0]
        tem_folder_path = annotation_file_folder_path + "/customized_annotation_file"
        if not os.path.isdir(tem_folder_path):
            os.mkdir(tem_folder_path)
        create_Known_splice_site_file_output = tem_folder_path + "/customized_annotation.splice_site.txt"
        logfile = open(create_Known_splice_site_file_output, 'w')
        logfile.write(stdout_value.decode('utf-8'))
        logfile.close()
        Known_splice_site = create_Known_splice_site_file_output
        GTF_file = annotation_file
        create_ExonLength_file_output = tem_folder_path + "/customized_annotation.ExonLength"
        create_ExonLength_file = subprocess.Popen(["Rscript", "./src/annotation_file/ExonLength_forGene_fromgtf.R", annotation_file, create_ExonLength_file_output])
        create_ExonLength_file.communicate()
        ExonLength_file = create_ExonLength_file_output


    def job_submission(sample_name, strandness, directory, Known_splice_site, GTF_file, ExonLength_file, sequence_type, R1_input="NA", R2_input="NA", unpaired_input = "NA"):

        output_directory = directory + "/output"

        bash_string_0 = '''#! /bin/bash
    
            #BSUB -L /bin/bash
            #BSUB -q short
            #BSUB -n 20 -W 8:00
            #BSUB -o myjob.out
            #BSUB -e myjob.err
            #BSUB -R span[hosts=1]
            #BSUB -R rusage[mem=5000]
    
            '''
        bash_string_1 = 'python {pyscript_path} --sample_name {sample_name} --R1_input {R1_input} --R2_input {R2_input} --unpaired_input {unpaired_input} --strandness {strandness} --directory {directory} --Known_splice_site {Known_splice_site} --GTF_file {GTF_file} --ExonLength_file {ExonLength_file} --sequence_type {sequence_type}'

        pyscript_path= directory + "/src/cluster_job_submission_FPKM_revision.py"

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





    if from_PreProcessing:

        input_df = pd.read_csv(inputCSV, sep=",")

        sample_name_list = input_df["sample_name"].values.tolist()
        Label_list = input_df["Label"].values.tolist()

        Strandness_list = input_df["Strandness"].values.tolist()

        if sequence_type == "Paired":
            R1_input_list = input_df["R1_input"].values.tolist()
            R2_input_list = input_df["R2_input"].values.tolist()
        elif sequence_type == "Single":
            unpaired_input_list = input_df["unpaired_input"].values.tolist()

        Jobid_list = []

        if from_PreProcessing_submission:
            if sequence_type == "Paired":
                for i in range(len(sample_name_list)):
                    sample_name = sample_name_list[i]
                    R1_input = R1_input_list[i]
                    R2_input = R2_input_list[i]
                    strandness = Strandness_list[i]
                    JobID = job_submission(sample_name=sample_name, R1_input=R1_input, R2_input=R2_input, strandness=strandness,
                                           directory=directory, Known_splice_site=Known_splice_site, GTF_file=GTF_file,
                                           ExonLength_file=ExonLength_file, sequence_type=sequence_type)
                    Jobid_list.append(JobID)
            elif sequence_type == "Single":
                for i in range(len(sample_name_list)):
                    sample_name = sample_name_list[i]
                    unpaired_input = unpaired_input_list[i]
                    strandness = Strandness_list[i]
                    JobID = job_submission(sample_name=sample_name, unpaired_input=unpaired_input, strandness=strandness,
                                           directory=directory, Known_splice_site=Known_splice_site, GTF_file=GTF_file,
                                           ExonLength_file=ExonLength_file, sequence_type=sequence_type)
                    Jobid_list.append(JobID)




            JobID_df = pd.DataFrame({"sample" : sample_name_list, "JobID" : Jobid_list})

            JobID_path = output_directory + "/JobID_information.txt"

            JobID_df.to_csv(JobID_path, sep='\t', index=False)

        if not from_PreProcessing_submission:
            JobID_path = output_directory + "/JobID_information.txt"
            JobID_df = pd.read_csv(JobID_path, sep='\t', index_col=None)

        for sample_name in sample_name_list:
            file_path = output_directory + "/" + str(sample_name) + "/done.txt"
            print(str(sample_name) + " is checking")
            while not os.path.exists(file_path):
                time.sleep(30)

        HTSEQ_count_path_list = []
        Gene_FPKM_path_list = []

        err_number = 0
        for sample_name in sample_name_list:
            HTSEQ_count_path = output_directory + "/" + str(sample_name) + "/dump/" + str(sample_name) + "/HTSEQ_count.txt"
            Gene_FPKM_path = output_directory + "/" + str(sample_name) + "/dump/" + str(sample_name) + "/Gene_FPKM.txt"
            if not os.path.isfile(HTSEQ_count_path):
                print(HTSEQ_count_path + " not exist")
                err_number = err_number + 1
            if not os.path.isfile(Gene_FPKM_path):
                print(Gene_FPKM_path + " not exist")
                err_number = err_number + 1
            HTSEQ_count_path_list.append(HTSEQ_count_path)
            Gene_FPKM_path_list.append(Gene_FPKM_path)

        if err_number == 0:
            print("all files are found")


        summary_step1_df = pd.DataFrame({"sample" : sample_name_list, "HTSEQ_count_path" : HTSEQ_count_path_list,
                                         "Gene_FPKM" : Gene_FPKM_path_list, "Label" : Label_list})

        summary_step1_path = output_directory + "/summary_step1.txt"

        summary_step1_df.to_csv(summary_step1_path, sep='\t', index=False, header=False)

        # for sample_name in sample_name_list:
        #     BAM_path_all = output_directory + "/" + str(sample_name) + "/dump/" + str(sample_name) + "/*bam*"
        #     tem_str = "rm -rf " + BAM_path_all
        #     os.system(tem_str)

        print("step_1 done")

    os.chdir(output_directory)

    if from_FPKM_ReadCount_generation:

        if annotation_file == "gencode_v22":  ## ["gencode_v22", "gencode_v29", "gencode_v37", any path to gtf]
            Pipeline_step_2_1_para_3 = annotation_file_folder_path + "/gencode_v22/gencode_v22_ENSG_Type.txt"
        elif annotation_file == "gencode_v29":
            Pipeline_step_2_1_para_3 = annotation_file_folder_path + "/gencode_v29/gencode_v29_ENSG_Type.txt"
        elif annotation_file == "gencode_v37":
            Pipeline_step_2_1_para_3 = annotation_file_folder_path + "/gencode_v37/gencode_v37_ENSG_Type.txt"
        elif annotation_file == "LncBook_Version2.0_all":
            Pipeline_step_2_1_para_3 = annotation_file_folder_path + "/LncBook_Version2.0_all/LncBook_Version2.0_all_ENSG_Type.txt"
        else:
            Pipeline_step_2_1_1_tem = subprocess.Popen(['python', annotation_file_folder_path + "/ENSG_Type.py", annotation_file], stdout=subprocess.PIPE,)
            stdout_value = Pipeline_step_2_1_1_tem.communicate()[0]
            tem_folder_path = annotation_file_folder_path + "/customized_annotation_file"
            if not os.path.isdir(tem_folder_path):
                os.mkdir(tem_folder_path)
            Pipeline_step_2_1_1_tem_output = tem_folder_path + "/customized_ENSG_Type.txt"
            logfile = open(Pipeline_step_2_1_1_tem_output, 'w')
            logfile.write(stdout_value.decode('utf-8'))
            logfile.close()
            Pipeline_step_2_1_para_3 = Pipeline_step_2_1_1_tem_output


        if "microarray" not in gene_type_list:
            Pipeline_step_2_1_para_1 = src_path + '/Expression_Matrix_RNAType_updated.py'
            Pipeline_step_2_1_para_2 = summary_step1_path

            run_Pipeline_step_2_1 = subprocess.Popen(['python', Pipeline_step_2_1_para_1, Pipeline_step_2_1_para_2, Pipeline_step_2_1_para_3])
            run_Pipeline_step_2_1.communicate()

        for gene_type in gene_type_list:
            if gene_type == "protein_coding":
                step2_output_path = summary_step1_path + ".FPKM.ProteinCoding.DEAoutput"
                step2_output_reformat_path = output_directory + "/ProteinCoding_output.txt"
            elif gene_type == "lincRNA":
                step2_output_path = summary_step1_path + ".FPKM.lincRNA.DEAoutput"
                step2_output_reformat_path = output_directory + "/lincRNA_output.txt"

            def step_2_2(summary_step1_path, src_path, gene_type):

                if gene_type == "protein_coding":
                    gene_type_para_a = """'NR==1 || $1=="protein_coding" {print}'"""
                    gene_type_para_b = "ProteinCoding"
                elif gene_type == "lincRNA":
                    gene_type_para_a = """'NR==1 || $1=="lincRNA" {print}'"""
                    gene_type_para_b = "lincRNA"

                run_Pipeline_step_2_2_1_para = summary_step1_path + ".ReadCount"

                tem_str = 'awk '+ gene_type_para_a + ' FS="\t" OFS="\t" ' + run_Pipeline_step_2_2_1_para
                run_Pipeline_step_2_2_1 = subprocess.Popen(tem_str, shell=True, stdout=subprocess.PIPE, )

                run_Pipeline_step_2_2_2 = subprocess.Popen(['cut', "-f", "2-"], stdin=run_Pipeline_step_2_2_1.stdout, stdout=subprocess.PIPE,  )
                stdout_value = run_Pipeline_step_2_2_2.communicate()[0]
                run_Pipeline_step_2_2_2_output = summary_step1_path + ".ReadCount."+ gene_type_para_b
                logfile = open(run_Pipeline_step_2_2_2_output, 'w')
                logfile.write(stdout_value.decode('utf-8'))
                logfile.close()

                run_Pipeline_step_2_2_3_para = summary_step1_path + ".FPKM"

                tem_str = 'awk ' + gene_type_para_a + ' FS="\t" OFS="\t" ' + run_Pipeline_step_2_2_3_para
                run_Pipeline_step_2_2_3 = subprocess.Popen(tem_str, shell=True, stdout=subprocess.PIPE, )

                run_Pipeline_step_2_2_4 = subprocess.Popen(['cut', "-f", "2-"], stdin=run_Pipeline_step_2_2_3.stdout, stdout=subprocess.PIPE,  )
                stdout_value = run_Pipeline_step_2_2_4.communicate()[0]
                run_Pipeline_step_2_2_4_output = summary_step1_path + ".FPKM."+ gene_type_para_b
                logfile = open(run_Pipeline_step_2_2_4_output, 'w')
                logfile.write(stdout_value.decode('utf-8'))
                logfile.close()

                run_Pipeline_step_2_2_5_para_1 = src_path + "/sparse_check.py"
                run_Pipeline_step_2_2_5_para_2 = summary_step1_path + ".ReadCount." + gene_type_para_b
                run_Pipeline_step_2_2_5 = subprocess.Popen(['python', run_Pipeline_step_2_2_5_para_1, run_Pipeline_step_2_2_5_para_2])
                run_Pipeline_step_2_2_5.communicate()


            if gene_type != "microarray":
                # '($3 > 1 || $3 < -1) && $6<0.05 {print $1}' Popen awk do not need " ' "
                # awk_condition = "'($3 > " + step2_output_foldchange_threshold + ' || $3 < -' + step2_output_foldchange_threshold + ') && $' + step2_output_pvalue_column + '<' + step2_output_pvalue_threshold + " {print $1}'"
                # print(awk_condition)
                step_2_2(summary_step1_path, src_path, gene_type)
                print(gene_type + " step_2 done")

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i',"--inputCSV", default=None, type=str, required=True)
    parser.add_argument('-t',"--biomarker_target_gene_type", choices=["protein_coding", "lincRNA", "both"], default="protein_coding", required=True)
    parser.add_argument('-s',"--sequence_type", default=None, type=str, required=True)
    parser.add_argument('-a',"--annotation_file", default="LncBook_Version2.0_all", type=str, required=True)





    args = parser.parse_args()

    # RNASeq_preprocessing(gene_type_list=args.gene_type_list, only_filter_num=args.only_filter_num,
    #                       directory=args.directory, DE_type=args.DE_type,
    #                       filter_requirment_table_path=args.filter_requirment_table_path)

    RNASeq_preprocessing(biomarker_target_gene_type=args.biomarker_target_gene_type, sequence_type=args.sequence_type,
                         annotation_file=args.annotation_file,
                         inputCSV=args.inputCSV_path)

    # stab_selection_main(folder_path=args.folder_path, file_name=args.train_file_name, testfile_name=args.test_file_name, threshold=args.threshold)



if __name__ == '__main__':
    if not sys.warnoptions:
        warnings.simplefilter("ignore")
        os.environ["PYTHONWARNINGS"] = "ignore"  # Also affect subprocesses
    main()
