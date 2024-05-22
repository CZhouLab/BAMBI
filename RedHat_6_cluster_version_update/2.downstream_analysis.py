
import os
import pandas as pd
from io import StringIO
import time
import re
import argparse
from docx import Document
from docx.shared import Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from pdf2image import convert_from_path
import subprocess
import os.path
import shutil
import sys
import warnings



def downstream_analysis(biomarker_target_gene_type="protein_coding", target_pvalue_type = "padj",
                         target_pvalue_threshold = 0.05, target_foldchange_threshold = 1.0,
                         target_maxmin_remove_threshold = 1.0, target_overlap_area_threshold = 0.01,
                        dataset_name="Customized_Dataset", resubmit=False, HPC_parallel = False):
    if not HPC_parallel:
        directory = "/usr/src/app"
        output_directory = os.getcwd() + "/output"
    else:
        directory = os.getcwd()
        output_directory = directory + "/output"

    src_path = directory + "/src"

    statitical_based_feature_selection_gene_filter_script_path = src_path + "/downstream_analysis_scripts/2.statitical_based_feature_selection_gene_filter.py"
    machine_learning_based_feature_selection_10CV_scirpt_path = src_path + "/downstream_analysis_scripts/3_machine_learning_based_feature_selection_10CV.py"
    result_collection_scirpt_path = src_path + "/downstream_analysis_scripts/4_result_collection.py"

    run_Pipeline_tem = subprocess.Popen(['python', statitical_based_feature_selection_gene_filter_script_path,
                                         '--directory', os.getcwd(), #directory,
                                         '--biomarker_target_gene_type', biomarker_target_gene_type,
                                         '--target_pvalue_type', target_pvalue_type,
                                         '--target_pvalue_threshold', str(target_pvalue_threshold),
                                         '--target_foldchange_threshold', str(target_foldchange_threshold),
                                         '--target_maxmin_remove_threshold', str(target_maxmin_remove_threshold),
                                         '--target_overlap_area_threshold', str(target_overlap_area_threshold)])

    run_Pipeline_tem.communicate()

    run_Pipeline_tem = subprocess.Popen(['python', machine_learning_based_feature_selection_10CV_scirpt_path,
                                         '--directory', os.getcwd(), #directory,
                                         '--biomarker_target_gene_type', biomarker_target_gene_type,
                                         '--dataset_name', dataset_name])

    run_Pipeline_tem.communicate()

    run_Pipeline_tem = subprocess.Popen(['python', result_collection_scirpt_path,
                                         '--directory', os.getcwd(), #directory,
                                         '--biomarker_target_gene_type', biomarker_target_gene_type,
                                         '--dataset_name', dataset_name])

    run_Pipeline_tem.communicate()


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t',"--biomarker_target_gene_type", choices=["protein_coding", "lincRNA", "microarray"],
                        default="protein_coding", required=True)
    parser.add_argument('-p',"--target_pvalue_type", default=None, type=str, required=True)
    parser.add_argument('-q', "--target_pvalue_threshold", default=None, type=float, required=True)
    parser.add_argument('-f', "--target_foldchange_threshold", default=None, type=float, required=True)
    parser.add_argument('-m', "--target_maxmin_remove_threshold", default=None, type=float, required=True)
    parser.add_argument('-o', "--target_overlap_area_threshold", default=None, type=float, required=True)
    parser.add_argument('-n', "--dataset_name", default="Customized_Dataset", type=str, required=False)
    parser.add_argument("--r", type=str, default="False", help="whether resubmit termiated job" )
    parser.add_argument('-h', "--HPC", type=str, default="FALSE", required=False)
    args = parser.parse_args()

    def t_or_f(fs):
        ua = str(fs).upper()
        if 'TRUE'.startswith(ua):
            return True
        elif 'FALSE'.startswith(ua):
            return False
        else:
            return True

    # resubmit_requirment = t_or_f(args.r)
    HPC_parallel_opt = t_or_f(args.HPC)
    # Jobstatus_Check(biomarker_target_gene_type=args.biomarker_target_gene_type,
    #                 target_pvalue_type="padj",
    #                 target_pvalue_threshold=0.05, target_foldchange_threshold=1.0,
    #                 target_maxmin_remove_threshold=1.0, target_overlap_area_threshold=0.01,
    #                 dataset_name=args.dataset_name, resubmit=resubmit_requirment)

    downstream_analysis(biomarker_target_gene_type=args.biomarker_target_gene_type,
                        target_pvalue_type=args.target_pvalue_type,
                        target_pvalue_threshold=args.target_pvalue_threshold,
                        target_foldchange_threshold=args.target_foldchange_threshold,
                        target_maxmin_remove_threshold=args.target_maxmin_remove_threshold,
                        target_overlap_area_threshold=args.target_overlap_area_threshold,
                        dataset_name=args.biomarker_target_gene_type, resubmit=False, HPC_parallel=HPC_parallel_opt)

    # pipeline_step_1(sample_name=args.sample_name, R1_input=args.R1_input, R2_input=args.R2_input, strandness=args.strandness, directory=args.directory)

if __name__ == '__main__':
    if not sys.warnoptions:
        warnings.simplefilter("ignore")
        os.environ["PYTHONWARNINGS"] = "ignore"  # Also affect subprocesses
    main()

