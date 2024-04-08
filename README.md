# BAMBI: Integrative biostatistical and artificial-intelligence models discover coding and noncoding RNA genes as biomarkers

## Introduction

BAMBI (Biostatistics and Artificial-Intelligence integrated methods for Biomarker detection), to identify coding and non-coding RNA biomarkers for diagnosis and prognosis. BAMBI is a comprehensive RNA biomarker detection pipeline which can directly process RNA-seq data or microarray data without pre-processing and manually imputation. Compared to other methods, BAMBI and identify a small number of RNA biomarkers while still have a higher prediction performance than other methods. Thus, the putative biomarkers identified by BAMBI would have more easily adapted in clinical practice.

Authors: Peng Zhou (peng.zhou@umassmed.edu), Chan Zhou (chan.zhou@umassmed.edu)

Maintainer: Peng Zhou


## Prerequisites

To use BAMBI, you will need the following programs in your PATH:

•       python3 (>=3.7.3) 

•       python packages:
•       	numpy
•       	pandas
•       	scikit-learn (>=0.22.1)
•       	kdepy (conda install -c dmentipl kdepy)
•       	shap (conda install -c conda-forge shap)

•       R (>=3.7.3) 

•       R packages:
•       	edgeR
•       	limma
•       	impute
•       	GenomicFeatures (if you need to use customize annotation GTF file)

    
•       Softwares required for RNA-seq data preprocessing
•       htseq (>=0.6.1)
•       pysam (>=0.8.4)
•       hisat2 (>=2.0.5)
•       sambamba (>=0.6.8)

•       OS: high performance computing cluster in Linux (suggested)

## Dependency Docker Option

We also provide the Docker option for BAMBI usage. After you install docker, you can pull our image by:

```bash
docker pull billyzh/bambi_test_0:latest
```

And you can use run this docker image to build a temporary enviroment(container of Docker) for BAMBI by:

```bash
docker run -it --mount "type=bind,src=%cd%,target=/usr/src/data" billyzh/bambi_test_0 /bin/bash
```

After you built the temporary enviroment, you can run BAMBI same as non-Docker option like follow "Usage" Chapter.

In the docker version, all scripts have included in the docker image. 

When executing Python scripts, please use the full path to the script(/usr/src/app):

For example:

### Examples :
```bash
#General version
python 1.statitical_based_feature_selection_info_generation.py

# Docker version
python /usr/src/app/1.statitical_based_feature_selection_info_generation.py
```

Before you run the RNA-Seq raw data proprocess step(Step 0), use the follow command to activate a relative environment:

```bash
conda activate test_env_2
```

Before you run the statitical based feature selection or downstream analysis(Step 1 & 2), use the follow command to activate a relative environment:

```bash
conda activate test_env
```

## Usage

•       If you want to use your own RNA-Seq table or microarray table for BAMBI biomarker detection, you can skip the Step 0 RNA-Seq Preprocess

•       Example RNA-Seq table and microarray table for testing are provided in "example_datasets" folder 

## Step 0: RNASeq_Preprocessing

This script proprocess RNA-Seq raw data into FPKM and ReadCount table

it needs to input a table which includes the sequecning files information in "--inputCSV", relative sample files provided: 

•       "sample_name"

•       Catergory information: "Label" ("C" for Control and "T" for "Treatment")

•       Sequencing file path ("unpaired_input" for Single-End, and "R1_input" & "R2_input" for Paired-End)

•       Strandness information "Strandness" (""first", "second" or "unstrand")

Output files(under "output" dir): FPKM table(summary_step1.txt.FPKM.ProteinCoding or summary_step1.txt.FPKM.lincRNA), ReadCount table(summary_step1.txt.ReadCount.ProteinCoding or summary_step1.txt.ReadCount.lincRNA)

### Remark: 

•       Because RNA-Seq data preprocess is very time consuming, this step only work in the high performance computing cluster in Linux

•       If you need to do the RNASeq Preprocessing, you need to download the annoation folder from follow path, and save it under the src folder
	https://drive.google.com/drive/folders/1534bNkl0DalPEzxiuYDwA_SR0cW4T7UA?usp=sharing

```bash
python 0.RNASeq_preprocessing.py --inputCSV INFO_TABLE_PATH --biomarker_target_gene_type {protein_coding, lincRNA} --sequence_type {Single, Paired} --annotation_file ANNOTATION_NAME            

Arguments:

	--inputCSV [-i]			 # path to your sequecning files information table
	--biomarker_target_gene_type[-t] # target biomarker gene type, "protein_coding" or "lincRNA"
	--sequence_type [-s]		 # sequence files type, "Single" or "Paired"
	--annotation_file [-a] 		 # annotation file usage, ("LncBook_Version2.0_all","gencode_v22", "gencode_v29", "gencode_v37", or any path to your customized gtf)

```

### Examples :
```bash
python 0.RNASeq_preprocessing.py --inputCSV ./0.RNASeq_preprocessing_input_sample_Paired-End.csv --biomarker_target_gene_type protein_coding --sequence_type Paired --annotation_file LncBook_Version2.0_all 
```



## Step 1: statitical based feature selection info generation

it will generate a statitical metrics table for each genes, include: 

•       "pvalue": pvalue from selected differential expression analysis

•       "padj": adjust pvalue from selected differential expression analysis

•       "abs(log2(fold-change))": fold change from Fold Change analysis 

•       "max_val": the maximum expression in all sample, used for Low-Expression Genes Filter  

•       "distribution_overlap_area": between-group estimate distribution overlap area

User can use this table to select thresholds for different metrics for downstream analysis, this table saved in "./Gene_info.xlsx'.

Output files(under current dir):statitical metrics table named protein_coding_Gene_info.csv or lincRNA_Gene_info.csv

### Remark: 

•       If you want to use your own RNA-Seq table, you need to provide both FPKM and ReadCount tables

•       If you want to use your own microarray table, you need to provide microarray tables


```bash
python 1.statitical_based_feature_selection_info_generation.py --biomarker_target_gene_type {protein_coding, lincRNA, microarray} [optional options]           

Arguments:

	--biomarker_target_gene_type [-t]	# target biomarker gene type, "protein_coding" or "lincRNA" or "microarray"
Options:
	--RNASeq_FPKM_table_path [-F]		# if you want to use your own RNA-Seq table, you need to provide FPKM table path here
	--RNASeq_ReadCount_table_path [-R]	# if you want to use your own RNA-Seq table, you need to provide ReadCount table path here
	--microarray_table_path	[-m]		# if you want to use your own microarray table, you need to provide microarray table path here

```

### Examples :
```bash
python 1.statitical_based_feature_selection_info_generation.py --biomarker_target_gene_type protein_coding --RNASeq_FPKM_table_path ./sample_data/FPKM_table.csv --RNASeq_ReadCount_table_path ./sample_data/ReadCount_table.csv
```




## Step 2: downstream analysis

After you selected the thresholds for statitical based feature selection, BAMBI will automatically do the follow steps:

•       selected genes based on provide thresholds, and generate relative update gene table

•       machine learning based feature selection

•       collect results and provide suggested candidate biomarkers

Output files(under "result_summary" dir): individual gene biomarker information(summary_GeneType_high_frequency_gene.csv), gene panel biomarker information(summary_GeneType_selected_models_stat_between_partitions.csv)

```bash
python 2.downstream_analysis.py --biomarker_target_gene_type {protein_coding, lincRNA, microarray} [optional options]           

Arguments:

	--biomarker_target_gene_type [-t]		# target biomarker gene type, "protein_coding" or "lincRNA" or "microarray"
	--target_pvalue_type [-p]			# target differential expression pvalue type for gene filter, "pvalue" or "padj"
	--target_pvalue_threshold [-q]			# target differential expression pvalue threshold for gene filter, (float type, [0, 1], suggest <= 0.05)
	--target_foldchange_threshold [-f]		# target foldchange threshold for gene filter, (float type, [0, inf), suggest 1 or 0.585)
	--target_maxmin_remove_threshold [-m]		# target threshold for Low-Expression Genes Filter, (float type, suggest 1.0 for protein coding, 0.01 for lncRNA)
	--target_overlap_area_threshold [-o]		# target threshold for high distribution overlap Genes Filter, (float type, [0, 1])
	--dataset_name [-n]				# self-defined dataset name


```

### Examples :
```bash
python 2.downstream_analysis.py --biomarker_target_gene_type protein_coding --target_pvalue_type padj --target_pvalue_threshold 0.05 --target_foldchange_threshold 1 --target_maxmin_remove_threshold 1 --target_overlap_area_threshold 0.1 --dataset_name customized_name
```


