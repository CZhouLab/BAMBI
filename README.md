# BAMBI: Integrative biostatistical and artificial-intelligence models discover coding and noncoding RNA genes as biomarkers

## Introduction

BAMBI (Biostatistics and Artificial-Intelligence integrated methods for Biomarker detection), to identify coding and non-coding RNA biomarkers for diagnosis and prognosis. BAMBI is a comprehensive RNA biomarker detection pipeline which can directly process RNA-seq data or microarray data without additional pre-processing and manual imputation. Compared to other methods, BAMBI can identify a small number of RNA biomarkers while still have a higher prediction performance than other methods. Thus, the putative biomarkers identified by BAMBI would be more easily adapted to clinical practice.

**Authors**: Peng Zhou ([peng.zhou@umassmed.edu](mailto:peng.zhou@umassmed.edu)), Chan Zhou ([chan.zhou@umassmed.edu](mailto:chan.zhou@umassmed.edu))\
**Maintainer**: Peng Zhou

---

## Prerequisites

To use BAMBI, you will need the following programs in your PATH:

- **Python3 (>=3.7.3)**
  - Python packages:
    - `numpy`
    - `pandas`
    - `scikit-learn (>=0.22.1)`
    - `kdepy` (install with: `conda install -c dmentipl kdepy`)
    - `shap` (install with: `conda install -c conda-forge shap`)
- **R (>=3.7.3)**
  - R packages:
    - `edgeR`
    - `limma`
    - `impute`
    - `GenomicFeatures` (if using a customized annotation GTF file)
- **RNA-seq preprocessing software**:
  - `htseq (>=0.6.1)`
  - `pysam (>=0.8.4)`
  - `hisat2 (>=2.0.5)`
  - `sambamba (>=0.6.8)`
- **Operating System**: High-performance computing cluster in Linux (recommended)

---

## Dependency Docker Option

We also provide the Docker option for BAMBI usage. After you install docker, you can pull our image by:

```bash
docker pull billyzh/bambi_public:latest
```

And you can use run this docker image to build a temporary environment(container of Docker) for BAMBI by:

```bash
docker run -it --mount "type=bind,src=%cd%,target=/usr/src/data" billyzh/bambi_public /bin/bash
```

Once the temporary environment is built, you can run BAMBI as described in the "Usage" section.

In the Docker version, all scripts are included in the Docker image.
 
When executing Python scripts, use their full path (e.g., `/usr/src/app`):

For example:

### Examples:

```bash
# General version
python 1.statitical_based_feature_selection_info_generation.py

# Docker version
python /usr/src/app/1.statitical_based_feature_selection_info_generation.py
```
Before you run the RNA-Seq raw data preprocess step(Step 0), use the follow command to activate a relative environment:

```bash
conda activate proprocess_env
```

Before you run the statistical based feature selection or downstream analysis(Step 1 & 2), use the follow command to activate a relative environment:

```bash
conda activate analysis_env
```

## Usage

### General Notes:

- If you want to use your own RNA-Seq or microarray data for BAMBI biomarker detection, you can skip Step 0 (RNA-Seq Preprocessing).
- Example RNA-Seq and microarray tables for testing are provided in the `example_datasets` folder.

---

## Step 0: RNASeq_Preprocessing

This script preprocesses RNA-Seq raw data into FPKM and ReadCount tables.

**Input Requirements:**
A table that includes sequencing file information via the `--inputCSV` argument. The input table should contain:

- `sample_name`
- `Label`: "C" (Control) or "T" (Treatment)
- Sequencing file paths:
  - Single-end: `unpaired_input`
  - Paired-end: `R1_input` and `R2_input`
- `Strandness`: "first", "second", or "unstrand"

**Output Files:** (saved under `output` directory):

- FPKM table: `summary_step1.txt.FPKM.ProteinCoding` or `summary_step1.txt.FPKM.lincRNA`
- ReadCount table: `summary_step1.txt.ReadCount.ProteinCoding` or `summary_step1.txt.ReadCount.lincRNA`
  
**Notes:**

- Since RNA-Seq data preprocessing is time-consuming, this step only works in the high performance computing cluster in Linux
- If you need to preprocess RNA-Seq data, you need to download the annotation folder from follow path, and save it in the `src` folder.
	https://drive.google.com/drive/folders/1534bNkl0DalPEzxiuYDwA_SR0cW4T7UA?usp=sharing

```bash
python 0.RNASeq_preprocessing.py --inputCSV INFO_TABLE_PATH --biomarker_target_gene_type {protein_coding, lincRNA} --sequence_type {Single, Paired} --annotation_file ANNOTATION_NAME            

Arguments:

	--inputCSV [-i]			 # path to your sequencing files information table
	--biomarker_target_gene_type[-t] # target biomarker gene type, "protein_coding" or "lincRNA"
	--sequence_type [-s]		 # sequence files type, "Single" or "Paired"
	--annotation_file [-a] 		 # annotation file usage, ("LncBook_Version2.0_all","gencode_v22", "gencode_v29", "gencode_v37", or any path to your customized gtf)

```

### Examples :
```bash
python 0.RNASeq_preprocessing.py --inputCSV ./0.RNASeq_preprocessing_input_sample_Paired-End.csv --biomarker_target_gene_type protein_coding --sequence_type Paired --annotation_file LncBook_Version2.0_all 
```



## Step 1: statistical based feature selection info generation

This step generates a statistical metrics table for each gene, including:

- `pvalue`: P-value from differential expression analysis.
- `padj`: Adjusted p-value.
- `abs(log2(fold-change))`: Fold-change from analysis.
- `max_val`: Maximum expression value across samples (used for low-expression gene filtering).
- `distribution_overlap_area`: Between-group estimate distribution overlap area.


User can use this table to select thresholds for different metrics for downstream analysis, This table is saved as "./Gene_info.xlsx'.

Output files(under current dir):statistical metrics table named protein_coding_Gene_info.csv or lincRNA_Gene_info.csv

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

After you selected the thresholds for statistical based feature selection, BAMBI will automatically do the follow steps:

•       selected genes based on provided thresholds, and generate corresponding updated gene table

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


