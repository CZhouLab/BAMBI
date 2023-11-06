import pandas as pd
import math
import os
import sys



# DE_input_path="summary_step1.txt.ReadCount.lincRNA.DESeq2.output"
# FPKM_input_path="summary_step1.txt.FPKM.lincRNA"


DE_input_path=sys.argv[1]
FPKM_input_path=sys.argv[2]
DE_output_path=sys.argv[3]
microarray=sys.argv[4]


DE_df = pd.read_csv(DE_input_path, sep="\t",index_col=0)
FPKM_df = pd.read_csv(FPKM_input_path, sep="\t",index_col=0)

pseudo_foldchange = float(1e-6)

FPKM_df_transposition = FPKM_df.T
Label_list = []

for i in range(len(FPKM_df_transposition)):
    if str(FPKM_df_transposition.index.values[i]).endswith('_C'):
        Label_list.append(0)
    elif str(FPKM_df_transposition.index.values[i]).endswith('_T'):
        Label_list.append(1)
    else:
        print("Label not match")

if microarray == str(1):
    column_name_list = FPKM_df_transposition.columns.values.tolist()
    for column_name in column_name_list:
        FPKM_df_transposition[column_name] = 2**FPKM_df_transposition[column_name]

FPKM_df_transposition['label'] = Label_list
FPKM_df_transposition_0 = FPKM_df_transposition[FPKM_df_transposition['label'].isin([0])]
FPKM_df_transposition_1 = FPKM_df_transposition[FPKM_df_transposition['label'].isin([1])]

FPKM_df_transposition_0_transform = FPKM_df_transposition_0.drop(['label'], axis=1).T
FPKM_df_transposition_1_transform = FPKM_df_transposition_1.drop(['label'], axis=1).T

FPKM_mean_C_series = FPKM_df_transposition_0_transform.mean(axis=1)
FPKM_mean_T_series = FPKM_df_transposition_1_transform.mean(axis=1)

FPKM_median_C_series = FPKM_df_transposition_0_transform.median(axis=1)
FPKM_median_T_series = FPKM_df_transposition_1_transform.median(axis=1)

FPKM_df_C_info = FPKM_df_transposition_0_transform.copy()
FPKM_df_T_info = FPKM_df_transposition_1_transform.copy()

FPKM_df_C_info['mean_Control'] = FPKM_mean_C_series
FPKM_df_T_info['mean_Treatment'] = FPKM_mean_T_series

FPKM_df_C_info['median_Control'] = FPKM_median_C_series
FPKM_df_T_info['median_Treatment'] = FPKM_median_T_series


# FPKM_df_transposition_0_transform["mean"] = FPKM_df_transposition_0_transform.mean(axis=1)
# FPKM_df_transposition_1_transform["mean"] = FPKM_df_transposition_1_transform.mean(axis=1)
#
#
# FPKM_df_C_mean = FPKM_df_transposition_0_transform[["mean"]].rename(columns={'mean': 'mean_Control'})
# FPKM_df_T_mean = FPKM_df_transposition_1_transform[["mean"]].rename(columns={'mean': 'mean_Treatment'})





FPKM_df_with_info = pd.merge(FPKM_df_C_info, FPKM_df_T_info, left_index=True, right_index=True)

log2_mean_list = []
for i in range(len(FPKM_df_with_info)):
    mean_C = FPKM_df_with_info.iloc[i].at['mean_Control']
    mean_T = FPKM_df_with_info.iloc[i].at['mean_Treatment']
    if mean_C == 0 or mean_T == 0:
        mean_C = mean_C + pseudo_foldchange
        mean_T = mean_T + pseudo_foldchange
    tem_log2_foldchange = math.log(mean_T / mean_C, 2)
    log2_mean_list.append(tem_log2_foldchange)

log2_median_list = []
for i in range(len(FPKM_df_with_info)):
    median_C = FPKM_df_with_info.iloc[i].at['median_Control']
    median_T = FPKM_df_with_info.iloc[i].at['median_Treatment']
    if median_C == 0 or median_T == 0:
        median_C = median_C + pseudo_foldchange
        median_T = median_T + pseudo_foldchange
    tem_log2_foldchange = math.log(median_T / median_C, 2)
    log2_median_list.append(tem_log2_foldchange)

FPKM_df_with_info["log2(fold-change)"] = log2_median_list
FPKM_df_with_info["log2(fold-change)_mean"] = log2_mean_list

info_df = pd.merge(DE_df, FPKM_df_with_info, left_index=True, right_index=True, how='left')

DE_df["log2FoldChange"] = info_df["log2(fold-change)"].tolist()
DE_df["log2FoldChange_mean"] = info_df["log2(fold-change)_mean"].tolist()

tem_output_file = DE_output_path + ".tem"

DE_df.to_csv(tem_output_file, sep="\t", na_rep='NA')

source_fp = open(tem_output_file, 'r')
target_fp = open(DE_output_path, 'w')
first_row = True
for row in source_fp:
    if first_row:
        row = 'baseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\tlog2FoldChange_mean\n'
        first_row = False
    target_fp.write(row)

tem_str = "rm " + tem_output_file
os.system(tem_str)