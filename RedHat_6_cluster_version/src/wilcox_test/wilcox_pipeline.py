import pandas as pd
import os
from io import StringIO
import time
import re
import argparse
import subprocess
import os.path
import sys
from collections import defaultdict
from collections import OrderedDict


# Rscript_path = "/nl/umw_chan_zhou/Billy/Biomarker_Detection_Paper/Wilcoxon_test/test/wilcox_test.R"
# Readcount_input = "/nl/umw_chan_zhou/Billy/Biomarker_Detection_Paper/Wilcoxon_test/summary_step1.txt.ReadCount.lincRNA"


Readcount_input=sys.argv[1]
Rscript_path=sys.argv[2]
microarray_Rscript_path=sys.argv[3]
microarray_two_power_input=sys.argv[4]
microarray=sys.argv[5]

if microarray == str(0):
    df = pd.read_csv(Readcount_input, sep="\t", index_col=0)
elif microarray == str(1):
    df = pd.read_csv(microarray_two_power_input, sep="\t", index_col=0)
column_name_list = df.columns.tolist()
condition_list = []
for column_name in column_name_list:
    if str(column_name).endswith('_C'):
        condition_list.append("Control")
    elif str(column_name).endswith('_T'):
        condition_list.append("Treatment")
    else:
        print("not match")
df_output = pd.DataFrame({"condition": condition_list}).T
condition_file_path = Readcount_input + ".condition.csv"
df_output.to_csv(condition_file_path, sep="\t", header=False, index=False)
wilcox_output_tem_path = Readcount_input + ".DE_tem.output.tem.0"
if microarray == str(0):
    wilcox_test_R = subprocess.Popen(
        ['Rscript', Rscript_path, Readcount_input, condition_file_path, wilcox_output_tem_path])
elif microarray == str(1):
    wilcox_test_R = subprocess.Popen(
        ['Rscript', microarray_Rscript_path, microarray_two_power_input, condition_file_path, wilcox_output_tem_path])

wilcox_test_R.communicate()

df = pd.read_csv(wilcox_output_tem_path, sep = "\t", index_col=0)
df.insert(0, 'baseMean', "NA")
df.insert(2, 'lfcSE', "NA")
df.insert(3, 'stat', "NA")

tem_output_file = Readcount_input + ".DE_tem.output.tem.1"
df.to_csv(tem_output_file, sep="\t")

final_output = Readcount_input + ".DE_tem.output"

source_fp = open(tem_output_file, 'r')
target_fp = open(final_output, 'w')
first_row = True
for row in source_fp:
    if first_row:
        row = 'baseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n'
        first_row = False
    target_fp.write(row)

tem_str = "rm " + tem_output_file
os.system(tem_str)

tem_str = "rm " + wilcox_output_tem_path
os.system(tem_str)
