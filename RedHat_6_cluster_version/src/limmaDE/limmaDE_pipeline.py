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

microarray_input=sys.argv[1]
Rscript_path=sys.argv[2]

df = pd.read_csv(microarray_input, sep="\t", index_col=0)
column_name_list = df.columns.tolist()
condition_list = []
for column_name in column_name_list:
    if str(column_name).endswith('_C'):
        condition_list.append(0)
    elif str(column_name).endswith('_T'):
        condition_list.append(1)
    else:
        print("not match")
df_output = pd.DataFrame({"condition": condition_list})
condition_file_path = microarray_input + ".condition.csv"
df_output.to_csv(condition_file_path, sep="\t", header=True, index=False)

limma_output_path = microarray_input + ".DE_tem.output"
limma_test_R = subprocess.Popen(['Rscript', Rscript_path, microarray_input, condition_file_path, limma_output_path])
limma_test_R.communicate()
