import pandas as pd
import os
import sys

# def spare_check(ReadCount_path):
ReadCount_path = sys.argv[1]
df = pd.read_csv(ReadCount_path, index_col="Gene",sep="\t")

num_0_column = 0
for i in range(len(df)):
    for j in range(len(df.iloc[0, :])):
        if df.iloc[i, :][j] == 0:
            num_0_column += 1
            break
# print(num_0_column)

if num_0_column == len(df):
    revised_df = df + 1
    tem = "mv "+ ReadCount_path + " " + ReadCount_path + "_Ori"
    os.system(tem)
    revised_df.to_csv(ReadCount_path, sep="\t")

