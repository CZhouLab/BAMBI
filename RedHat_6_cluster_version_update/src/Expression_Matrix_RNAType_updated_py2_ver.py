#! /usr/bin/python
import sys
import re
from collections import defaultdict
from collections import OrderedDict

input_file=sys.argv[1]
annotation_file=sys.argv[2]

ReadCount_output_file=input_file+".ReadCount"
FPKM_output_file=input_file+".FPKM"

ReadCount_Dict=defaultdict(dict)
FPKM_Dict=defaultdict(dict)
SampleGroup=OrderedDict()

ENSG_Type_Dict={}
f=file(annotation_file,"r")
#f=open("/project/umw_chan_zhou/Data/LncExpDB/LncBook_Version2.0_Gene_Type.txt","r")
for line in f:
        line=line.strip()
        element=line.split("\t")
        ENSG_Type_Dict[element[0]]=element[1]
f.close

fi=file(input_file,"r")
for line in fi:
        line=line.strip()
        element=line.split("\t")
        Sample=element[0]
        ReadCount_file=element[1]
        FPKM_file=element[2]
        Group=element[3]

        Sample_Group=Sample+"_"+Group
        SampleGroup[Sample_Group]=1

        f1=file(ReadCount_file,"r")
        for line_1 in f1:
                line_1=line_1.strip()
                element_1=line_1.split("\t")
                ReadCount_Dict[element_1[0]][Sample_Group]=element_1[1]
        f1.close()

        f2=file(FPKM_file,"r")
        for line_2 in f2:
                line_2=line_2.strip()
                element_2=line_2.split("\t")
                FPKM_Dict[element_2[0]][Sample_Group]=element_2[1]
        f2.close()
fi.close()

fo1=file(ReadCount_output_file,"w")
fo1.write("GeneType\tGene\t"+"\t".join(SampleGroup)+"\n")
fo2=file(FPKM_output_file,"w")
fo2.write("GeneType\tGene\t"+"\t".join(SampleGroup)+"\n")

num_match = 0
num_NA = 0

for g in ReadCount_Dict:
        try:
                g_type=ENSG_Type_Dict[g]
                num_match = num_match + 1
        except Exception:
                g_type="NA"
                num_NA = num_NA + 1
        ReadCount_Value=[]
        FPKM_Value=[]
        for sample in SampleGroup:
                try:
                        readcount_value=ReadCount_Dict[g][sample]
                except Exception:
                        readcount_value=str(0)
                        
                try:
                        fpkm_value=FPKM_Dict[g][sample]
                except Exception:
                        fpkm_value=str(0)

                ReadCount_Value.append(readcount_value)
                FPKM_Value.append(fpkm_value)

        fo1.write(g_type+"\t"+g+"\t"+"\t".join(ReadCount_Value)+"\n")
        fo2.write(g_type+"\t"+g+"\t"+"\t".join(FPKM_Value)+"\n")
fo1.close()
fo2.close()

print("Total matched: " + str(num_match))
print("Total mis_matched: " + str(num_NA))
percentage = str(float(num_match)/(float(num_match) + float(num_NA)))
print("Total match percentage: " + percentage)