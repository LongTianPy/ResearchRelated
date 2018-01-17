#!/usr/bin/python
"""
"""

# IMPORT
import pandas as pd
import sys
import string


# FUNCTIONS
def report_LINgroup(df,positions):
    LINgroups = {}
    for Genome_ID in df.index:
        for i in range(len(positions)-1,0,-1):
            this_LINgroup = ",".join(df.loc[Genome_ID,"LIN"].split(",")[:i])
            if this_LINgroup not in LINgroups:
                LINgroups[this_LINgroup] = [Genome_ID]
            else:
                LINgroups[this_LINgroup].append(Genome_ID)
    return LINgroups


# MAIN
df_file = sys.argv[1]
df = pd.read_table(df_file,sep=",",header=0,index_col=0)
df.sort_values(by=["LIN"])
LINs = [i.split(",") for i in df["LIN"]]
positions = [i for i in string.ascii_uppercase[:20]]
# for i in positions:
#     df[i] = ["0"]*len(df.index)
positions_dict = {i:[] for i in positions}
for i in range(len(LINs)):
    for j in range(len(positions)):
        positions_dict[positions[j]].append(LINs[i][j])
for i in positions:
    df[i] = positions_dict[i]
LINgroups = report_LINgroup(df,positions)
LINgroups = {i:LINgroups[i] for i in LINgroups.keys() if len(LINgroups[i])>1}
LINgroup_to_Description = {}
Description_to_LINgroups = {}
for LINgroup in LINgroups.keys():
    sub_df = df.loc[LINgroups[LINgroup],]
    sub_Lineages = [i.split(";") for i in sub_df["Lineages"]]
    for i in range(1,len(sub_Lineages[0])):
        sub_sub_lineages = [";".join(one_lineage[:i]) for one_lineage in sub_Lineages]
        set_sub_sub_lineages = set(sub_sub_lineages)
        if len(set_sub_sub_lineages) == 1 and i != len(sub_Lineages[0])-1:
            continue
        elif len(set_sub_sub_lineages) == 1 and i == len(sub_Lineages[0])-1:
            described_lineage = ",".join(sub_Lineages[0])
            break
        else:
            described_lineage = ",".join(sub_Lineages[0][:i-1])
            break
    if described_lineage not in Description_to_LINgroups:
        Description_to_LINgroups[described_lineage] = [LINgroup]
    else:
        Description_to_LINgroups[described_lineage].append(LINgroup)
    LINgroup_to_Description[LINgroup] = described_lineage







