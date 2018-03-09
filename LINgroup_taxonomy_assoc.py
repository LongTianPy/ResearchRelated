#!/usr/bin/python
"""
"""

# IMPORT
from MySQLdb import Connect
import pandas as pd
import string

# VARIABLES
positions = list(string.ascii_uppercase)[:20]
ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species','strain']

# FUNCTIONS
def connect_to_db(db):
    conn = Connect("localhost", "root")
    c = conn.cursor()
    c.execute("use {0}".format(db))
    return conn, c

def extract_LINs(c):
    c.execute("select Genome_ID,LIN from LIN WHERE Scheme_ID=4")
    tmp = c.fetchall()
    Genome_ID = [int(i[0]) for i in tmp]
    LIN = [i[1].split(",") for i in tmp]
    df = pd.DataFrame('0',columns=positions,index=Genome_ID)
    for i in range(len(Genome_ID)):
        genome_id = Genome_ID[i]
        lin = LIN[i]
        for j in range(len(positions)):
            df.loc[genome_id,positions[j]] = lin[j]
    return df

def find_LINgroups(LIN_df,idx,previous_level_LINgroups,LINgroup_total):
    this_level_LINgroups = {}
    for each_previous_LINgroup in previous_level_LINgroups.keys():
        sub_df = LIN_df.loc[previous_level_LINgroups[each_previous_LINgroup]]
        this_level_ids = list(set(sub_df[positions[idx]]))
        for each_id in this_level_ids:
            this_level_LINgroups[str(each_previous_LINgroup)+","+str(each_id)] = list(sub_df[sub_df[positions[idx]]==each_id].index)
    for each_key in this_level_LINgroups.keys():
        LINgroup_total[each_key] = this_level_LINgroups[each_key]
    return this_level_LINgroups,LINgroup_total

def find_shared_lineage(LINgroup_key,LINgroups,taxonomy_df,described_LINgroups):
    members = LINgroups[LINgroup_key]
    sub_df = taxonomy_df.loc[members]
    common_lineage = []
    for rank in ranks:
        if len(set(list(sub_df[rank])))==1:
            common_lineage.append(sub_df.iloc[0,1])
        else:
            break
    common_lineage = ";".join(common_lineage)
    if common_lineage not in described_LINgroups:
        described_LINgroups[common_lineage] = [LINgroup_key]
    else:
        if len(described_LINgroups[common_lineage][0]) < len(LINgroup_key):
            described_LINgroups[common_lineage] = [LINgroup_key]
        else:
            if len(described_LINgroups[common_lineage][0]) == len(LINgroup_key):
                described_LINgroups[common_lineage].append(LINgroup_key)
    return described_LINgroups


def main():
    conn,c = connect_to_db("LINdb_NCBI_RefSeq_test")
    taxonomy_df = pd.read_table("taxonomy_to_load.txt",sep="\t",header=0,index_col=0)
    LIN_df = extract_LINs(c)
    LINgroups_total = {}
    LINgroups = {}
    A_level_ids = list(set(LIN_df["A"]))
    for each_id in A_level_ids:
        LINgroups[each_id] = list(LIN_df[LIN_df["A"]==each_id].index)
    for each_key in LINgroups.keys():
        LINgroups_total[each_key] = LINgroups[each_key]
    for i in range(1,len(positions)):
        LINgroups,LINgroups_total = find_LINgroups(LIN_df,i,LINgroups,LINgroups_total)
    LINgroups_mc2 = {each_key:LINgroups_total[each_key] for each_key in LINgroups_total.keys() if len(LINgroups_total[each_key])>1}
    described_LINgroups = {}
    for each_key in LINgroups_mc2:
        described_LINgroups = find_shared_lineage(each_key,LINgroups_mc2,taxonomy_df,described_LINgroups)
    with open("taxonomy_to_LINgroup.txt","w") as f:
        f.write("Lineage\tLINgroup(s)\n")
        for each_key in described_LINgroups:
            f.write(each_key)
            f.write("\t")
            f.write(" & ".join(described_LINgroups[each_key]))
            f.write("\n")

# MAIN
if __name__ == '__main__':
    main()