#!/usr/bin/python
"""
"""

# IMPORT
from MySQLdb import Connect
import pandas as pd
import string
from Bio import Entrez

# VARIABLES
positions = list(string.ascii_uppercase)[:20]
ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species','strain']
Entrez.email = "aaa@bb.cc"

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
            if len(list(sub_df[sub_df[positions[idx]]==each_id].index))>1:
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

def find_species_by_LINgroups(LINgroup_key,LINgroups,taxonomy_df,described_LINgroups):
    members = LINgroups[LINgroup_key]
    sub_df = taxonomy_df.loc[members]
    species_column = sub_df["species"]
    if len(set(list(species_column))) == 1:
        genus = sub_df.iloc[0,5]
        species = sub_df.iloc[0,6]
        taxon = genus + " " + species
        if taxon not in described_LINgroups:
            described_LINgroups[taxon] = [LINgroup_key]
        else:
            described_LINgroups[taxon].append(LINgroup_key)
    return described_LINgroups

def summarize_described_LINgroups_shortest(described_LINgroups):
    described_LINgroups_shortest = {}
    for each_key in described_LINgroups:
        lingroups = described_LINgroups[each_key]
        shortest = lingroups[0]
        for i in lingroups[1:]:
            if len(i) < len(shortest):
                shortest = i
        described_LINgroups_shortest[each_key]=shortest
    return described_LINgroups_shortest

def extract_taxonomy_by_taxid(tax_id):
    handler = Entrez.efetch(db='taxonomy',id=str(tax_id),retmode='xml')
    record = Entrez.read(handler)[0]
    species_name = record["ScientificName"]
    return species_name

def summarize_described_LINgroups_longest(c):
    c.execute("SELECT LIN.Genome_ID, NCBI_Tax_ID.NCBI_Tax_ID, LIN.LIN FROM LIN,NCBI_Tax_ID,Taxonomy "
              "WHERE LIN.Scheme_ID=4 AND LIN.Genome_ID IN "
              "(SELECT DISTINCT(Taxonomy.Genome_ID) FROM Taxonomy WHERE Rank_ID=7) "
              "AND LIN.Genome_ID=Taxonomy.Genome_ID AND NCBI_Tax_ID.Rank_ID=7 "
              "AND Taxonomy.NCBI_Tax_ID=NCBI_Tax_ID.NCBI_Tax_ID "
              "AND NCBI_Tax_ID.NCBI_Tax_ID IN "
              "(SELECT NCBI_Tax_ID FROM Taxonomy WHERE Taxonomy.Rank_ID=7 GROUP BY Taxonomy.NCBI_Tax_ID HAVING count(Taxonomy.NCBI_Tax_ID)>1) "
              "ORDER BY NCBI_Tax_ID.Taxon ASC")
    tmp = c.fetchall()
    Genome_ID = [int(i[0]) for i in tmp]
    species_ID = [int(i[1]) for i in tmp]
    LIN = [i[2].split(",") for i in tmp]
    df = pd.DataFrame()
    df_out = pd.DataFrame()
    df["species_ID"] = species_ID
    df.index=Genome_ID
    for i in range(len(positions)):
        df[positions[i]]=[lin[i] for lin in LIN]
    unique_species_id = list(set(species_ID))
    df_out["species_id"]=unique_species_id
    LINgroup = []
    species_names = []
    for unique_id in unique_species_id:
        species_name = extract_taxonomy_by_taxid(unique_id)
        species_names.append(species_name)
        sub=df[df["species_ID"]==unique_id]
        lingroup=[]
        for i in range(len(positions)):
            position = positions[i]
            if len(set(sub[position])) > 1:
                break
            else:
                lingroup.append(sub.iloc[0,i+1])
        lingroup = ",".join(lingroup)
        LINgroup.append(lingroup)
    df_out["Species"] = species_names
    df_out["LINgroup"]=LINgroup
    return df_out



def main():
    conn,c = connect_to_db("LINdb_NCBI_RefSeq_test")
    taxonomy_df = pd.read_table("taxonomy_to_load.txt",sep="\t",header=0,index_col=0)
    # LIN_df = extract_LINs(c)
    LIN_df = pd.read_table("LIN_df.txt",sep="\t",header=0,index_col=0)
    LINgroups_total = {}
    LINgroups = {}
    A_level_ids = list(set(LIN_df["A"]))
    for each_id in A_level_ids:
        if len(list(LIN_df[LIN_df["A"]==each_id].index))>1:
            LINgroups[each_id] = list(LIN_df[LIN_df["A"]==each_id].index)
    for each_key in LINgroups.keys():
        LINgroups_total[each_key] = LINgroups[each_key]
    for i in range(1,len(positions)):
        LINgroups,LINgroups_total = find_LINgroups(LIN_df,i,LINgroups,LINgroups_total)
    # LINgroups_mc2 = {each_key:LINgroups_total[each_key] for each_key in LINgroups_total.keys() if len(LINgroups_total[each_key])>1}
    LINgroups_keys = sorted(LINgroups_total.keys())
    described_LINgroups = {}
    for each_key in LINgroups_keys:
        described_LINgroups = find_species_by_LINgroups(each_key,LINgroups_total,taxonomy_df,described_LINgroups)
    described_LINgroups_shortest = summarize_described_LINgroups_shortest(described_LINgroups)
    with open("taxonomy_to_LINgroup_shortest.txt","w") as f:
        f.write("Lineage\tLINgroup(s)\n")
        for each_key in described_LINgroups_shortest:
            f.write(each_key)
            f.write("\t")
            f.write(described_LINgroups_shortest[each_key])
            f.write("\n")
    described_LINgroups_longest = summarize_described_LINgroups_longest(c)
    described_LINgroups_longest.to_csv("taxonomy_to_LINgroup_longest.txt",sep="\t",header=True,index=True,index_label="Genome_ID")

# MAIN
if __name__ == '__main__':
    main()