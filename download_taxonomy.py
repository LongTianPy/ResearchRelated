#!/usr/bin/python
"""
"""

# IMPORT
from MySQLdb import Connect
import pandas as pd
from Bio import Entrez
import sys

# VARIABLES
Entrez.email = "aaa@bb.cc"
ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species','strain']
ranks_dict = {'superkingdom':1, 'phylum':2, 'class':3, 'order':4, 'family':5, 'genus':6, 'species':7, 'strain':9}

# FUNCTIONS
def connect_to_db(db):
    conn = Connect("localhost", "root")
    c = conn.cursor()
    c.execute("use {0}".format(db))
    return conn, c

def extract_all_existing(c):
    df = pd.DataFrame()
    c.execute("Select Genome_ID from Genome")
    tmp = c.fetchall()
    Genome_IDs = [int(i[0]) for i in tmp]
    c.execute("Select AttributeValue_ID, Genome_ID,Attribute_ID,AttributeValue from AttributeValue")
    tmp = c.fetchall()
    AttributeValue_ID = [int(i[0]) for i in tmp]
    Genome_ID = [int(i[1]) for i in tmp]
    Attribute_ID = [int(i[2]) for i in tmp]
    AttributeValue = [i[3] for i in tmp]
    df["Genome_ID"] = Genome_ID
    df["Attribute_ID"] = Attribute_ID
    df["AttributeValue"] = AttributeValue
    df.index = AttributeValue_ID
    df.to_csv("RefSeq_meta.txt", sep="\t",header=True)
    return Genome_IDs, df

def extract_taxonomy_by_taxid(tax_id,db):
    name_list = {rank:[] for rank in ranks}
    handler = Entrez.efetch(db='taxonomy',id=str(tax_id),retmode='xml')
    record = Entrez.read(handler)[0]
    lineage_list = record["LineageEx"]
    if record["Rank"] == "no rank": # Case 1: This tax id points to the individual strain, then
        strain_name_full = record["ScientificName"] # Full name: genus species strain, we dont need it be so long
        for taxon in lineage_list:
            if taxon["Rank"] in ranks:
                name_list[taxon["Rank"]] = [taxon["ScientificName"],taxon["TaxId"]]
        species_name_full = name_list["species"][0]
        genus_name = name_list["genus"][0]
        species_name_simple = species_name_full[len(genus_name)+1:]
        strain_name_simple = strain_name_full[len(species_name_full)+1:]
        name_list["species"][0] = species_name_simple
        name_list["strain"] = [strain_name_simple,tax_id]
    else: # If this given tax id does not represent individual strain (Although I hope this case only happens in RefSeq)
        name_list[record["Rank"]] = [record["ScientificName"], tax_id]
        for taxon in lineage_list:
            if taxon["Rank"] in ranks:
                name_list[taxon["Rank"]] = [taxon["ScientificName"],taxon["TaxId"]]
        if name_list["species"] != []:
            species_name_full = name_list["species"][0]
            genus_name = name_list["genus"][0]
            species_name_simple = species_name_full[len(genus_name) + 1:]
            name_list["species"][0] = species_name_simple
        else:
            species_name_simple = str(db[db["Attribute_ID"==2]]["AttributeValue"])
            name_list["species"] = [species_name_simple,"N/A"]
        try:
            strain_name_simple = str(db[db["Attribute_ID"==4]]["AttributeValue"])
        except:
            strain_name_simple = "N/A"
        name_list["strain"] = [strain_name_simple, "N/A"]
    for i in name_list.keys():
        if name_list[i] == []:
            name_list[i] = ["N/A","N/A"]
    return name_list

def extract_taxonomy_by_entry(db):
    name_list = {rank: [] for rank in ranks}
    given_lineage = db[db["Attribute_ID"]<5]
    genus = given_lineage.iloc[0,2]
    species = given_lineage.iloc[1,2]
    strain = given_lineage.iloc[3,2]
    handler_strain = Entrez.esearch(db='taxonomy',term=' '.join([genus,species,strain]))
    record_strain = Entrez.read(handler_strain)
    if record_strain["Count"] == '0':
        handler_species = Entrez.esearch(db='taxonomy',term = " ".join([genus,species]))
        record_species = Entrez.read(handler_species)
        if record_species["Count"] == '0':
            handler_genus = Entrez.esearch(db='taxonomy',term=genus)
            record_genus = Entrez.read(handler_genus)
            if record_genus["Count"] == '0':
                name_list["genus"] = [genus,"N/A"]
                for i in name_list.keys():
                    if name_list[i] == []:
                        name_list[i] = ["N/A","N/A"]
            else:
                tax_id_genus = record_genus["IdList"][0]
                handler_genus = Entrez.efetch(db='taxonomy',id=tax_id_genus,retmode='xml')
                record_genus = Entrez.read(handler_genus)[0]
                genus = record_genus['ScientificName']
                lineage_list = record_genus["LineageEx"]
                for taxon in lineage_list:
                    if taxon["Rank"] in name_list:
                        name_list[taxon["Rank"]] = [taxon["ScientificName"],taxon["TaxId"]]
                name_list["genus"] = [genus,tax_id_genus]
            name_list["species"] = [species, "N/A"]
            name_list["strain"] = [strain, "N/A"]
        else:
            tax_id_species = record_species["IdList"][0]
            handler_species = Entrez.efetch(db='taxonomy',id=tax_id_species,retmode='xml')
            record_species = Entrez.read(handler_species)[0]
            lineage_list = record_species["LineageEx"]
            for taxon in lineage_list:
                if taxon["Rank"] in name_list:
                    name_list[taxon["Rank"]] = [taxon["ScientificName"], taxon["TaxId"]]
            name_list["species"] = [species,tax_id_species]
        name_list["strain"] = [strain, "N/A"]
    else:
        tax_id_strain = record_species["IdList"][0]
        handler_strain = Entrez.efetch(db='taxonomy',id=tax_id_strain,retmode='xml')
        record_strain = Entrez.read(handler_strain)[0]
        lineage_list = record_strain["LineageEx"]
        for taxon in lineage_list:
            if taxon["Rank"] in name_list:
                name_list[taxon["Rank"]] = [taxon["ScientificName"], taxon["TaxId"]]
        species_name_full = name_list["species"][0]
        genus_name = name_list["genus"][0]
        species_name_simple = species_name_full[len(genus_name) + 1:]
        name_list["species"][0] = species_name_simple
        name_list["strain"] = [strain,tax_id_strain]
    for i in name_list.keys():
        if name_list[i] == []:
            name_list[i] = ["N/A","N/A"]
    return name_list

def write_taxonomy_to_db(c,conn,lineage,genome_id):
    for rank in ranks:
        rank_id = ranks_dict[rank]
        [taxon,id] = lineage[rank]
        if id != 'N/A':
            c.execute("select exists(select NCBI_Tax_ID from NCBI_Tax_ID WHERE NCBI_Tax_ID={0})".format(int(id)))
            if c.fetchone()[0] == 0:
                c.execute("insert into NCBI_Tax_ID (NCBI_Tax_ID,Taxon,Rank_ID) VALUES ({0},'{1}',{2})".format(int(id),taxon,rank_id))
                conn.commit()
            c.execute("insert into Taxonomy (Genome_ID, Rank_ID,NCBI_Tax_ID) VALUES ({0},{1},{2})".format(genome_id,rank_id,int(id)))
            conn.commit()
        else:
            c.execute("insert into Taxonomy (Genome_ID,Rank_ID,Taxon) VALUES ({0},{1},'{2}')".format(genome_id,rank_id,taxon))
            conn.commit()

def main():
    db_original = "LINdb_NCBI_RefSeq"
    db_new = "LINdb_NCBI_RefSeq_test"
    conn_original, c_original = connect_to_db(db_original)
    conn_new, c_new = connect_to_db(db_new)
    Genome_IDs, meta = extract_all_existing(c_original)
    df = pd.DataFrame(columns=ranks,index=Genome_IDs)
    for Genome_ID in Genome_IDs:
        sub_meta = meta[meta["Genome_ID"]==Genome_ID]
        try:
            tax_id_uploaded = int(sub_meta[sub_meta["Attribute_ID"]==15]["AttributeValue"])
            lineage = extract_taxonomy_by_taxid(tax_id_uploaded,sub_meta)
        except:
            lineage = extract_taxonomy_by_entry(sub_meta)
        for rank in ranks:
            try:
                df.loc[Genome_ID,rank] = lineage[rank][0]
            except:
                print(rank)
                print(lineage[rank])
                sys.exit()
        df.to_csv("taxonomy_to_load.txt",sep='\t',header=True,index=True)
        # write_taxonomy_to_db(c_new,conn_new,lineage,Genome_ID)


# MAIN
if __name__ == '__main__':
    main()