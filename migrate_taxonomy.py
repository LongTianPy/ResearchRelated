#!/usr/bin/python
"""
"""

# IMPORT
from MySQLdb import Connect
import pandas as pd
from Bio import Entrez
import requests
from bs4 import BeautifulSoup as bs
import multiprocessing as mp

# VARIABLES
Entrez.email = "aaa@bbb.ccc"
ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species','strain']
ranks_dict = {'superkingdom':1, 'phylum':2, 'class':3, 'order':4, 'family':5, 'genus':6, 'species':7, 'strain':9}

# FUNCTIONS
def connect_to_db(db):
    conn = Connect("localhost", "root")
    c = conn.cursor()
    c.execute("use {0}".format(db))
    return conn, c

def extract_taxonomy_by_tax_id(tax_id,db):
    df = {rank:[] for rank in ranks}
    handler = Entrez.efetch(db='taxonomy',id=tax_id)
    record = Entrez.read(handler)[0]
    if record["Rank"] == "no rank": # Then it is very very likely to be a strain
        lineage_list = record["LineageEx"]
        for each_lineage in lineage_list:
            if each_lineage["Rank"] in df:
                df[each_lineage["Rank"]] = [each_lineage["ScientificName"],each_lineage["TaxId"]]
        df['species'][0] = df['species'][0][len(df['genus'][0])+1:]
        strain_full = record["ScientificName"]
        strain = strain_full[len(df['species'][0])+1:]
        df["strain"] = strain
    else:
        df[record["Rank"]] = [record["ScientificName"],record["TaxId"]]
        lineage_list = record["LineageEx"]
        for each_lineage in lineage_list:
            if each_lineage["Rank"] in df:
                df[each_lineage["Rank"]] = [each_lineage["ScientificName"], each_lineage["TaxId"]]
        strain = db[db["Attribute_ID"]==4]["AttributeValue"][0]
        df["Strain"] = [strain,"N/A"]
    return df

def extract_taxonomy_by_taxon(genome_id,c,db):
    df = {rank:[] for rank in ranks}
    c.execute("select AttributeValue from Attributevalue where Genome_ID={0} and Attribute_ID<5")
    tmp = c.fetchall()
    taxonomy = [i[0] for i in tmp]
    if not taxonomy[1].startswith(taxonomy[0]):
        species = " ".join(taxonomy[:2])
    else:
        species = taxonomy[1]
    handler = Entrez.esearch(db='taxonomy',term=species,retmode='xml')
    record = Entrez.read(handler)
    if record["Count"] != '0':
        species_tax_id = record["IdList"][0]
        # df["species"] = [species,species_tax_id]
        handler2 = Entrez.efetch(db='taxonomy',id=species_tax_id,retmode='xml')
        record2 = Entrez.read(handler2)[0]
        lineage_list = record2["LineageEx"]
        for each_lineage in lineage_list:
            if each_lineage["Rank"] in df:
                df[each_lineage["Rank"]] = [each_lineage["ScientificName"],each_lineage["TaxId"]]
        just_species = df['species'][0][len(df['genus'][0])+1:]
        df['species'][0] = just_species
        strain = db[db["Attribute_ID"]==4]["AttributeValue"][0]
        df['strain'] = [strain,'N/A']
    else:
        genus = taxonomy[0]
    return df

def lookup_tax_id_by_taxon(taxon):
    handler = Entrez.esearch(db='taxonomy',term=taxon,retmode='xml')
    record = Entrez.read(handler)
    return int(record["IdList"][0])

def load_into_db(genome_id,lineage,c,conn):
    """
    Except strain name
    :param genome_id:
    :param lineage:
    :param c:
    :param conn:
    :return:
    """
    for rank in ranks:
        rank_id = ranks_dict[rank]
        value = lineage[rank][0]
        tax_id = lookup_tax_id_by_taxon(value)
        c.execute("select exists(select NCBI_Tax_ID from NCBI_Tax_ID where NCBI_Tax_ID={0})".format(tax_id))
        if c.fetchone()[0] == 0:
            c.execute("insert into NCBI_Tax_ID (NCBI_Tax_ID,Taxon,Rank_ID) VALUES ({0}, '{1}', {2})".format(tax_id,value,rank_id))
            conn.commit()
        c.execute("insert into Taxonomy (Genome_ID,Rank_ID,NCBI_Tax_ID) VALUES ({0},{1},{2})".format(genome_id,rank,tax_id))
        conn.commit()


def get_all_existing(c):
    c.execute("select AttributeValue.Genome_ID,AttributeValue.Attribute_ID,Attribute.AttributeName,AttributeValue.AttributeValue from "
              "AttributeValue,Attribute WHERE AttributeValue.Attribute_ID=Attribute.Attribute_ID")
    tmp = c.fetchall()
    Genome_ID = [int(i[0]) for i in tmp]
    Attribute_ID = [int(i[1]) for i in tmp]
    AttributeName = [i[3] for i in tmp]
    AttributeValue = [str(i[4]) for i in tmp]
    df = pd.DataFrame(columns=[Genome_ID,Attribute_ID,AttributeName,AttributeValue])
    c.execute("select Genome_ID from Genome")
    tmp = c.fetchall()
    genome_id = [int(i[0]) for i in tmp]
    return df, genome_id

def main():
    conn,c = connect_to_db("LINdb_NCBI_RefSeq")
    conn_new,c_new = connect_to_db("LINdb_NCBI_RefSeq_test")
    db_main,Genome_IDs = get_all_existing(c)
    for Genome_ID in Genome_IDs:
        sub_df = db_main[db_main["Genome_ID"]==Genome_ID]
        if sub_df[sub_df["Attribute_ID"]==15]["AttributeValue"][0] != "N/A":
            this_tax_id = int(sub_df[sub_df["Attribute_ID"]==15]["AttributeValue"][0])
            c_new.execute("select exists(select NBI_Tax_ID from NCBI_Tax_ID if NCBI_Tax_ID={0})".format(this_tax_id))
            if c.fetchone()[0] == 0: # While there is still chance that this is not a strain-level tax id
                lineage_df,no_ranks = extract_taxonomy_by_tax_id(this_tax_id)
                load_into_db(Genome_ID,lineage_df,c_new,conn_new)
            strainname = sub_df[sub_df["Attribute_ID"]==4]["AttributeValue"][0]



# def extract_taxonomy(c):
#     c.execute("select AttributeValue.Genome_ID,AttributeValue.Attribute_ID,Attribute.AttributeName,AttributeValue.AttributeValue from "
#               "AttributeValue,Attribute WHERE AttributeValue.Attribute_ID=Attribute.Attribute_ID")
#     tmp = c.fetchall()
#     Genome_ID = [int(i[0]) for i in tmp]
#     Attribute_ID = [int(i[1]) for i in tmp]
#     AttributeName = [i[3] for i in tmp]
#     AttributeValue = [str(i[4]) for i in tmp]
#     df = pd.DataFrame(columns=[Genome_ID,Attribute_ID,AttributeName,AttributeValue])
#     return df

def load_taxonomy(df,c,conn):
    c.execute("select Genome_ID from Genome")
    tmp = c.fetchall()
    Genome_IDs = [int(i[0]) for i in tmp]
    for Genome_ID in Genome_IDs:
        sub_df = df[df["Genome_ID"]==Genome_ID]
        if sub_df[sub_df["Attribute_ID"]==15]["AttributeValue"] != "N/A":
            strain_tax_id = sub_df[sub_df["Attribute_ID"]==15]["AttributeValue"]
            handler = Entrez.efetch(db='taxonomy', id=strain_tax_id,rettype='xml')
            record = Entrez.read(handler)









                if sub_df.loc[i,"Attribute_ID"] != 3:
                    c.execute("select Rank_ID from Taxonomic_ranks WHERE Rank='{0}'".format(str(sub_df.loc[i,"AttributeName"]).lower()))
                    Rank_ID = c.fetchone()[0]
                    taxon = sub_df.loc[i,"AttributeValue"]
                    c.execute("select exists(select NCBI_Tax_ID from NCBI_Tax_ID where Taxon='{0}')".format(taxon))
                    tmp = c.fetchone()[0]
                    if tmp == 0:
                        handler=Entrez.esearch(db='taxonomy',term=taxon,retmode='xml')
                        record=Entrez.read(handler)
                        try:
                            this_tax_id = record["IdList"][0]
                            c.execute("insert into NCBI_Tax_ID (NCBI_Tax_ID,Taxon,Rank_ID) values ({0},'{1}',{2})".format(this_tax_id,taxon,Rank_ID))
                            conn.commit()
                            c.execute("insert into Taxonomy (Genome_ID,Rank_ID,NCBI_Tax_ID) VALUES ({0},{1},{2})".format(
                                Genome_ID, Rank_ID, this_tax_id))
                            conn.commit()
                        except:
                            c.execute("insert into Taxonomy (Genome_ID,Rank_ID,Taxon) values ({0},{1},'{2}')".format(Genome_ID,Rank_ID,taxon))
                    else:
                        c.execute("select NCBI_Tax_ID from NCBI_Tax_ID WHERE Taxon='{0}'".format(taxon))
                        this_tax_id = c.fetchone()[0]
                        c.execute("insert into Taxonomy (Genome_ID,Rank_ID,NCBI_Tax_ID) VALUES ({0},{1},{2})".format(
                                Genome_ID, Rank_ID, this_tax_id))
                        conn.commit()





# MAIN