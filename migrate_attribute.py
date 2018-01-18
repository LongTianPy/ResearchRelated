#!/usr/bin/python
"""
"""

# IMPORT
from MySQLdb import Connect
import pandas as pd
import numpy as np
from Bio import Entrez

# FUNCTIONS
def connect_to_db():
    conn = Connect("localhost", "root")
    c = conn.cursor()
    c.execute("use LINdb_NCBI_RefSeq_test")
    return conn, c



# MAIN
if __name__ == '__main__':
    Entrez.email = "aaa@bbb.ccc"
    conn, c = connect_to_db()
    df = pd.read_table("/home/linproject/Workspace/attribute_table.txt",sep="\t",index_col=0,header=0)
    c.execute("SELECT Rank FROM Taxonomic_ranks")
    tmp = c.fetchall()
    Ranks = [i[0] for i in tmp]
    new_df = pd.DataFrame("N/A",index=df.index,columns=Ranks)
    Genome_IDs = df.index
    for Genome_ID in Genome_IDs:
        if not np.isnan(df.loc[Genome_ID,"Taxonomy ID"]): # If NCBI taxonomy ID provided
            tax_id = str(df.loc[Genome_ID,"Taxonomy ID"])
            handler = Entrez.efetch(db="taxonomy", id=tax_id, rettype="xml")
            record = Entrez.read(handler)[0]
            ParentTaxId = record["ParentTaxId"]
            LineageEx = record["LineageEx"]
            for each_rank in LineageEx[1:]:
                rank_name = each_rank["ScientificName"]
                rank_tax_id = each_rank["TaxId"]
                rank = each_rank["Rank"]
                new_df.loc[Genome_ID,rank] = rank_name
                # c.execute("SELECT EXISTS(SELECT Rank_ID FROM Taxonomic_ranks WHERE Rank='{0}')".format(rank))
                # if c.fetchone()[0] == 1:
                #     c.execute("SELECT Rank_ID, Rank_order FROM Taxonomic_ranks WHERE Rank='{0}'".format(rank))
                #     tmp = c.fetchone()
                #     rank_id = int(tmp[0])
                #     rank_order = int(tmp[1])
                #     c.execute("SELECT EXISTS(SELECT NCBI_Tax_ID FROM NCBI_Tax_ID WHERE NCBI_Tax_ID={0}) AND Rank_order={1}}".format(int(rank_tax_id)),rank_order)
                #     if c.fetchone()[0] == 0:
                #         c.execute("INSERT INTO NCBI_Tax_ID (NCBI_Tax_ID, Taxon, Rank_order) VALUES ({0}, '{1}', {2})".format(rank_tax_id, rank_name, rank_order))
                #         conn.commit()
                #     sql = "INSERT INTO Taxonomy (Genome_ID, Rank_ID, NCBI_Tax_ID) VALUES ({0}, {1}, {2})"
                #     c.execute(sql.format(Genome_ID, rank_id, rank_tax_id))
                #     conn.commit()
            # strain name is always entered by the user him-/her-self
            strain_raw = df.loc[Genome_ID, "Strain"].split(" ")
            strain = []
            for i in strain_raw:
                if i not in strain:
                    strain.append(i)
            strain = " ".join(strain)
            # c.execute("INSERT INTO Taxonomy (Genome_ID, Rank_ID, NCBI_Tax_ID) VALUES ({0},9,{1})".format(Genome_ID,int(tax_id)))
            new_df.loc[Genome_ID,"strain"] = strain
        else:
            genus = df.loc[Genome_ID,"Genus"]
            species = str(genus) + "+" + str(df.loc[Genome_ID,"Species"])
            c.execute("SELECT NCBI_Tax_ID FROM NCBI_Tax_ID WHERE Taxon='{0}'".format(species))
            tmp = c.fetchone()
            if len(tmp) > 0:
                species_tax_id = str(tmp[0])
                handler = Entrez.efetch(db="taxonomy", id=tax_id, rettype="xml")
                record = Entrez.read(handler)[0]
                LineageEx = record["LineageEx"]
                for each_rank in LineageEx[1:]:
                    rank_name = each_rank["ScientificName"]
                    rank_tax_id = each_rank["TaxId"]
                    rank = each_rank["Rank"]
                    # c.execute("SELECT Rank_ID, Rank_order FROM Taxonomic_ranks WHERE Rank='{0}'".format(rank))
                    # tmp = c.fetchone()
                    # rank_id = int(tmp[0])
                    # rank_order = int(tmp[1])
                    # sql = "INSERT INTO Taxonomy (Genome_ID, Rank_ID, NCBI_Tax_ID) VALUES ({0}, {1}, {2})"
                    # c.execute(sql.format(Genome_ID, rank_id, rank_tax_id))
                    # conn.commit()
                    new_df.loc[Genome_ID, rank] = rank_name
            else:
                handler = Entrez.esearch(term=species, db="taxonomy",retmode="xml")
                record = Entrez.read(handler)
                species_tax_id = record["IdList"][0]
                handler = Entrez.efetch(db="taxonomy", id=species_tax_id, rettype="xml")
                record = Entrez.read(handler)[0]
                LineageEx = record["LineageEx"]
                for each_rank in LineageEx[1:]:
                    rank_name = each_rank["ScientificName"]
                    rank_tax_id = each_rank["TaxId"]
                    rank = each_rank["Rank"]
                    # c.execute("SELECT EXISTS(SELECT Rank_ID FROM Taxonomic_ranks WHERE Rank='{0}')".format(rank))
                    # if c.fetchone()[0] == 1:
                    #     c.execute("SELECT Rank_ID, Rank_order FROM Taxonomic_ranks WHERE Rank='{0}'".format(rank))
                    #     tmp = c.fetchone()
                    #     rank_id = int(tmp[0])
                    #     rank_order = int(tmp[1])
                    #     c.execute(
                    #         "SELECT EXISTS(SELECT NCBI_Tax_ID FROM NCBI_Tax_ID WHERE NCBI_Tax_ID={0}) AND Rank_order={1}}".format(
                    #             int(rank_tax_id)), rank_order)
                    #     if c.fetchone()[0] == 0:
                    #         c.execute(
                    #             "INSERT INTO NCBI_Tax_ID (NCBI_Tax_ID, Taxon, Rank_order) VALUES ({0}, '{1}', {2})".format(
                    #                 rank_tax_id, rank_name, rank_order))
                    #         conn.commit()
                    #     sql = "INSERT INTO Taxonomy (Genome_ID, Rank_ID, NCBI_Tax_ID) VALUES ({0}, {1}, {2})"
                    #     c.execute(sql.format(Genome_ID, rank_id, rank_tax_id))
                    #     conn.commit()
                    new_df.loc[Genome_ID, rank] = rank_name
        strain_raw = df.loc[Genome_ID,"strain"].split(" ")
        strain = []
        for i in strain_raw:
            if i not in strain:
                strain.append(i)
        strain = " ".join(strain)
        # c.execute("INSERT INTO Taxonomy (Genome_ID, Rank_ID, Taxon) VALUES ({0}, 9, '{1}')".format(Genome_ID, strain))
        new_df.loc[Genome_ID,"strain"] = strain
    new_df.to_csv("new_db_taxonomy.csv",sep="\t")