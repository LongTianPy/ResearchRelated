#!/usr/bin/python
"""Based on the formula provided by
"High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries",
when we select 0.001 as the minimal jaccard similarity returned from Sourmash and k=9,
the approximate ANI value (FastANI) (0.7000662131) is very close to our A position
threshold 0.7. Here I'm going to test if this approximate value really correspond to the
BLAST-based ANI value.
What I'm going to do is to select groups of genome where in each group, they only share
the same A position.
"""

# IMPORT
from MySQLdb import Connect
import string
import os
from os.path import isfile, isdir, join
import pandas as pd
import sys
import shutil
from scipy.cluster import hierarchy
import json
import operator

# VARIABLES
sourmash_dir = "/home/linproject/Workspace/Sourmash2.0/all_sketches/"

# FUNCTIONS
def connect_to_db():
    conn = Connect("127.0.0.1", "LINbase", "Latham@537")
    c = conn.cursor()
    c.execute("use LINdb_NCBI_RefSeq_test")
    return conn, c

def fetch_genomes_by_LINgroups(position,c):
    """
    Get LINgroups of genomes
    :param position: length of LINgroup
    :param c: database cursor
    :return: dictionary {LINgroup:[genomes]}
    """
    position = position.upper()
    c.execute("SELECT LIN.Genome_ID, LIN.LIN, Genome.FilePath FROM Genome,LIN WHERE "
              "LIN.Scheme_ID=4 AND LIN.Genome_ID=Genome.Genome_ID")
    tmp = c.fetchall()
    Genome_ID = [int(i[0]) for i in tmp]
    LIN = [i[1] for i in tmp]
    FilePath = [i[2] for i in tmp]
    df = pd.DataFrame()
    df["LIN"] = LIN
    df["FilePath"] = FilePath
    df.index = Genome_ID
    positions = list(string.ascii_uppercase)[:20]
    LINgroup_length = positions.index(position) + 1
    LINgroups = {}
    for i in range(len(Genome_ID)):
        lingroup = ",".join(LIN[i].split(",")[:LINgroup_length+1])
        if lingroup not in LINgroups:
            LINgroups[lingroup] = Genome_ID[i]
    final_LINgroups = {}
    for each_LINgroup in LINgroups.keys():
        desired_lingroup = ",".join(each_LINgroup.split(",")[:-1])
        if desired_lingroup not in final_LINgroups:
            final_LINgroups[desired_lingroup] = [LINgroups[each_LINgroup]]
        else:
            final_LINgroups[desired_lingroup].append(LINgroups[each_LINgroup])
    return df, final_LINgroups

def create_signatures(working_dir, df, final_LINgroups):
    """
    create signatures in target folders at specific k
    :param working_dir:
    :param final_LINgroups:
    :return:
    """
    if not isdir(working_dir):
        os.mkdir(working_dir)
    cmd = "sourmash compute {0} -k 11,15,21,31 -n 2000 -o {1} -q"
    for each in final_LINgroups:
        each_working_dir = join(working_dir,each)
        if not isdir(each_working_dir):
            os.mkdir(each_working_dir)
        for each_genome in final_LINgroups[each]:
            filepath = str(df.loc[each_genome,"FilePath"])
            sig_path = join(each_working_dir, "{0}.sig".format(each_genome))
            os.system(cmd.format(filepath, sig_path))
            shutil.copy(sig_path,working_dir)

def compare_each_LINgroup(working_dir, final_LINgroups):
    for each in final_LINgroups:
        each_working_dir = join(working_dir, each)
        cmd = "sourmash compare {0}/*.sig -k 21 -o {0}/output.txt --csv {0}/output.csv -q"
        os.system(cmd.format(each_working_dir))


def cluster_by_threshold(working_dir,df, threshold):
    print("Clustering genomes at jaccard similarity threshold of {0}".format(threshold))
    samples = list(df.columns)
    sample_distances = {}
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            sample_distances[(samples[i], samples[j])] = 1 - df.loc[samples[i], samples[j]]
    keys = [sorted(k) for k in sample_distances.keys()]
    values = sample_distances.values()
    sorted_keys, distances = zip(*sorted(zip(keys, values)))
    Z = hierarchy.linkage(distances)
    labels = sorted(set([key[0] for key in sorted_keys] + [sorted_keys[-1][-1]]))
    # hierarchy.dendrogram(Z, labels=labels)
    cutree = hierarchy.cut_tree(Z, height=threshold)
    clusters = {}
    for i in range(len(cutree)):
        if str(cutree[i][0]) not in clusters:
            clusters[str(cutree[i][0])] = [labels[i]]
        else:
            clusters[str(cutree[i][0])].append(labels[i])
    with open(join(working_dir,"clusters_{0}.json".format(threshold)),"w") as f:
        json.dump(clusters,f)
    return clusters

def describe_clusters(clusters):
    clusters_sorted = sorted(clusters,key=lambda x:clusters[x])
    for each in clusters_sorted:
        print(each)
        genomes = clusters[each]
        for genome in genomes:
            c.execute("select Genome.Genome_ID,LIN.LIN,NCBI_Tax_ID.Taxon from Genome,LIN,Taxonomy,NCBI_Tax_ID where "
                      "Genome.Genome_ID=LIN.Genome_ID and "
                      "Taxonomy.Genome_ID=Genome.Genome_ID and "
                      "Taxonomy.NCBI_Tax_ID=NCBI_Tax_ID.NCBI_Tax_ID and "
                      "Taxonomy.Rank_ID=6 and "
                      "Genome.FilePath='{0}'".format(genome))
            tmp = c.fetchone()
            genome_id = int(tmp[0])
            lin = tmp[1]
            genus = tmp[2]
            print("{0}\t{2}\t{1}".format(genome_id,lin,genus))

def retrieve_meta(each_cluster,c):
    genome_ids = []
    LINs = []
    for i in each_cluster:
        c.execute('SELECT Genome.Genome_ID,LIN.LIN FROM Genome,LIN WHERE '
                  'Genome.Genome_ID=LIN.Genome_ID AND '
                  'Genome.FilePath="{0}"'.format(i))
        tmp = c.fetchone()
        genome_id = str(tmp[0])
        lin = tmp[1]
        genome_ids.append(genome_id)
        LINs.append(lin)
    for i in range(1,8):
        c.execute('SELECT Taxonomic_ranks.Rank, Taxonomy.NCBI_Tax_ID,NCBI_Tax_ID.Taxon FROM Taxonomy, NCBI_Tax_ID,Taxonomic_ranks WHERE '
                  'Taxonomy.Genome_ID in ({0}) AND '
                  'NCBI_Tax_ID.NCBI_Tax_ID=Taxonomy.NCBI_Tax_ID AND '
                  'Taxonomy.Rank_ID=Taxonomic_ranks.Rank_ID AND '
                  'Taxonomy.Rank_ID={1}'.format(",".join(genome_ids), i))
        tmp = c.fetchall()
        tax_id = [res[1] for res in tmp]
        if len(set(tax_id)) == 1:
            lca = [tmp[0][0],tmp[0][2]]
        else:
            break
    LINs = [i.split(",") for i in LINs]
    i=0
    common_LINgroup = ''
    while i<8:
        LINgroups = [','.join(lin[:i]) for lin in LINs]
        if len(set(LINgroups)) == 1:
            common_LINgroup = LINgroups[0]
            i += 1
        else:
            break
    return lca, common_LINgroup

def validate_clusters(working_dir, df):
    threshold = 0.89
    step = 0.0005
    positions = list(string.ascii_uppercase)[:20]
    while threshold > 0:
        clusters = cluster_by_threshold(working_dir,df,threshold)
        with open(join(working_dir, 'cluster_{0}.txt').format(1-threshold),'w') as f:
            f.write('Threshold\tRank\tTaxon\tLINgroup\tPosition\tSize\n')
            for each in clusters:
                lca, common_LINgroup = retrieve_meta(clusters[each],c)
                position = positions[len(common_LINgroup.split(','))-1]
                # print('{0}\t{1}\t{2}\t{3}\n'.format(1-threshold, lca[0], lca[1], common_LINgroup, position))
                f.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(1-threshold, lca[0], lca[1], common_LINgroup, position,len(clusters[each])))
        threshold = threshold - step





# MAIN
if __name__ == '__main__':
    position = sys.argv[1]
    working_dir = sys.argv[2]
    conn, c = connect_to_db()
    df, final_LINgroups = fetch_genomes_by_LINgroups(position,c)
    create_signatures(working_dir,df,final_LINgroups)
    compare_each_LINgroup(working_dir,final_LINgroups)