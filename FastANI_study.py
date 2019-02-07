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
              "LIN.Scheme_ID=4")
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
    cmd = "sourmash compute {0} -k 9 -n 1000 -o {1} -q"
    for each in final_LINgroups:
        each_working_dir = join(working_dir,each)
        os.mkdir(each_working_dir)
        for each_genome in final_LINgroups[each]:
            filepath = str(df.loc[each_genome,"FilePath"])
            sig_path = join(each_working_dir, "{0}.sig".format(each_genome))
            os.system(cmd.format(filepath, sig_path))

def compare_each_LINgroup(working_dir, final_LINgroups):
    for each in final_LINgroups:
        each_working_dir = join(working_dir, each)
        cmd = "sourmash compare {0}/*.sig -k 9 -o {0}/output.txt --csv {0}/output.csv -q"
        os.system(cmd.format(each_working_dir))


# MAIN
if __name__ == '__main__':
    position = sys.argv[1]
    working_dir = sys.argv[2]
    conn, c = connect_to_db()
    df, final_LINgroups = fetch_genomes_by_LINgroups(position,c)
    create_signatures(working_dir,df,final_LINgroups)
    compare_each_LINgroup(working_dir,final_LINgroups)