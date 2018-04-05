#!/usr/bin/python
"""
"""

# IMPORT
import os

# VARIABLES
genome_dir = "/home/linproject/Workspace/jellyfish_genomes/"
count_dir = "/home/linproject/Workspace/jellyfish_count/"
countfile_dir = "/home/linproject/Workspace/jellyfish_countfile/"


# FUNCTIONS
def jellyfish_count(fasta):
    prefix = fasta.split(".")[0]
    cmd = "jellyfish count {0}{1}.fasta -m 31 -t 4 -s 1M -C -o {2}{1}.jf".format(genome_dir, prefix, count_dir)
    os.system(cmd)
    return prefix

def jellyfish_dump(prefix):
    cmd = "jellyfish dump {0}{1}.jf -t -c -o {2}{1}.Kmer".format(count_dir,prefix,countfile_dir)
    os.system(cmd)

# MAIN
if __name__ == '__main__':
    genomes = [file for file in os.listdir(genome_dir) if file.endswith("fasta")]
    for genome in genomes:
        prefix = jellyfish_count(genome)
        jellyfish_dump(prefix)