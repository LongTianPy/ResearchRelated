#!/usr/bin/python
"""
"""

# IMPORT
import pandas as pd
import matplotlib.pyplot as plt
from pandas.stats.api import ols
import sys
import matplotlib as mpl
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amssymb}",
                                     r"\usepackage{amsmath}"]
import multiprocessing as mp
import seaborn as sns; sns.set(color_codes=True)
import os

# FUNCTIONS
def fill_df(job_pair,df):
    idx_count = job_pair[0]
    idx = job_pair[1]
    col = job_pair[2]
    df.loc[idx_count,"est_ANI"] = est_ANI.get_value(idx,col)
    df.loc[idx_count,"wkid"] = wkid.get_value(idx,col)
    if str(idx) == col:
        df.loc[idx_count,"ANI"] = 1
        df.loc[idx_count,"cov"] = 1
        df.loc[idx_count,"aln"] = 1
    else:
        df.loc[idx_count,"ANI"] = pyani_ANI.get_value(idx,col)
        df.loc[idx_count,"cov"] = pyani_cov.get_value(idx,col)
        df.loc[idx_count,"aln"] = pyani_aln.get_value(idx,col)
    return df

# MAIN
if __name__ == '__main__':
    working_dir = "/Users/longtian/Desktop/LIN/bbmap_eval"
    os.chdir(working_dir)
    est_ANI = pd.read_csv("est_ANI.csv", header=0, index_col=0)
    # wkid = pd.read_csv("wkid.csv", header=0, index_col=0)
    # pyani_ANI = pd.read_csv("pyani_ANI.csv", header=0, index_col=0)
    # pyani_cov = pd.read_csv("pyani_cov.csv", header=0, index_col=0)
    # pyani_aln = pd.read_csv("pyani_aln.csv", header=0, index_col=0)
    bbmap_idxs = est_ANI.index
    bbmap_cols = est_ANI.columns
    # pyani_idxs = pyani_ANI.index
    # pyani_cols = pyani_ANI.columns
    total = len(bbmap_cols)^2
    job_pairs = []
    idx_count = 1
    for i in bbmap_idxs:
        for j in bbmap_cols:
            job_pairs.append([idx_count, i, j])
            idx_count += 1
    # pyani_df = pd.DataFrame(columns=["ANI", "cov", "aln", "est_ANI", "wkid"], index=[i + 1 for i in range(total)])
    # for i in job_pairs:
    #     pyani_df = fill_df(i, pyani_df)
    # pyani_df.to_csv("pyani_df.csv")
    # df_pyani_notna = pyani_df[pyani_df["ANI"].notnull()]
    # df_non_empty = df_pyani_notna[df_pyani_notna["wkid"] != 0]
    # df_non_empty.to_csv("pyani_bbmap_non_empty.csv")
    # df_non_empty = pd.read_table("pyani_bbmap_non_empty.csv",index_col=0,header=0,sep=",")
    # ols_ANI_estANI = ols(x=df_non_empty["ANI"], y=df_non_empty["est_ANI"])
    # beta1_ANI_estANI = ols_ANI_estANI.beta.x
    # beta0_ANI_estANI = ols_ANI_estANI.beta.intercept
    # r2_ANI_estANI = ols_ANI_estANI.r2
    # r2_ANI_estANI
    # sns.regplot(x="ANI", y="est_ANI", x_ci='ci',ci=95, data=df_non_empty, line_kws={'color':'red'})
    # plt.xlabel("ANI")
    # plt.ylabel("Estimated ANI transformed from k-mer identity")
    # title = r"Correlation between ANI and estimated ANI, $R^2$={0}".format(r2_ANI_estANI)
    # plt.title(title)
    # plt.savefig("ANI_EstANI.pdf")

    # Look at pairs whose cov > 70%
    # df_cov70 = df_non_empty[df_non_empty["cov"] >= 0.7]
    # ols_ANI_estANI_cov70 = ols(x=df_cov70["ANI"], y=df_cov70["est_ANI"])
    # beta1_ANI_estANI_cov70 = ols_ANI_estANI_cov70.beta.x
    # beta0_ANI_estANI_cov70 = ols_ANI_estANI_cov70.beta.intercept
    # r2_ANI_estANI_cov70 = ols_ANI_estANI_cov70.r2
    # r2_ANI_estANI_cov70
    # print("y={0}x+{1}".format(beta1_ANI_estANI_cov70, beta0_ANI_estANI_cov70))
    # print(r2_ANI_estANI_cov70)
    # sns.regplot(x="ANI",y="est_ANI",data=df_cov70,line_kws={'color':'red'})
    # plt.xlabel("ANI")
    # plt.ylabel("Estimated ANI transformed from k-mer identity")
    # title = r"Correlation between ANI and estimated ANI, overall genome alignment percentage $\geqslant70\%$ , $R^2$={0}".format(r2_ANI_estANI_cov70)
    # plt.title(title)
    # plt.savefig("ANI_EstANI_cov70.pdf")
    #
    pyani_df = pd.read_table("pyani_df.csv",header=0,index_col=0,sep=",")
    metadata = pd.read_table("metadata_include_family.csv", sep=",", index_col=0, header=0)
    index_pyani_df = pyani_df.index
    not_same_family = pd.DataFrame()
    table_idx = []
    ANI = []
    cov = []
    aln = []
    EST_ANI = []
    wkid = []
    query_family = [metadata.loc[i[1], "Family"] for i in job_pairs]
    subject_family = [metadata.loc[int(i[2]), "Family"] for i in job_pairs]
    pyani_df["query_family"] = query_family
    pyani_df["subject_family"] = subject_family
    pyani_df.to_csv("pyani_df_w_family.csv")
    # query_family_chosen = []
    # subject_family_chosen = []
    # for i in job_pairs:
    #     idx = i[0]
    #     if idx in index_pyani_df:
    #         query = i[1]
    #         subject = int(i[2])
    #         if pyani_df.loc[idx, "ANI"] > 0.7 and metadata.loc[query, "Family"] != metadata.loc[subject, "Family"]:
    #             if metadata.loc[query, "Genus"] != metadata.loc[subject, "Genus"]:
    #                 table_idx.append(idx)
    #                 ANI.append(pyani_df.loc[idx,"ANI"])
    #                 cov.append(pyani_df.loc[idx,"cov"])
    #                 aln.append(pyani_df.loc[idx,"aln"])
    #                 EST_ANI.append(pyani_df.loc[idx,"est_ANI"])
    #                 wkid.append(pyani_df.loc[idx,"wkid"])
    #                 query_family_chosen.append(pyani_df.loc[idx,"query_family"])
    #                 subject_family_chosen.append(pyani_df.loc[idx,"subject_family"])
    # not_same_family["ANI"] = ANI
    # not_same_family["cov"] = cov
    # not_same_family["aln"] = aln
    # not_same_family["est_ANI"] = EST_ANI
    # not_same_family["wkid"] = wkid
    # not_same_family["query_family"] = query_family_chosen
    # not_same_family["subject_family"] = subject_family_chosen
    # not_same_family.index = table_idx
    # not_same_family.to_csv("not_same_family.csv")



