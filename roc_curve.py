#!/usr/bin/python
"""
"""

# IMPORT
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import numpy as np
from sklearn.model_selection import cross_val_score
from sklearn import metrics
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
import random
import matplotlib.pyplot as plt
import graphviz
from sklearn import tree
import os
import operator
import sys
import os
import numpy as np
from scipy import interp
import matplotlib.pyplot as plt
from itertools import cycle
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold,KFold
# FUNCTIONS


# MAIN
if __name__ == '__main__':
    depth = int(sys.argv[1])
    identifier = sys.argv[2]
    os.chdir("/Users/longtian/Desktop/LIN/bbmap_eval")
    pyani_df = pd.read_table("pyani_df_notnull.csv", sep=",", index_col=0, header=0)
    pyani_df = pyani_df[(pyani_df["Jaccard_similarity"] == 0)]
    num_same_family = len(pyani_df[pyani_df[identifier] == 1].index)
    num_notsame_family = len(pyani_df[pyani_df[identifier] == 0].index)
    print(num_same_family, num_notsame_family)
    if num_same_family >= num_notsame_family:
        sample_size = num_notsame_family
    else:
        sample_size = num_same_family
    pyani_notfamily_5000 = pyani_df[pyani_df[identifier] == 0].sample(n=sample_size, random_state=0)
    pyani_samefamily_5000 = pyani_df[pyani_df[identifier] == 1].sample(n=sample_size, random_state=0)
    # pyani_df_notnull = pyani_df
    # pyani_df_notnull = pyani_df[pyani_df["Jaccard_similarity"]==0]
    pyani_df_notnull = pyani_notfamily_5000.append(pyani_samefamily_5000)
    features = ["ANI70", "ANI75", "ANI", "cov", "ANI60"]
    features = ["cov", "ANI", "wkid"]
    data_size = len(pyani_df_notnull.index)
    pool = pyani_df_notnull.index
    random_state = np.random.RandomState(0)
    classifier = RandomForestClassifier(n_estimators=250, max_depth=depth, random_state=0)
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0,1,100)
    for i in range(10):
        train_set_index = random.sample(list(pool), round(data_size / 10))
        test_set_index = list(set(pool) - set(train_set_index))
        train_set = pyani_df_notnull.ix[train_set_index]
        test_set = pyani_df_notnull.ix[test_set_index]
        probas_ = classifier.fit(train_set[features],train_set[identifier]).predict_proba(test_set[features])
        fpr, tpr, thresholds = roc_curve(test_set[identifier], probas_[:,1])
        tprs.append(interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_curve)
        plt.plot(fpr, tpr, lw=1, alpha=0.3, label="ROC fold %d (AUC = %0.2f)" % (i, roc_auc))
    plt.plot([0,1],[0,1],linestyle='--',lw=2,color='r',label="Luck",alpha=0.8)
    try:
        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        std_auc = np.std(aucs,axis=0)
        plt.plot(mean_fpr, mean_tpr, color='b',label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),lw=2, alpha=.8)
        std_tpr = np.std(tprs, axis=0)
        tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
        tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
        plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,label=r'$\pm$ 1 std. dev.')
    except:
        pass
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('{0} level Random Forest Model Evaluation, No Jaccard similarity detected, depth={1}'.format("Order",depth))
    plt.legend(loc="lower right")
    plt.savefig("ROC curve {0} depth = {1}_wo_Jaccard".format(identifier,depth))




# X = pyani_df_notnull[features]
# y = pyani_df_notnull[identifier]
# # n_samples, n_features = X.shape
# random_state = np.random.RandomState(0)
# cv = KFold(n_splits=10)
# classifier = RandomForestClassifier(n_estimators=250, max_depth=depth, random_state=0)
# tprs = []
# aucs = []
# mean_fpr = np.linspace(0,1,100)
# i = 0
# for train, test in cv.split(X,y):
#     probas_ = classifier.fit(X.loc[[train],],y.loc[[train],]).predict_proba(X.loc[[test],])
#     fpr, tpr, thresholds = roc_curve(y.loc[[test],],probas_[:,1])
#     tprs.append(interp(mean_fpr, fpr, tpr))
#     tprs[-1][0] = 0.0
#     roc_auc = auc(fpr,tpr)
#     aucs.append(roc_curve)
#     plt.plot(fpr,tpr,lw=1,alpha=0.3,label="ROC fold %d (AUC = %0.2f)"%(i,roc_auc))
#     i += 1
# plt.plot([0,1],[0,1],linestyle='--',lw=2,color='r',label="Luck",alpha=0.8)
# mean_tpr = np.mean(tprs, axis=0)
# mean_tpr[-1] = 1.0
# mean_auc = auc(mean_fpr, mean_tpr)
# std_auc = np.std(aucs)
# plt.plot(mean_fpr, mean_tpr, color='b',label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),lw=2, alpha=.8)
# std_tpr = np.std(tprs, axis=0)
# tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
# tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
# plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,label=r'$\pm$ 1 std. dev.')
# plt.xlim([-0.05, 1.05])
# plt.ylim([-0.05, 1.05])
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('Receiver operating characteristic example')
# plt.legend(loc="lower right")
# plt.savefig("Family Level Random Forest Model Evaluation, depth = {0}".format(depth))