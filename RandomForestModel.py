#!/usr/bin/python3
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
from sklearn import preprocessing

# FUNCTIONS
def visualize_classifier(model, X, y, ax=None, cmap='rainbow'):
    ax = ax or plt.gca()
    # Plot the training points
    ax.scatter(X[:, 0], X[:, 1], c=y, s=30, cmap=cmap,clim=(y.min(), y.max()), zorder=3)
    ax.axis('tight')
    ax.axis('off')
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    # fit the estimator
    model.fit(X, y)
    xx, yy = np.meshgrid(np.linspace(*xlim, num=200),np.linspace(*ylim, num=200))
    Z = model.predict(np.c_[xx.ravel(), yy.ravel()]).reshape(xx.shape)
    # Create a color plot with the results
    n_classes = len(np.unique(y))
    contours = ax.contourf(xx, yy, Z, alpha=0.3,levels=np.arange(n_classes + 1) - 0.5,cmap=cmap, clim=(y.min(), y.max()),zorder=1)
    ax.set(xlim=xlim, ylim=ylim)

# MAIN
if __name__ == '__main__':
    depth = int(sys.argv[1])
    identifier = sys.argv[2]
    os.chdir("/Users/longtian/Desktop/LIN/bbmap_eval")
    pyani_df = pd.read_table("pyani_df_notnull.csv", sep=",", index_col=0, header=0)
    pyani_df = pyani_df[(pyani_df["Jaccard_similarity"] == 0) & (pyani_df["ANI"]!=0)]
    num_same_family = len(pyani_df[pyani_df[identifier]==1].index)
    num_notsame_family = len(pyani_df[pyani_df[identifier]==0].index)
    print(num_same_family,num_notsame_family)
    if num_same_family >= num_notsame_family:
        sample_size = num_notsame_family
    else:
        sample_size = num_same_family
    pyani_notfamily_5000 = pyani_df[pyani_df[identifier]==0].sample(n=sample_size,random_state=0)
    pyani_samefamily_5000 = pyani_df[pyani_df[identifier]==1].sample(n=sample_size,random_state=0)
    # pyani_df_notnull = pyani_df
    # pyani_df_notnull = pyani_df[pyani_df["Jaccard_similarity"]==0]
    pyani_df_notnull = pyani_notfamily_5000.append(pyani_samefamily_5000)
    features = ["ANI70", "ANI", "cov", "ANI60", "wkid"]
    features = ["cov","ANI","wkid"]
    features = ["ANI","cov"]
    # features = ["cov"]
    data_size = len(pyani_df_notnull.index)
    pool = pyani_df_notnull.index
    enc = preprocessing.OneHotEncoder()
    # pyani_df_notnull = pd.DataFrame(index=pool)
    # for i in features:
    #     pyani_df_notnull[i] = enc.fit_transform(pyani_df[i].reshape(-1,1)).toarray()
    #     # pyani_df_notnull = pyani_df_notnull[pyani_df_notnull["ANI"]<=0.762]

    train_set_index = random.sample(list(pool), round(data_size / 10))
    test_set_index = list(set(pool) - set(train_set_index))
    train_set = pyani_df_notnull.ix[train_set_index]
    test_set = pyani_df_notnull.ix[test_set_index]

    train, test = train_set, test_set
    clf = RandomForestClassifier(n_estimators=250, max_depth=depth, random_state=0)
    clf = clf.fit(train[features],train[identifier])
    scores = [metrics.roc_auc_score(test[identifier], i.predict(test[features])) for i in clf.estimators_]
    index, value = max(enumerate(scores), key=operator.itemgetter(1))
    print(value)

    f = open("NoJaccard_dotfile_randomforest_depth{0}_{1}.dot".format(depth, identifier), "w")
    tree.export_graphviz(clf.estimators_[index], out_file=f,feature_names=features, filled=True,rounded=True)
    f.close()
    os.system("dot {0} -Tpng -o {1}".format("NoJaccard_dotfile_randomforest_depth{0}_{1}.dot".format(depth, identifier),
                                            "NoJaccard_dotfile_randomforest_depth{0}_{1}.png".format(depth, identifier)))