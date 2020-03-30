#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 16:09:57 2019

@author: xhernandez
"""

# Load modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.lines as mlines

#%% Load data
species = pd.read_csv("results/virus_RCUs_RelbyProt.csv",index_col=0)
tropism = pd.read_csv("data/virus_list.tsv", sep="\t" ,index_col=0)

#%% PCA virus
data = species.dropna()

pca = PCA(n_components=2)
x = StandardScaler().fit_transform(species)
principalComponents = pca.fit_transform(StandardScaler().fit_transform(x))
expl_var = pca.explained_variance_ratio_
features = pd.DataFrame(pca.components_.transpose(), columns=["PCA1","PCA2"], index = species.columns)

## Plot
principalDf = pd.DataFrame(principalComponents,columns=["PCA1","PCA2"],index=data.index)
 
# Plot pca
fig = plt.figure(figsize = (20,15))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel(str('Dim 1 (%1.2f%%)' % (expl_var[0]*100)), fontsize = 15)
ax.set_ylabel(str('Dim 2 (%1.2f%%)' % (expl_var[1]*100)), fontsize = 15)
ax.set_title('Dimension Reduction', fontsize = 20)

# Color based on tropism
labels = list(set(tropism.loc[:,"tropism"]))
cmap = plt.cm.get_cmap('jet') # set color map
colors = cmap(np.arange(len(labels))/(len(labels)-1))
handles = [mlines.Line2D([], [], marker='o', color=c, alpha=0.5) for c in colors]

for label, color in zip(labels,colors):
    idx = list([s in tropism.loc[tropism.loc[:,"tropism"]==label,:].index for s in principalDf.index])
    ax.scatter(principalDf.loc[idx, 'PCA1'], principalDf.loc[idx, 'PCA2'], c = color, s = 50, alpha=0.5)
    for lab in principalDf.index[idx]:
       ax.annotate(lab,[principalDf.loc[lab, 'PCA1'],principalDf.loc[lab, 'PCA2']],fontsize=8)
fig.legend(handles,labels,loc="center right")
ax.grid()


#%% Save
principalDf.to_csv("results/PCA_relbyprot.csv")
features.to_csv("results/PCA_relbyprot_features.csv")