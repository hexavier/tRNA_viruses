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
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.preprocessing import StandardScaler
import matplotlib.lines as mlines

#%% Load data
species = pd.read_csv("results/virus_RCUs_RelbyProt.csv",index_col=0)
tropism = pd.read_csv("data/virus_list.tsv", sep="\t" ,index_col=0)

#%% Dimension reduction
data = species.dropna()
keep = [s for s in tropism.dropna(subset=["tropism"]).index if s in data.index]

## Discriminant analysis
targets = tropism.loc[keep,"tropism"]
data = species.loc[keep,:]
   
## Discriminant analysis
disc = LinearDiscriminantAnalysis(n_components=2, store_covariance=True)
x = StandardScaler().fit_transform(data)
principalComponents = disc.fit_transform(x,targets)
class_score = disc.score(x,targets)
expl_var = disc.explained_variance_ratio_
normcoef = disc.coef_
coef = pd.DataFrame(normcoef.transpose(), index = data.columns)

## Plot
principalDf = pd.DataFrame(principalComponents,columns=["PCA1","PCA2"],index=data.index)
 
# Plot pca
fig = plt.figure(figsize = (20,15))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel(str('Dim 1 (%1.2f%%)' % (expl_var[0]*100)), fontsize = 15)
ax.set_ylabel(str('Dim 2 (%1.2f%%)' % (expl_var[1]*100)), fontsize = 15)
ax.set_title(str('Dimension Reduction (score=%1.3f)' % class_score), fontsize = 20)

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
principalDf.to_csv("results/DiscriminantAnalysis_relbyprot_tropism.csv")
coef.to_csv("results/DiscriminantAnalysis_relbyprot_tropism_features.csv")