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
import re
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.lines as mlines

def transformdata(data,transf):
    if transf=="log":
        outdata = data.apply(np.log)
        # Remove inf values
        outdata.replace([np.inf, -np.inf], np.nan,inplace=True)
    elif transf=="arcsinh":
        outdata = data.apply(np.arcsinh)
    elif transf=="sqrt":
        outdata = data.apply(np.sqrt)
    elif transf=="rel":
        # Compute relative data
        outdata = pd.DataFrame(columns=data.columns,index=data.index)
        aa = list(set([s[0:3] for s in outdata.index]))
        for n in aa:
            idx = [n in s for s in [l[0:3] for l in data.index]]
            total = data.loc[idx,:].sum()
            outdata.loc[data.index[idx],:] = data.loc[data.index[idx],:]/total
            iszero = (total==0)
            if any(iszero):
                outdata.loc[data.index[idx],iszero] = 1.0/sum(idx)
        outdata.iloc[:,:] = np.float64(outdata)
    else:
        outdata=data
        
    return outdata

#%% Load data
codus = pd.read_csv("data/CU_SARScoronaviruses2019ALL.tsv", sep="\t",index_col=0)
codons = pd.read_csv("data/codons_table.tab", sep="\t", index_col = 0)
metadata = pd.read_csv("data/SARScoronaviruses2019ALL.csv", index_col=0)

#%% Compute average CU over the genome
codustable = codus.iloc[:,2:].transpose()
codustable.index = [codons.loc[s,"AA"]+s for s in codustable.index]
relbyprot = transformdata(codustable,"rel")
# Average pre virus
acclist = list(set(codus.Accession))
data = pd.DataFrame([relbyprot.loc[:,codus.index[codus.Accession==acc]].mean(axis=1) for acc in acclist], index = acclist)

#%% PCA
pca = PCA(n_components=2)
x = StandardScaler().fit_transform(data)
principalComponents = pca.fit_transform(x)
expl_var = pca.explained_variance_ratio_
features = pd.DataFrame(pca.components_.transpose(), columns=["PCA1","PCA2"], index = data.columns)

## Plot
labels = list(set(metadata.loc[acclist,"Locality"]))
principalDf = pd.DataFrame(principalComponents,columns=["PCA1","PCA2"],index=data.index)
principalDf["source"] = metadata.loc[acclist,"Locality"]
principalDf["date"] = metadata.loc[acclist,"Collection Date"]

#%% Plot cancer places
fig = plt.figure(figsize = (12,10))
ax = plt.subplot() 
ax.set_xlabel(str('Dim 1 (%1.2f%%)' % (expl_var[0]*100)), fontsize = 15)
ax.set_ylabel(str('Dim 2 (%1.2f%%)' % (expl_var[1]*100)), fontsize = 15)
ax.set_title('Dimension Reduction', fontsize = 20)

cmap = plt.cm.get_cmap('jet') # set color map
colors = cmap(np.arange(len(labels))/(len(labels)-1))
handles = [mlines.Line2D([], [], marker='o', color=c, alpha=0.5) for c in colors]

for label, color in zip(labels,colors):
    idx = principalDf["source"]==label
    ax.scatter(principalDf.loc[idx, 'PCA1'], principalDf.loc[idx, 'PCA2'], c = color, s = 50, alpha=0.5)
    for lab in principalDf.index[idx]:
        ax.annotate(lab,[principalDf.loc[lab, 'PCA1'],principalDf.loc[lab, 'PCA2']],fontsize=8)
fig.legend(handles,labels,loc="center right")
ax.grid()

#%% Plot dates
fig = plt.figure(figsize = (15,10))
ax = plt.subplot() 
ax.set_xlabel(str('Dim 1 (%1.2f%%)' % (expl_var[0]*100)), fontsize = 15)
ax.set_ylabel(str('Dim 2 (%1.2f%%)' % (expl_var[1]*100)), fontsize = 15)
ax.set_title('Dimension Reduction', fontsize = 20)

plt.scatter(principalDf.loc[:, 'PCA1'], principalDf.loc[:, 'PCA2'], c = principalDf.loc[:,"date"].rank(), s = 50, alpha=0.5, cmap="jet")
cb = plt.colorbar()

labels = list(set(principalDf["source"]))
for label in labels:
    idx = principalDf["source"]==label
    for lab in principalDf.index[idx]:
        ax.annotate(lab,[principalDf.loc[lab, 'PCA1'],principalDf.loc[lab, 'PCA2']],fontsize=8)
ax.grid()