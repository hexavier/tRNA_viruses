#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 11:10:15 2020

@author: xhernandez
"""
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import auc

#%% Create dataframe with all features and targets
# Upload data
SDA = pd.read_csv("results/RTE_tissueMeans.csv", index_col=0)
tropism = pd.read_csv("data/virus_list.tsv", sep="\t" ,index_col=0)

# Keep only accession in the ViralZone list
viralzone = pd.read_csv("data/ViralZone_human_viruses.tsv", sep="\t")
VZaccessions = [x for s in viralzone.Genome for x in s.split(", ") if x!="Not available"]
tokeep = [s in VZaccessions for s in SDA.Accession]
SDA = SDA.loc[tokeep,:]

# Create input dataframe with features
indf = pd.DataFrame(SDA.iloc[:,0:23])
# Associate each row to target values
indf.loc[:,"tropism"] = [tropism.loc[s,"tropism"] for s in SDA.Species]

# Remove non-annotated viruses
indf.dropna(subset=["tropism"], inplace=True)

#%% For each tropism, create a RF model that distinguishes one tropism from others based on SDA

tropisms = list(set(indf.tropism))
# Create results object
randomROC_DF = pd.DataFrame(index = range(0,100),columns = tropisms)
randomPR_DF = pd.DataFrame(index = range(0,100),columns = tropisms)
for l in range(0,100):
    indf.loc[:,'tropism'] = list(indf.loc[:,'tropism'].sample(frac=1))
    for t in tropisms:
        X = indf.iloc[:,0:23]  # Features
        y = indf.loc[:,'tropism']==t  # Labels: 1 for t, 0 for other tropisms
        # Check number of target groups
        print("There are %i proteins in %s, and %i proteins in non-%s" % (sum(y),t,(len(y)-sum(y)),t))
        # To make groups comparable, separate equal samples of target and non-target
        roc_aucs = []; prec_aucs = []; 
        n_limiting = min(sum(y),(len(y)-sum(y))) # Check whether limiting group is target or non-target
        for n in range(100):
            # Separe targets from non-targets and sample from them
            targets = indf.loc[y,:].sample(n=n_limiting)
            nontargets = indf.loc[~y,:].sample(n=n_limiting)
            tempdf = targets.append(nontargets)
            # Define features and targets
            Xtemp = tempdf.iloc[:,0:23]
            ytemp = tempdf.loc[:,'tropism']==t
            # Create a Random Forest Classifier
            clf=RandomForestClassifier(n_estimators=100)
            # Cross validation with 5 splits
            cv = StratifiedKFold(n_splits=5)
            roc_auctemp = []; prec_auctemp = []
            for i, (train, test) in enumerate(cv.split(Xtemp, ytemp)):
                clf.fit(Xtemp.iloc[train,:], ytemp.iloc[train])
                # ROC curve and Precision-Recall curve
                # predict probabilities
                clf_probs = clf.predict_proba(Xtemp.iloc[test,:])
                # keep probabilities for the positive outcome only
                clf_probs = clf_probs[:, 1]
                precision, recall, _ = precision_recall_curve(ytemp.iloc[test], clf_probs)
                # Save iteration
                roc_auctemp.append(roc_auc_score(ytemp.iloc[test], clf_probs))
                prec_auctemp.append(auc(recall,precision))
            
            # Model Accuracy, how often is the classifier correct?
            roc_aucs.append(np.mean(roc_auctemp))
            prec_aucs.append(np.mean(prec_auctemp))
        
        # Compute averages
        randomROC_DF.loc[l,t] = np.mean(roc_aucs)
        randomPR_DF.loc[l,t] = np.mean(prec_aucs)

#%% Save
randomROC_DF.to_csv("results/permutationROC_randomforest_VZ.csv")
randomPR_DF.to_csv("results/permutationPR_randomforest_VZ.csv")

