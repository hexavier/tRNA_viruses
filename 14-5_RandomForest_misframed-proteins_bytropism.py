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
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import seaborn as sns

#%% Create dataframe with all features and targets
# Upload data
SDA = pd.read_csv("results/frame1RTE_tissueMeans.csv", index_col=0)
tropism = pd.read_csv("data/virus_list.tsv", sep="\t" ,index_col=0)

# Create input dataframe with features
indf = pd.DataFrame(SDA.iloc[:,0:23])
# Associate each row to target values
indf.loc[:,"tropism"] = [tropism.loc[s,"tropism"] for s in SDA.Species]

# Remove non-annotated viruses
indf.dropna(inplace=True)

#%% For each tropism, create a RF model that distinguishes one tropism from others based on SDA

tropisms = list(set(indf.tropism))
# Create results object
colnames = ["roc_auc","accuracy"]; colnames.extend(list(indf.columns[0:23]))
RFmeans = pd.DataFrame(index = tropisms, columns = colnames)

for t in tropisms:
    X = indf.iloc[:,0:23]  # Features
    y = indf.loc[:,'tropism']==t  # Labels: 1 for t, 0 for other tropisms
    # Check number of target groups
    print("There are %i proteins in %s, and %i proteins in non-%s" % (sum(y),t,(len(y)-sum(y)),t))
    # To make groups comparable, separate equal samples of target and non-target
    accuracies = []; roc_aucs = []; prec_aucs = []; 
    features = pd.DataFrame(index=indf.columns[0:23])
    n_limiting = min(sum(y),(len(y)-sum(y))) # Check whether limiting group is target or non-target
    # Initiate structures for ROC curve
    aucs = []
    linear_xaxis = np.linspace(0, 1, 100)
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
        acctemp = []; roc_auctemp = []; prec_auctemp = []
        for i, (train, test) in enumerate(cv.split(Xtemp, ytemp)):
            clf.fit(Xtemp.iloc[train,:], ytemp.iloc[train])
            # ROC curve and Precision-Recall curve
            # predict probabilities
            clf_probs = clf.predict_proba(Xtemp.iloc[test,:])
            # keep probabilities for the positive outcome only
            clf_probs = clf_probs[:, 1]
            fpr, tpr, _ = roc_curve(ytemp.iloc[test], clf_probs)
            precision, recall, _ = precision_recall_curve(ytemp.iloc[test], clf_probs)
            # Interpolate tpr values to match a fixed x axis
            interp_tpr = np.interp(linear_xaxis, fpr, tpr) 
            interp_tpr[0] = 0.0
            interp_precision = np.interp(linear_xaxis, precision, recall) 
            interp_precision[0] = 1.0
            # Save iteration
            roc_auctemp.append(roc_auc_score(ytemp.iloc[test], clf_probs))
            prec_auctemp.append(auc(recall,precision))
            # Accuracy
            y_pred=clf.predict(Xtemp.iloc[test,:])
            acctemp.append(accuracy_score(ytemp.iloc[test],y_pred))
        
        # Model Accuracy, how often is the classifier correct?
        roc_aucs.append(np.mean(roc_auctemp))
        prec_aucs.append(np.mean(prec_auctemp))
        accuracies.append(np.mean(acctemp))
    
    # Record results
    trop_means = [np.mean(roc_aucs), np.mean(accuracies)]; trop_means.extend(features.mean(axis=1))
    RFmeans.loc[t,:] = trop_means
    
#%% Save
RFmeans.to_csv("results/mean_randomforest_frame2proteins_bytropism.csv")
