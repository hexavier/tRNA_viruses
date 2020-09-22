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
colnames = ["roc_auc","accuracy"]; colnames.extend(list(indf.columns[0:23]))
RFmeans = pd.DataFrame(index = tropisms, columns = colnames)
RFstd = pd.DataFrame(index = tropisms, columns = colnames)

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
    tprs = []; aucs = [] ; precisions = []
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
        acctemp = []; roc_auctemp = []; feattemp = []; tprtemp = []; prectemp = []; prec_auctemp = []
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
            tprtemp.append(interp_tpr)
            prectemp.append(interp_precision)
            roc_auctemp.append(roc_auc_score(ytemp.iloc[test], clf_probs))
            prec_auctemp.append(auc(recall,precision))
            # Accuracy
            y_pred=clf.predict(Xtemp.iloc[test,:])
            acctemp.append(accuracy_score(ytemp.iloc[test],y_pred))
            # Features
            feattemp.append(clf.feature_importances_)
        
        # Model Accuracy, how often is the classifier correct?
        roc_aucs.append(np.mean(roc_auctemp))
        prec_aucs.append(np.mean(prec_auctemp))
        tprs.append(pd.DataFrame(tprtemp).mean(axis = 0))
        precisions.append(pd.DataFrame(prectemp).mean(axis = 0))
        accuracies.append(np.mean(acctemp))
        # Analyze features
        features.loc[:,n] = pd.DataFrame(feattemp).mean(axis=0).values
    
    # Record results
    trop_means = [np.mean(roc_aucs), np.mean(accuracies)]; trop_means.extend(features.mean(axis=1))
    RFmeans.loc[t,:] = trop_means
    trop_std = [np.std(roc_aucs), np.std(accuracies)]; trop_std.extend(features.std(axis=1))
    RFstd.loc[t,:] = trop_std
    
#%% #Creating a bar plot of feature importances
    fig = plt.figure()
    features = features.reindex(features.mean(axis=1).sort_values(ascending=False).index, axis=0)
    sns.barplot(x=features.mean(axis=1), y=features.index, ci=features.std(axis=1))
    plt.errorbar(x= features.mean(axis=1), y= features.index, fmt='none', xerr=features.std(axis=1), ecolor="k")
    # Add labels to your graph
    plt.xlabel('Feature Importance Score')
    plt.ylabel('Features')
    plt.title(str("Visualizing Important Features in %s" % t))
    plt.show()
    fig.savefig(str("plots/RandomForest_VZ/features_proteins_%s.pdf" % t), bbox_inches='tight')
#%% #ROC curve
    fig = plt.figure()
    mean_tpr = np.mean(tprs, axis=0); std_tpr = np.std(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = np.mean(roc_aucs); std_auc = np.std(roc_aucs)
    
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)
    plt.plot(linear_xaxis, mean_tpr, marker='.', color="b", lw = 2,
             label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc))
    plt.fill_between(linear_xaxis, tprs_lower, tprs_upper, color='grey', alpha=.2,
                    label=r'$\pm$ 1 std. dev.')
    # axis labels
    plt.title(str("ROC curve of %s" % t))
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    # show the legend
    plt.legend(loc="lower right")
    # show the plot
    plt.show()
    fig.savefig(str("plots/RandomForest_VZ/roc_proteins_%s.pdf" % t), bbox_inches='tight')
    
#%% #Precision-Recall Curve
    fig = plt.figure()
    mean_prec = np.mean(precisions, axis=0); std_prec = np.std(precisions, axis=0)
    mean_prec[-1] = 0.0
    mean_auc = np.mean(prec_aucs); std_auc = np.std(prec_aucs)
    
    precs_upper = np.minimum(mean_prec + std_prec, 1)
    precs_lower = np.maximum(mean_prec - std_prec, 0)
    
    chance = len(y[y==1]) / len(y)
    plt.plot([0, 1], [chance, chance], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)
    plt.plot(linear_xaxis, mean_prec, marker='.', color="b", lw = 2,
             label=r'Mean Precision-Recall (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc))
    plt.fill_between(linear_xaxis, precs_lower, precs_upper, color='grey', alpha=.2,
                    label=r'$\pm$ 1 std. dev.')
    # axis labels
    plt.title(str("Precision-Recall Curve curve of %s" % t))
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    # show the legend
    plt.legend(loc="lower left")
    # show the plot
    plt.show()
    fig.savefig(str("plots/RandomForest_VZ/prec-recall_proteins_%s.pdf" % t), bbox_inches='tight')

#%% Save
RFmeans.to_csv("results/mean_randomforestVZ_proteins_bytropism.csv")
RFstd.to_csv("results/std_randomforestVZ_proteins_bytropism.csv")
