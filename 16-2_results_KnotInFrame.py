# -*- coding: utf-8 -*-

import re
import pandas as pd
import numpy as np

#%% Load data
metadata = pd.read_csv("data/refseq_humanvirus_CoCoPUT.tsv",index_col=0, sep="\t")

# Fasta sequences parsing
text_file = open("results/KnotInFrame_humanviruses.txt", "r")
allfasta = text_file.read().split('>|')[1:]
text_file.close()
seqs = dict(zip([s.split("\n")[0] for s in allfasta],[s.split("\n")[1:] for s in allfasta]))

#%% Write seqs without "weird" nucleotides
proteins = [re.findall("(?<=_cds_)[A-Z0-9_]*(?<=..)",s)[0] for s in seqs.keys()]

output = pd.DataFrame(columns=["knotted_mfe","nested_mfe","delta_mfe"], index = proteins)
for n,p in enumerate(seqs.keys()):
    if len(seqs[p])>2:
        output.iloc[n,0] = float(re.findall("[0-9..-]+",seqs[p][6])[0])
        output.iloc[n,1] = float(re.findall("[0-9..-]+",seqs[p][7])[0])
        output.iloc[n,2] = output.iloc[n,0] - output.iloc[n,1]
    else:
        output.iloc[n,0] = np.nan
        output.iloc[n,1] = np.nan
        output.iloc[n,2] = np.nan

#%% Gather data of randomized sequences
rand_out = pd.DataFrame(columns=["rand1","rand2","rand3","rand4","rand5"], index = proteins)
for rand in range(5):
    text_file = open(str("results/KnotInFrame_humanviruses_rand%i.txt" % rand), "r")
    allfasta = text_file.read().split('>|')[1:]
    text_file.close()
    seqs = dict(zip([s.split("\n")[0] for s in allfasta],[s.split("\n")[1:] for s in allfasta]))

    for n,p in enumerate(seqs.keys()):
        if len(seqs[p])>2:
            knotted = float(re.findall("[0-9..-]+",seqs[p][6])[0])
            nested = float(re.findall("[0-9..-]+",seqs[p][7])[0])
            rand_out.iloc[n,rand] = knotted - nested
        else:
            rand_out.iloc[n,rand] = np.nan

#%% Save
outputALL = pd.concat([metadata.loc[output.index,:].iloc[:,0:9],output,rand_out],axis=1)
outputALL.to_csv("results/KnotInFrame_humanviruses.csv")