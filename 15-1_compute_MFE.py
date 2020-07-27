# -*- coding: utf-8 -*-

import pandas as pd
import RNA
import numpy as np
import re

def get_codons(seq,frame):
    cutseq = seq[frame:]
    seqcodons = [cutseq[n:n+3] for n in range(0,len(cutseq),3) if len(cutseq)>=(n+3)]
    return seqcodons

#%% Load data
metadata = pd.read_csv("data/refseq_humanvirus_CoCoPUT.tsv",index_col=0, sep="\t")
codontab = pd.read_csv("data/codons_table.tab", sep="\t", index_col=0)

# Fasta sequences parsing
text_file = open("data/refseq_humanvirus.fasta", "r")
allfasta = text_file.read().split('>')
text_file.close()
seqs = dict(zip([s.split("\n")[0] for s in allfasta],["".join(s.split("\n")[1:]) for s in allfasta]))

#%% Create output table
# Initialize structures
mfe_df = pd.DataFrame(columns = ["WT"], index = metadata.index)
mfe_random_df = pd.DataFrame(columns = np.arange(0,10), index = metadata.index)

# Count codons
for p in metadata.index:
    # Find sequence
    idx = [s for s in seqs.keys() if p in s]
    seq = seqs[idx[0]]
    
    # Compute MFE of WT
    (ss, mfe) = RNA.fold(seq)
    mfe_df.loc[p,"WT"] = mfe
    
    # Randomize sequence preserving codons and protein
    seqcodons = get_codons(seq,0)
    seqaa = np.array([codontab.loc[c,"AA"] if c in codontab.index else "nonAA" for c in seqcodons])
    for n in range(0,10):
        seqrand = np.array(seqcodons)
        for aa in list(set(codontab.AA)):
            seqrand[seqaa==aa] = np.random.permutation(seqrand[seqaa==aa])
        (ss, mfe) = RNA.fold("".join(seqrand))
        mfe_random_df.loc[p,n] = mfe

#%% Save output
output_frame = pd.concat([metadata.iloc[:,0:9],mfe_df,mfe_random_df],axis=1)
output_frame.to_csv("results/MFE_humanvirus2.tsv", sep="\t")
