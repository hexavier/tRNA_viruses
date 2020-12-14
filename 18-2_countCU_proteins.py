# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import re

def get_codons(seq,frame):
    cutseq = seq[frame:]
    seqcodons = (cutseq[n:n+3] for n in range(0,len(cutseq),3) if len(cutseq)>=(n+3))
    return seqcodons

#%% Load data
codontab = pd.read_csv("data/codons_table.tab", sep="\t", index_col=0)

# Fasta sequences parsing
text_file = open("data/proteomics/protein_sequences.fa", "r")
allfasta = text_file.read().split('>')[1:]
text_file.close()
seqs = dict(zip([s.split("\n")[0] for s in allfasta],["".join(s.split("\n")[1:]) for s in allfasta]))

#%% Create output table
# Find gene names
genes = {g.split("cds_")[1].split(".")[0]:g for g in seqs.keys()}

# Initialize structures
frame = pd.DataFrame(columns = codontab.index, index = genes.keys())

# Count codons
for p in frame.index:
    # Find sequence
    seq = seqs[genes[p]]
    
    # Count codons
    seqcodons = get_codons(seq,0)
    frame.loc[p,:] = 0
    for c in seqcodons:
        if c in frame.columns:
            frame.loc[p,c] += 1

#%% Save output
frame.to_csv("data/proteomics/refseq_proteomics_codoncount.tsv", sep="\t")
