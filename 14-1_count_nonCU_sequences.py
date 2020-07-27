# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import re

def get_codons(seq,frame):
    cutseq = seq[frame:]
    seqcodons = (cutseq[n:n+3] for n in range(0,len(cutseq),3) if len(cutseq)>=(n+3))
    return seqcodons

def get_dinucleotides(seq):
    seqdinucl = (seq[n:n+2] for n in range(0,len(seq),2) if len(seq)>=(n+2))
    return seqdinucl

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
frame_plus1 = pd.DataFrame(columns = codontab.index, index = metadata.index)
frame_plus2 = pd.DataFrame(columns = codontab.index, index = metadata.index)
dinucleotides = pd.DataFrame(columns = ["AT","AG","AC","AA",
                                        "TT","TG","TC","TA",
                                        "GT","GG","GC","GA",
                                        "CT","CG","CC","CA"], index = metadata.index)

# Count codons
for p in metadata.index:
    # Find sequence
    idx = [s for s in seqs.keys() if p in s]
    seq = seqs[idx[0]]

    # Count codons for frame +1
    seqcodons_plus1 = get_codons(seq,1)
    frame_plus1.loc[p,:] = 0
    for c in seqcodons_plus1:
        if c in frame_plus1.columns:
            frame_plus1.loc[p,c] += 1
    
    # Count codons for frame +2
    seqcodons_plus2 = get_codons(seq,2)
    frame_plus2.loc[p,:] = 0
    for c in seqcodons_plus2:
        if c in frame_plus2.columns:
            frame_plus2.loc[p,c] += 1
    
    # Count dinucleotides
    seqdinucl = get_dinucleotides(seq)
    dinucleotides.loc[p,:] = 0
    for c in seqdinucl:
        if c in dinucleotides.columns:
            dinucleotides.loc[p,c] += 1

#%% Save output
output_frame1 = pd.concat([metadata.iloc[:,0:9],frame_plus1],axis=1)
output_frame1.to_csv("data/refseq_humanvirus_frame1.tsv", sep="\t")

output_frame2 = pd.concat([metadata.iloc[:,0:9],frame_plus2],axis=1)
output_frame2.to_csv("data/refseq_humanvirus_frame2.tsv", sep="\t")

output_dinucl = pd.concat([metadata.iloc[:,0:9],dinucleotides],axis=1)
output_dinucl.to_csv("data/refseq_humanvirus_dinucleotides.tsv", sep="\t")
