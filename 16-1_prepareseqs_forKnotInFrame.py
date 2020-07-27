# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
def get_codons(seq,frame):
    cutseq = seq[frame:]
    seqcodons = [cutseq[n:n+3] for n in range(0,len(cutseq),3) if len(cutseq)>=(n+3)]
    return seqcodons

#%% Load data
codontab = pd.read_csv("data/codons_table.tab", sep="\t", index_col=0)
# Fasta sequences parsing
text_file = open("data/refseq_humanvirus.fasta", "r")
allfasta = text_file.read().split('>lcl')
text_file.close()
seqs = dict(zip([s.split("\n")[0] for s in allfasta],["".join(s.split("\n")[1:]) for s in allfasta]))

#%% Write seqs without "weird" nucleotides
fastaout = open("data/refseq_humanvirus_forKnotInFrame.fasta","w") 

for s in seqs.keys():
    if not any(["Y" in seqs[s],"M" in seqs[s],"R" in seqs[s],"W" in seqs[s]]):
        # \n is placed to indicate EOL (End of Line) 
        fastaout.write(str(">%s \n" % s))
        fastaout.write(str("%s \n" % seqs[s]))
fastaout.close()

#%% Randomize sequence preserving codons and protein

for n in range(0,5):
    fastaout = open(str("data/refseq_humanvirus_forKnotInFrame_rand%s.fasta" % n),"w") 
    for s in seqs.keys():
        seq = seqs[s]
        if not any(["Y" in seq,"M" in seq,"R" in seq,"W" in seq]):
            # \n is placed to indicate EOL (End of Line) 
            fastaout.write(str(">%s \n" % s))
            
            # Randomize
            seqcodons = get_codons(seq,0)
            seqaa = np.array([codontab.loc[c,"AA"] if c in codontab.index else "nonAA" for c in seqcodons])
            seqrand = np.array(seqcodons)
            for aa in list(set(codontab.AA)):
                seqrand[seqaa==aa] = np.random.permutation(seqrand[seqaa==aa])
            fastaout.write(str("%s \n" % "".join(seqrand)))
    fastaout.close()

