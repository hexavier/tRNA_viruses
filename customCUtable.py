#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Create CU table
import pandas as pd
import re
# Load files
codoncount = pd.read_csv("data/CU_SARScoronaviruses2019ALL.csv", header = None, index_col=0)

text_file = open("data/CU_SARScoronavirus2019ALL.fa", "r")
fasta = text_file.read().split('\n')
text_file.close()
codingseq = [s for s in fasta if s[0:1]==">"]

codonorder = ["TTT","TTC","TTA","TTG","CTT","CTC","CTA","CTG","ATT","ATC","ATA",
              "ATG","GTT","GTC","GTA","GTG","TAT","TAC","TAA","TAG","CAT","CAC",
              "CAA","CAG","AAT","AAC","AAA","AAG","GAT","GAC","GAA","GAG","TCT",
              "TCC","TCA","TCG","CCT","CCC","CCA","CCG","ACT","ACC","ACA","ACG",
              "GCT","GCC","GCA","GCG","TGT","TGC","TGA","TGG","CGT","CGC","CGA",
              "CGG","AGT","AGC","AGA","AGG","GGT","GGC","GGA","GGG"]

#%% Create codon usage table
# Extract metadata
protein_id = [re.findall("(?<=\[protein_id=).+(?=\.[1-9]\])",s)[0] for s in codingseq]
description = [re.findall("(?<=\[protein=)[0-9A-Za-z \-\(\)]+(?=\] \[)",s)[0] for s in codingseq]
accession = [re.findall("(?<=\|).+(?=\.[1-9]_cds_)",s)[0] for s in codingseq]
# Start structure
codus = pd.DataFrame(index=protein_id,columns=["description","Accession"]+codonorder)
codus.loc[:,"description"] = description
codus.loc[:,"Accession"] = accession

for n,p in enumerate(protein_id):
    protcodons = codoncount.iloc[(64*n+n+1):(64*n+n+65),0]
    codus.loc[p,codonorder] = protcodons[codonorder].values

#%% Write
codus.to_csv("data/CU_SARScoronaviruses2019ALL.tsv",sep="\t")