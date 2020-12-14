library(ggplot2)
library(ggpubr)

# Load tAIs and proteomics
TAIs = read.csv("results/RtAI_proteomics_viralinfection.csv", row.names = 1)
prot = read.delim("data/proteomics/proteomics.tsv", sep="\t")
prot = prot[!is.na(prot$RefSeq),]; rownames(prot) = prot$RefSeq

# Map time points between proteomics vs tRNAs
mapping = data.frame(row.names=seq(1,6))
mapping$virus = c("HSV1","HIV1","HIV1","HIV1","HCMV","HCMV")
mapping$trna = c("inf","inf05h","inf12h","inf24h","inf24h","inf72h")
mapping$prot = c("_12h","_6h","_12h","_24h","_24h","_72h")

# Build dataset
cor_dataset = c()
for (s in rownames(mapping)){
  trnatemp = TAIs[,grepl(mapping[s,"virus"],colnames(TAIs))&grepl(mapping[s,"trna"],colnames(TAIs)),drop=F]
  prottemp = prot[,grepl(mapping[s,"virus"],colnames(prot))&grepl(mapping[s,"prot"],colnames(prot)),drop=F]
  # Remove NA
  idx = rowSums(is.na(prottemp))<ncol(prottemp)
  prottemp = prottemp[idx,,drop=F]
  trnatemp = trnatemp[rownames(prottemp),,drop=F]
  
  # Compute correlations
  dataset_temp = data.frame(row.names = seq(1,ncol(prottemp)*ncol(trnatemp)))
  dataset_temp$correlation = as.numeric(cor(cbind(prottemp,trnatemp),method="spearman",use="pairwise")[1:ncol(prottemp),(ncol(prottemp)+1):(ncol(prottemp)+ncol(trnatemp))])
  dataset_temp$time = mapping[s,"trna"]
  dataset_temp$virus = mapping[s,"virus"]
  cor_dataset = rbind(cor_dataset,dataset_temp)
}

# Plot correlations
ggplot(cor_dataset, aes(x=time, y=correlation, color=time)) +
  facet_grid(. ~ virus, scales = "free") +
  geom_jitter() +
  labs(title="Correlations",y = "Spearman Correlation") + 
  theme_classic()
