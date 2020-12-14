library(gplots)
library(ggplot2)
library(ggpubr)

# Load data
sda = read.csv("results/proteomicsSDA.csv", row.names = 1)
prot = read.delim("data/proteomics/proteomics.tsv", sep="\t")
prot = prot[!is.na(prot$RefSeq),]; rownames(prot) = prot$RefSeq
mapping = c(HCMV="CESC;HNSC;SKCM;UCEC",HSV1="CESC;HNSC;SKCM;UCEC",HIV1="THYM",SARS2="LUAD;LUSC",IAV="LUAD;LUSC")

# Compute correlations
viruses = unique(prot$Species)
dataset = c()
for (v in viruses){
  # Proteomics data
  prot_temp = prot[prot$Species %in% v, grep(v,colnames(prot))]
  conds = sapply(colnames(prot_temp),function(x) strsplit(x,"_|\\.")[[1]][2])
  times = sapply(colnames(prot_temp),function(x) strsplit(x,"_|\\.")[[1]][3])
  # SDA data
  tissuecol = rowSums(sapply(strsplit(mapping[v],";")[[1]],function(x) grepl(x,colnames(sda))))>0
  sda_temp = sda[rownames(prot_temp),tissuecol]
  # Compute correlations
  dataset_temp = data.frame(row.names = seq(1,ncol(prot_temp)*ncol(sda_temp)))
  dataset_temp$correlation = as.numeric(cor(cbind(prot_temp,sda_temp),method="spearman",use="pairwise")[1:ncol(prot_temp),(ncol(prot_temp)+1):(ncol(prot_temp)+ncol(sda_temp))])
  dataset_temp$tcga = as.character(sapply(colnames(sda_temp),function(x) rep(x,ncol(prot_temp))))
  dataset_temp$tissue = as.character(sapply(colnames(sda_temp),function(x) rep(strsplit(x, "\\.")[[1]][1],ncol(prot_temp))))
  dataset_temp$time = rep(times,ncol(sda_temp))
  dataset_temp$cond = rep(conds,ncol(sda_temp))
  dataset_temp$virus = v
  # Merge
  dataset = rbind(dataset,dataset_temp)
}
dataset$time <- factor(dataset$time, levels = c("0h","2h","4h","6h","8h","10h","12h","16h","18h","24h","48h","72h"))


# Plot
ggplot(dataset[dataset$cond=="inf",], aes(x=time, y=correlation, fill=time)) +
  facet_wrap(~virus, scales = "free") +
  geom_boxplot(size=1,alpha=0.5) +
  stat_summary(position=position_dodge(1),fun=median, geom="point", shape=23, size=2) + 
  labs(title="SDA vs Prot",y = "Spearman Correlation") + 
  theme_classic()

