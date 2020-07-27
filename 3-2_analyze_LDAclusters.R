library(gplots)

# Load data
species = read.csv("results/virus_RCUs_RelbyProt.csv", row.names = 1)
metadata = read.csv("data/virus_list.tsv", sep="\t", row.names = 1)
classifier = "tropism"

# Keep tropism-specific viruses
data = species[(rownames(species) %in% rownames(metadata)[!is.na(metadata[,classifier])]),]
groups = unique(as.character(metadata[!is.na(metadata[,classifier]),classifier]))

#### Compute average RCU of tropism and cluster
avgRCU = sapply(groups, 
                function(x) apply(data[(rownames(data) %in% (rownames(metadata)[metadata[,classifier] %in% x])),],2,mean))
plot(hclust(dist(t(avgRCU))))
heatmap.2(cor(avgRCU, method = "spearman"),margins = c(11,11))

#### Perform clustering with SDAs, to evaluate correspondence
showtissues = c("GBM", "COAD", "READ", "LIHC", "CHOL","LUAD", "LUSC", "THYM", "CESC", "UCEC", "HNSC", "SKCM")
SDAs = t(read.csv("data/matched_AAcTE.csv",row.names = 1))
tissues = sapply(colnames(SDAs), function(x) strsplit(x, "\\.")[[1]][1])
trna_mean = data.frame(row.names = rownames(SDAs))
for (t in showtissues){
  trna_mean[,t] = rowMeans(SDAs[,tissues %in% t], na.rm = T)
}
plot(hclust(dist(t(trna_mean))))
heatmap.2(cor(trna_mean,method = "spearman"),margins = c(5,5))

#### Perform clustering with gene expression, to evaluate correspondence
showtissues = c("GBM", "COAD", "READ", "LIHC", "CHOL","LUAD", "LUSC", "THYM", "CESC", "UCEC", "HNSC")
genexp = read.csv("data/healthy_mrnaseq.csv",row.names = 1)
mapping = read.csv("data/TSS2catype.tsv", sep="\t", row.names = 1, na.strings = "")
tissues = as.character(sapply(colnames(genexp), function(x) mapping[strsplit(x, "\\.")[[1]][2],"Study.Name"]))
genexp_mean = data.frame(row.names = rownames(genexp))
for (t in showtissues){
  genexp_mean[,t] = rowMeans(genexp[,tissues %in% t], na.rm = T)
}
plot(hclust(dist(t(genexp_mean))))
heatmap.2(cor(genexp_mean, method = "spearman"),margins = c(5,5))
