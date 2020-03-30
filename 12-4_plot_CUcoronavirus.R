library(ggplot2)
library(ggpubr)

# Load SDAs of coronaviruses
showtypes = c("HNSC-oral tongue","HNSC-oral cavity","HNSC-floor of mouth","HNSC-base of tongue","HNSC-larynx","LUSC","LUAD")
SDA = read.csv("results/RTE_CUcoronavirus.csv", row.names = 1)
subclinicSDA = read.csv("results/subclinicRTE_CUcoronavirus.csv", row.names = 1)
allSDA = cbind(SDA,subclinicSDA)
geneexp = read.csv("results/receptors_expression_coronavirus.csv", row.names = 1)

# Create dataset structure
SDAdataset = c()
labels = unique(rownames(geneexp))
for (l in labels){
  # Dataset for SDA
  type = gsub(" |-",".",l)
  dataset_temp = data.frame(row.names = rownames(allSDA))
  dataset_temp$SDA = as.numeric(allSDA[,type])
  dataset_temp$virus = as.character(allSDA$Species)
  dataset_temp$catype = l
  SDAdataset = rbind(SDAdataset,dataset_temp)
}
# Fix order of tissues along the descending respiratory tract
SDAdataset$catype <- factor(SDAdataset$catype,levels = showtypes)

EXPdataset = c()
labels = unique(colnames(geneexp))
for (l in labels){
  # Dataset for gene expression
  dataset_temp = data.frame(row.names = showtypes)
  dataset_temp$expression = as.numeric(geneexp[showtypes,l])/max(as.numeric(geneexp[showtypes,l]))
  dataset_temp$catype = showtypes
  dataset_temp$gene = l
  EXPdataset = rbind(EXPdataset,dataset_temp)
}
# Fix order of tissues along the descending respiratory tract
EXPdataset$catype <- factor(EXPdataset$catype,levels = showtypes)

# Plot
# Set cancer types to plot
ggplot(SDAdataset[SDAdataset$catype %in% showtypes,], aes(x=catype, y=SDA, fill=virus)) +  
  geom_violin(position=position_dodge(1)) +
  theme_classic() +
  stat_compare_means(aes(label = ..p.signif..),na.rm=T, method="wilcox.test") +
  scale_fill_manual(values=c("#f8766d","#51ad00","#5647ea","#ffda00")) + 
  theme(axis.text.x=element_text(size = rel(1), angle=45, hjust=1, vjust =  1)) +
  labs(title="Coronaviruses",y = "SDA")
virussig = compare_means(SDA ~ virus, group.by = c("catype"), method = "wilcox.test", 
              data=SDAdataset[SDAdataset$catype %in% showtypes,],
              ref.group = "Severe acute respiratory syndrome coronavirus 2")
tissuesig = compare_means(SDA ~ catype, method = "wilcox.test", 
                         data=SDAdataset[(SDAdataset$catype %in% showtypes)&(SDAdataset$virus %in% "Severe acute respiratory syndrome coronavirus 2"),])

ggplot(EXPdataset[EXPdataset$catype %in% showtypes,], aes(x=catype, y=expression, group=gene, color=gene)) +  
  geom_line() +
  theme_classic() +
  geom_point() +
  theme(axis.text.x=element_text(size = rel(1), angle=45, hjust=1, vjust =  1)) +
  labs(title="Coronaviruses",y = "Expression")

## Plot the rest of tissues (only show average SDA)
cancer_types = c("HNSC","LUSC","LUAD","CESC","UCEC","BLCA","PRAD","SKCM","BRCA","THCA","ESCA","STAD","COAD","READ","PAAD","LIHC","CHOL","KIRC","KIRP","KICH","PCPG","GBM","THYM")
AVGdataset = c()
viruses = as.character(unique(SDA$Species))
for (l in cancer_types){
  # Dataset for SDA
  type = gsub(" |-",".",l)
  dataset_temp = data.frame(row.names = viruses)
  dataset_temp$SDA = as.numeric(sapply(viruses, function(x) mean(SDA[allSDA$Species %in% x,type])))
  dataset_temp$virus = viruses
  dataset_temp$catype = l
  AVGdataset = rbind(AVGdataset,dataset_temp)
}
# Fix order of tissues along the descending respiratory tract
AVGdataset$catype <- factor(AVGdataset$catype,levels = cancer_types)

ggplot(AVGdataset, aes(x=catype, y=SDA, group=virus, color=virus)) +  
  geom_line() +
  geom_point() +
  theme_classic() +
  scale_color_manual(values=c("#f8766d","#51ad00","#5647ea","#ffda00")) + 
  theme(axis.text.x=element_text(size = rel(1), angle=45, hjust=1, vjust =  1)) +
  labs(title="Coronaviruses",y = "SDA")
