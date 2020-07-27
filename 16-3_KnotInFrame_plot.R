library(ggplot2)
library(ggpubr)

# Load data
knotinframe = read.csv("results/KnotInFrame_humanviruses.csv", row.names = 1)
metadata = read.csv("data/virus_list.tsv",sep="\t",row.names = 1)

# Reestructure data in rows
dataset1 = data.frame(row.names = rownames(knotinframe))
dataset1$deltaMFE = knotinframe$delta_mfe
dataset1$sequence = "real"
dataset1$tropism = factor(as.character(metadata[knotinframe$Species,"tropism"]),
                         levels = c("Immune cells","Epithelial cells","Neurons","Respiratory tract","Hepatocytes","Intestine"))
# Random MFE
dataset2 = data.frame(row.names = rownames(knotinframe))
dataset2$deltaMFE = rowMeans(knotinframe[,c("rand1","rand2","rand3","rand4","rand5")],na.rm=T)
dataset2$sequence = "mean_random"
dataset2$tropism = factor(as.character(metadata[knotinframe$Species,"tropism"]),
                          levels = c("Immune cells","Epithelial cells","Neurons","Respiratory tract","Hepatocytes","Intestine"))
# Create full dataset
isna = (!is.na(dataset1$deltaMFE))&(!is.na(dataset2$deltaMFE))
dataset = rbind(dataset1[isna,],dataset2[isna,])
dataset$sequence = factor(as.character(dataset$sequence),c("real","mean_random"))

# Plot delta by tropism
diffexp = compare_means(deltaMFE~tropism, data=dataset1, method="wilcox.test") # group.by = c("cond","cell")
ggplot(dataset1, aes(x=tropism, y=deltaMFE, fill=tropism)) +  
  geom_violin(position=position_dodge(1),aes(fill = tropism)) +
  labs(title="KnotInFrame",y = "Delta MFE") + 
  scale_fill_manual(values=c("#0028ff80","#00d4ff80","#7dff7a80","#ffe60080","#ff470080","#80000080"), na.value="#dedede") + 
  stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
  theme_classic()

diffexp = compare_means(deltaMFE~sequence, data=dataset, method="wilcox.test", alternative="greater", paired = T, group.by = c("tropism"))
ggplot(dataset, aes(x=tropism, y=deltaMFE, fill=sequence)) + 
  geom_violin(position=position_dodge(1)) +
  labs(title="KnotInFrame",y = "Delta MFE") + 
  stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
  theme_classic()
