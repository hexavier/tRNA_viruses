library(ggplot2)
library(ggpubr)
library(scales)

# Load data
mfe = read.csv("results/MFE_humanvirus.tsv", sep="\t", row.names = 1)
metadata = read.csv("data/virus_list.tsv",sep="\t",row.names = 1)

# Reestructure data in rows
dataset1 = data.frame(row.names = rownames(mfe))
dataset1$MFE = mfe$WT
dataset1$sequence = "real"
dataset1$tropism = factor(as.character(metadata[mfe$Species,"tropism"]),
                         levels = c("Immune cells","Epithelial cells","Neurons","Respiratory tract","Hepatocytes","Intestine"))
# Random MFE
dataset2 = data.frame(row.names = rownames(mfe))
dataset2$MFE = rowMeans(mfe[,c("rand1","rand2","rand3","rand4","rand5","rand6","rand7","rand8","rand9","rand10")],na.rm=T)
dataset2$sequence = "mean_random"
dataset2$tropism = factor(as.character(metadata[mfe$Species,"tropism"]),
                          levels = c("Immune cells","Epithelial cells","Neurons","Respiratory tract","Hepatocytes","Intestine"))
# Create full dataset
isna = (!is.na(dataset1$MFE))&(!is.na(dataset2$MFE))
dataset = rbind(dataset1[isna,],dataset2[isna,])
dataset$sequence = factor(as.character(dataset$sequence),c("real","mean_random"))

# Plot delta by tropism
diffexp = compare_means(MFE~sequence, data=dataset, method="wilcox.test",alternative="greater" ,paired = T, group.by = c("tropism"))
ggplot(dataset, aes(x=tropism, y=MFE, fill=sequence)) + 
  geom_violin(position=position_dodge(1)) +
  labs(title="MFE",y = "MFE") + 
  stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
  scale_y_continuous(breaks=c(0,-10,-50,-100,-200,-300,-500,-1000,-2000,-3000,-5000,-8000),trans=trans_new("arcsinh",asinh,sinh)) +
  theme_classic()

# Plot delta MFE
deltadataset = dataset1
deltadataset$MFE = dataset1$MFE - dataset2$MFE
ggplot(deltadataset, aes(x=tropism, y=MFE, fill=tropism)) +  
  geom_violin(position=position_dodge(1),aes(fill = tropism)) +
  labs(title="MFE",y = "Delta MFE") + 
  scale_fill_manual(values=c("#0028ff80","#00d4ff80","#7dff7a80","#ffe60080","#ff470080","#80000080"), na.value="#dedede") + 
  stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
  theme_classic()
