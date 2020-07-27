#### COMPUTE GC ####
spe = "Severe acute respiratory syndrome coronavirus 2"
spe = "Severe acute respiratory syndrome-related coronavirus"
output = data.frame(row.names = rownames(codus)[codus$Species %in% spe])
for (prot in rownames(output)){
  gc1 = 0
  gc2 = 0
  gc3 = 0
  for (c in colnames(codus)[15:ncol(codus)]){
    if (substr(c,1,1) %in% c("G","C")){
      gc1 = gc1 + codus[prot,c]
    }
    if (substr(c,2,2) %in% c("G","C")){
      gc2 = gc2 + codus[prot,c]
    }
    if (substr(c,3,3) %in% c("G","C")){
      gc3 = gc3 + codus[prot,c]
    }
  }
  total = sum(codus[prot,15:ncol(codus)])
  output[prot,"GC1"] = gc1*100/total
  output[prot,"GC2"] = gc2*100/total
  output[prot,"GC3"] = gc3*100/total
}


#### PLOT GC12 vs GC3 ####
library(ggplot2)

# Load data
codus = read.csv("data/refseq_humanvirus_CoCoPUT.tsv", row.names = 1, sep="\t")
metadata = read.csv("data/virus_list.tsv", sep="\t", row.names = 1)

# Extract variables
species = t(sapply(as.character(unique(codus$Species)), 
                   function(x) colMeans(codus[codus$Species %in% x,c("GC1.","GC2.","GC3."),drop=F],na.rm=T)))

# Define structure
dataset = data.frame(row.names = rownames(species))
dataset$tropism = factor(as.character(sapply(rownames(species), function(x) metadata[x,"tropism"])),
                         c("Immune cells","Epithelial cells","Neurons","Respiratory tract","Hepatocytes","Intestine"))
dataset$gc3 = as.numeric(species[,"GC3."])
dataset$gc12 = as.numeric(rowMeans(species[,c("GC1.","GC2.")]))

# Plot
ggplot(dataset, aes(x=gc3, y=gc12, color=tropism)) +
  geom_point(size=2,alpha=0.7,na.rm=F) + 
  scale_color_manual(values=c("#0028ff80","#00d4ff80","#7dff7a80","#ffe60080","#ff470080","#80000080"), na.value="#dedede") + 
  theme_classic()


