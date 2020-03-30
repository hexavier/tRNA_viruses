library(ggplot2)

# Load SDAs of coronaviruses
metadata = read.csv("data/SARScoronaviruses2019ALL.csv", row.names = 1)
showtypes = c("HNSC","LUSC","LUAD","CESC","UCEC","BLCA","PRAD","SKCM","BRCA","THCA","ESCA","STAD","COAD","READ","PAAD","LIHC","CHOL","KIRC","KIRP","KICH","PCPG","GBM","THYM")
SDA = read.csv("results/RTE_SARSCoV2s.csv", row.names = 1)

# Set strains and proteins to compare
otherstrains = c("NC_004718","MN996532") #MG772933
sars2 = "NC_045512"
otherprot = as.character(SDA[SDA$Accession %in% otherstrains,"description"])
sars2prot = as.character(SDA[SDA$Accession %in% sars2,"description"])
allprot = c(unique(otherprot[sapply(otherprot, function(x) sum(otherprot==x)==length(otherstrains))]),sars2prot)
sharedprot = allprot[duplicated(allprot)]
# Keep only selected strains and proteins
SDA = SDA[(SDA$Accession %in% c(otherstrains,sars2))&(SDA$description %in% sharedprot),]

# Create dataset structure
SDAdataset = c()
for (l in showtypes){
  for (s in otherstrains){
    # Dataset for SDA
    dataset_temp = data.frame(row.names = sharedprot)
    dataset_temp$deltaSDA = as.numeric(SDA[SDA$Accession %in% sars2,l]) - as.numeric(SDA[SDA$Accession %in% s,l])
    dataset_temp$protein = sharedprot
    dataset_temp$strain = as.character(metadata[s,"Collection.Date"])
    dataset_temp$catype = l
    SDAdataset = rbind(SDAdataset,dataset_temp)
  }
}
# Fix order of tissues along the descending respiratory tract
SDAdataset$catype <- factor(SDAdataset$catype,levels = showtypes)

# Plot
# Set cancer types to plot

ggplot(SDAdataset, aes(x=catype, y=deltaSDA, fill=strain)) +  
  geom_violin(position=position_dodge(1)) +
  stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
  theme_classic() +
  theme(axis.text.x=element_text(size = rel(1), angle=45, hjust=1, vjust =  1)) +
  labs(title="SARS-CoV-2 - other-CoV",y = "deltaSDA")

showprot = c("orf1ab polyprotein","spike glycoprotein","nucleocapsid phosphoprotein","membrane glycoprotein")
ggplot(SDAdataset[SDAdataset$protein %in% showprot,], aes(x=catype, y=deltaSDA, group=interaction(protein,strain))) +  
  geom_line(aes(color=protein,linetype=strain),size=1) +
  theme_classic() +
  theme(axis.text.x=element_text(size = rel(1), angle=45, hjust=1, vjust =  1)) +
  labs(title="SARS-CoV-2 - other-CoV",y = "deltaSDA")



