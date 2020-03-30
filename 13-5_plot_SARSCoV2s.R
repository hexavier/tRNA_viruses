library(ggplot2)

# Load SDAs of coronaviruses
metadata = read.csv("data/SARScoronaviruses2019ALL.csv", row.names = 1)
showtypes = c("HNSC-oral tongue","HNSC-floor of mouth","HNSC-base of tongue","HNSC-oral cavity","HNSC-larynx","LUSC","LUAD")
SDA = read.csv("results/RTE_SARSCoV2s.csv", row.names = 1)
subclinicSDA = read.csv("results/subclinicRTE_SARSCoV2s.csv", row.names = 1)
allSDA = cbind(SDA,subclinicSDA)

# Create dataset structure
SDAdataset = c()
for (l in showtypes){
  # Dataset for SDA
  type = gsub(" |-",".",l)
  dataset_temp = data.frame(row.names = rownames(allSDA))
  dataset_temp$SDA = as.numeric(allSDA[,type])
  dataset_temp$protein = as.character(allSDA$description)
  dataset_temp$date = as.character(metadata[allSDA$Accession,"Collection.Date"])
  dataset_temp$place = as.character(metadata[allSDA$Accession,"Locality"])
  dataset_temp$catype = l
  SDAdataset = rbind(SDAdataset,dataset_temp)
}

# Plot
# Set cancer types to plot

ggplot(SDAdataset, aes(x=catype, y=SDA, fill=date)) +  
  geom_violin(position=position_dodge(1)) +
  stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
  theme_classic() +
  theme(axis.text.x=element_text(size = rel(1), angle=45, hjust=1, vjust =  1)) +
  theme(legend.position = "none") +
  labs(title="SARS-CoV-2",y = "SDA")

# Create dataset structure with the mean over proteins
dataset = c()
sequences = as.character(unique(allSDA$Accession))
for (l in showtypes){
  # Dataset for SDA
  type = gsub(" |-",".",l)
  dataset_temp = data.frame(row.names = sequences)
  dataset_temp$SDA = as.numeric(sapply(sequences, function(x) mean(allSDA[allSDA$Accession %in% x,type])))
  dataset_temp$date = sprintf("%s-%s",as.character(metadata[sequences,"Collection.Date"]),as.character(metadata[sequences,"Locality"]))
  dataset_temp$catype = l
  dataset = rbind(dataset,dataset_temp)
}

ggplot(dataset, aes(x=date, y=SDA, group=catype, color=catype)) +  
  geom_line() +
  geom_point() +
  theme_classic() +
  theme(axis.text.x=element_text(size = rel(1), angle=45, hjust=1, vjust =  1)) +
  labs(title="SARS-CoV-2",y = "SDA")
