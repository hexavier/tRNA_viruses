library(ggplot2)
library(ggpubr)

# Load proteomics
proteomics = read.delim("data/proteomics/proteomics.tsv", sep="\t")

# Build dataset
dataset = c()
for (s in colnames(proteomics)[8:74]){
  virus = strsplit(s,"_|\\.")[[1]][1]
  cond = strsplit(s,"_|\\.")[[1]][2]
  time = strsplit(s,"_|\\.")[[1]][3]; time = as.numeric(substr(time,1,nchar(time)-1))
  # Keep proteins of the virus used
  idx = (proteomics$Species %in% virus)
  temp_prot = proteomics[idx,s,drop=F]
  explabels = proteomics[idx,"expression"]
  proteins = proteomics[idx,"RefSeq"]
  # Keep data
  dataset_temp = data.frame(row.names = 1:nrow(temp_prot))
  dataset_temp$proteomics = as.numeric(unlist(temp_prot))
  dataset_temp$protein = proteins
  dataset_temp$expression = explabels
  dataset_temp$virus = virus
  dataset_temp$cond = cond
  dataset_temp$time = time
  dataset = rbind(dataset,dataset_temp)
}
# Get averages
avg_dataset = aggregate(dataset$proteomics,
                        by=list(protein=dataset$protein,virus=dataset$virus,time=dataset$time,cond=dataset$cond),
                        data=dataset,FUN=mean,na.action = na.omit)
avg_dataset$vog = sapply(as.character(avg_dataset$protein), function(x) as.character(unique(proteomics[proteomics$RefSeq %in% x,"vogdb_annotation"])))
avg_dataset$expression = sapply(as.character(avg_dataset$protein), function(x) as.character(unique(proteomics[proteomics$RefSeq %in% x,"expression"])))

# Order genes
avg_dataset$expression <- factor(avg_dataset$expression, levels = c("IE","E","L"))

# Plot
ggplot(avg_dataset[avg_dataset$cond=="inf",], aes(x=time, y=x)) +
  facet_wrap(~virus, scales = "free") +
  labs(title="Proteomics",y = "Protein abundance") + 
  stat_summary(aes(color = expression),fun=median, geom="line", size=2) + 
  stat_summary(aes(fill = expression),fun.data = median_hilow, geom = "ribbon",alpha = 0.3, fun.args=c(conf.int=0.2)) +
  theme_classic()
ggplot(avg_dataset[avg_dataset$cond=="inf",], aes(x=time, y=x, group=protein, color=expression)) +
  facet_wrap(~virus, scales = "free") +
  geom_line(size=1,alpha=0.5) +
  geom_point() +
  labs(title="Proteomics",y = "Protein abundance") + 
  theme_classic()

# Wrt max
rel_dataset = c()
for (p in unique(avg_dataset$protein)){
  tempdata = avg_dataset[avg_dataset$protein %in% p,]
  norm_factor = max(tempdata$x)
  if (length(norm_factor)>0){
    tempdata$x = tempdata$x/norm_factor
    rel_dataset = rbind(rel_dataset,tempdata)
  }
}

# Plot
ggplot(rel_dataset[rel_dataset$cond=="inf",], aes(x=time, y=x, color=expression)) +
  facet_wrap(~virus, scales = "free") +
  labs(title="Proteomics",y = "Protein abundance") + 
  stat_summary(fun=median, geom="line", size=2) + 
  stat_summary(fun.data = median_hilow, geom = "errorbar",width = 0.5, fun.args=c(conf.int=0.5)) +
  theme_classic()
