library(ggplot2)
library(ggpubr)

# Load tAIs
TAIs = read.csv("results/RtAI_viralinfection.csv", row.names = 1)

# Build dataset
dataset = c()
for (s in colnames(TAIs)[1:27]){
  cell = strsplit(s,"_|\\.")[[1]][1]
  virus = strsplit(s,"_|\\.")[[1]][2]
  cond = strsplit(s,"_|\\.")[[1]][3]
  # Keep proteins of the virus used
  idx = (TAIs$abbreviation %in% virus)
  temp_tai = TAIs[idx,s,drop=F]
  explabels = TAIs[idx,"expression"]
  isXr = grep("Xr",as.character(TAIs[idx,"vogdb_annotation"]))
  isXs = grep("Xs",as.character(TAIs[idx,"vogdb_annotation"]))
  # Keep data
  dataset_temp = data.frame(row.names = 1:nrow(temp_tai))
  dataset_temp$tai = as.numeric(unlist(temp_tai))
  dataset_temp$protein = rownames(temp_tai)
  dataset_temp$expression = explabels
  dataset_temp$vog = NA
  dataset_temp[isXr,"vog"] = "Xr"; dataset_temp[isXs,"vog"] = "Xs"
  dataset_temp$virus = virus
  dataset_temp$cond = cond
  dataset_temp$cell = cell
  dataset_temp = dataset_temp[!is.na(explabels),]
  dataset = rbind(dataset,dataset_temp)
}
# Order genes
dataset$expression <- factor(dataset$expression, levels = c("IE","E","L","polyprotein"))

# Compute p-values and correct FDR
diffexp1 = compare_means(tai~expression, data=dataset[!grepl("mock|uninf",dataset$cond),], group.by = c("cond","cell"), method="wilcox.test")
diffexp2 = compare_means(tai~cond, data=dataset[!grepl("mock|uninf",dataset$cond),], group.by = c("expression","cell"), method="wilcox.test",paired=T)
allcomp = merge(diffexp1,diffexp2,all=T)
allcomp$p.corrected = p.adjust(allcomp$p,method="fdr")
# Plot
ggplot(dataset[!grepl("mock|uninf",dataset$cond),], aes(x=cond, y=tai, fill=expression)) +
  facet_grid(. ~ cell+virus, scales = "free") +
  geom_violin(position=position_dodge(1)) +
  labs(title="tRNAs in infection",y = "tAI") + 
  stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
  theme_classic()
ggplot(dataset, aes(x=cond, y=tai, fill=expression)) +
  facet_grid(. ~ cell+virus, scales = "free") +
  geom_violin(position=position_dodge(1)) +
  labs(title="tRNAs in infection",y = "tAI") + 
  stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
  theme_classic()

## Compute delta wrt mock/uninfected
controls = dataset[grepl("mock|uninf",dataset$cond),]
infect = dataset[!grepl("mock|uninf",dataset$cond),]
# Compute average of controls
controls_time = sapply(controls$cond, function(x) substr(x,(nchar(x)-2),nchar(x)))
for (i in rownames(infect)){
  time = substr(infect[i,"cond"],(nchar(infect[i,"cond"])-2),nchar(infect[i,"cond"]))
  # Compute mock to subtract
  if (substr(time,3,3)=="h"){
    ctrl_tai = mean(controls[(controls$protein %in% infect[i,"protein"])&
                             (controls$virus %in% infect[i,"virus"])&
                             (controls$cell %in% infect[i,"cell"])&
                             (controls_time %in% time),"tai"])
    infect[i,"cond"] = time
  }else{
    ctrl_tai = mean(controls[(controls$protein %in% infect[i,"protein"])&
                             (controls$virus %in% infect[i,"virus"])&
                             (controls$cell %in% infect[i,"cell"]),"tai"])
  }
  # Make subtraction
  infect[i,"tai"] = infect[i,"tai"]-ctrl_tai
}

# Plot BOXPLOTS!!
diffexp = compare_means(tai~expression, data=infect, group.by = c("cond","cell"), method="wilcox.test")
ggplot(infect, aes(x=cond, y=tai, fill=expression)) +
  facet_grid(. ~ cell+virus, scales = "free") +
  geom_violin(position=position_dodge(1)) +
  labs(title="tRNAs in infection",y = "tAI") + 
  stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
  theme_classic()
