library(gplots)

# load data
TAIs = read.csv("results/RTE_tissueMeans.csv",row.names = 1)

# Compute average of each pathway
RTEspecies = t(sapply(unique(as.character(TAIs$Species)), 
                                      function(x) colMeans(TAIs[(TAIs$Species %in% x),colnames(TAIs)[1:23]],na.rm=T)))

### Analyze delta RTE ###
Xs_prot = TAIs[grep("Xs",TAIs$annotation),]
Xs = t(sapply(unique(as.character(Xs_prot$Species)), 
              function(x) colMeans(Xs_prot[(Xs_prot$Species %in% x),colnames(Xs_prot)[1:23]],na.rm=T)))
Xr_prot = TAIs[grep("Xr",TAIs$annotation),]
Xr = t(sapply(unique(as.character(Xr_prot$Species)), 
              function(x) colMeans(Xr_prot[(Xr_prot$Species %in% x),colnames(Xr_prot)[1:23]],na.rm=T)))

# Map species in Xr and Xs
idx = as.character(rownames(Xs)[rownames(Xs) %in% rownames(Xr)])
# Compute delta
deltaRTE = Xr[idx,] - Xs[idx,]
classificator = "tropism"
metadata = read.csv("data/virus_list.tsv",sep="\t",row.names = 1)
RTEmerged = sapply(as.character(unique(metadata[!is.na(metadata[,classificator]),classificator])), 
                   function(x) if (sum(rownames(deltaRTE) %in% rownames(metadata)[metadata[,classificator] %in% x])>1){
                     colMeans(deltaRTE[(rownames(deltaRTE) %in% rownames(metadata)[metadata[,classificator] %in% x]),],na.rm = T)}else{
                       deltaRTE[(rownames(deltaRTE) %in% rownames(metadata)[metadata[,classificator] %in% x]),]
                     })
# Heatmap
heatmap.2(RTEmerged,symm=F, margins = c(10,5))