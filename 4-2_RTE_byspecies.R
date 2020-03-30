library(gplots)

# load data
TAIs = read.csv("results/RTE_tissueMeans.csv",row.names = 1)

# Compute average of each pathway
RTEspecies = t(sapply(unique(as.character(TAIs$Species)), 
                                      function(x) colMeans(TAIs[(TAIs$Species %in% x),colnames(TAIs)[1:23]],na.rm=T)))

# Save output
write.csv(RTEspecies,"results/species_RTEs.csv")

### MERGE VIRUSES ###
classificator = "tropism"
metadata = read.csv("data/virus_list.tsv",sep="\t",row.names = 1)
RTEmerged = sapply(as.character(unique(metadata[!is.na(metadata[,classificator]),classificator])), 
                    function(x) if (sum(rownames(RTEspecies) %in% rownames(metadata)[metadata[,classificator] %in% x])>1){
                      colMeans(RTEspecies[(rownames(RTEspecies) %in% rownames(metadata)[metadata[,classificator] %in% x]),],na.rm = T)}else{
                        RTEspecies[(rownames(RTEspecies) %in% rownames(metadata)[metadata[,classificator] %in% x]),]
                      })

## PLOT ##
# Heatmap
heatmap.2(RTEmerged,symm=F, margins = c(10,5))
