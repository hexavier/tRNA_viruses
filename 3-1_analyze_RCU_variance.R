calc_SS <- function(df) sum(as.matrix(dist(df,method = "euclidean")^2)) / (2 * nrow(df))

# Save output
species = read.csv("results/virus_RCUs_RelbyProt.csv", row.names = 1)
metadata = read.csv("data/virus_list.tsv", sep="\t", row.names = 1)

# Initate structure
scoretab = matrix(ncol= ncol(metadata), nrow = 3)
colnames(scoretab) = colnames(metadata); rownames(scoretab) = c("WB","Silhouette","Dunn")

# Calculate sum of squares
for (classifier in colnames(metadata)){
  data = species[(rownames(species) %in% rownames(metadata)[!is.na(metadata[,classifier])]),]
  totss = calc_SS(data)
  groups = unique(as.character(metadata[!is.na(metadata[,classifier]),classifier]))
  withinss = sum(sapply(groups, 
                        function(x) calc_SS(data[(rownames(data) %in% (rownames(metadata)[metadata[,classifier] %in% x])),])),na.rm = T)
  # Calculate how much of the total SS is between groups (separation) or within groups (cohesion) (ie. the variance that is explained by the grouping)
  # Ref: WB-index: A sum-of-squares based index for cluster validity
  score = (length(groups)*withinss/(totss - withinss))
  scoretab["WB",classifier] = score
}

# Silhouette index
library(cluster)
for (classifier in colnames(metadata)){
  notnameta = metadata[!is.na(metadata[,classifier]),]
  data = species[(rownames(species) %in% rownames(notnameta)),]
  labels = as.character(sapply(rownames(data), function(x) notnameta[x,classifier]))
  # Change labels to numeric id
  groups = 1:length(unique(as.character(notnameta[,classifier])))
  names(groups) = unique(as.character(notnameta[,classifier]))
  cluster_id = sapply(labels, function(x) groups[x])
  # Calculate score
  scorebyclust = silhouette(cluster_id,dist(data))
  score = summary(scorebyclust)$avg.width
  scoretab["Silhouette",classifier] = score
}

# Dunn index
library(fpc)
for (classifier in colnames(metadata)){
  notnameta = metadata[!is.na(metadata[,classifier]),]
  data = species[(rownames(species) %in% rownames(notnameta)),]
  labels = as.character(sapply(rownames(data), function(x) notnameta[x,classifier]))
  # Change labels to numeric id
  groups = 1:length(unique(as.character(notnameta[,classifier])))
  names(groups) = unique(as.character(notnameta[,classifier]))
  cluster_id = sapply(labels, function(x) groups[x])
  # Calculate score
  scorebyclust = cluster.stats(dist(data),cluster_id)
  score = scorebyclust$dunn
  scoretab["Dunn",classifier] = score
}

# Save
write.csv(scoretab,"results/RCU_variance.csv")
