library(gplots)

# Load data
features = read.csv("results/mean_randomforest_cells_proteins_bytropism.csv", row.names = 1)

# Heatmap
heatmap.2(t(features[,3:9]),symm=F, margins = c(11,5), na.color = "grey")