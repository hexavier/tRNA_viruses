library(gplots)

# Load data
features = read.csv("results/mean_randomforest_proteins_bytropism.csv", row.names = 1)

# Heatmap
heatmap.2(t(features[,3:25]),symm=F, col=bluered(255), margins = c(11,5), na.color = "grey")
