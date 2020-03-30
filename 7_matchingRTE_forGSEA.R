
# load data
TAIs = read.csv("results/RTE_tissueMeans.csv",row.names = 1)
metadata = read.csv("data/virus_list.tsv",sep="\t",row.names = 1)

# Define matching
mapping = c(COAD="Intestine",READ="Intestine",GBM="Neurons",LIHC="Hepatocytes",CHOL="Hepatocytes",
            LUAD="Respiratory tract",LUSC="Respiratory tract",THYM="Immune cells",UCEC="Epithelial cells",
            CESC="Epithelial cells",HNSC="Epithelial cells",SKCM="Epithelial cells")

# Keep RTE of matching proteins
output = data.frame(row.names = rownames(TAIs))
output$RTE = apply(TAIs,1,function(x) mean(as.numeric(x[(colnames(TAIs) %in% names(mapping)[(mapping %in% metadata[x["Species"],"tropism"])])])))

# Write output
write.table(output,"results/matchingtissueRTEs.rnk", sep= "\t")

