# Upload proteins
codus = read.csv("data/refseq_humanvirus_CoCoPUT.tsv", sep="\t", row.names = 1)

# Viral annotations
metadata = read.csv("data/virus_list.tsv",sep="\t",row.names = 1)
tropisms = unique(as.character(metadata$tropism_general)); tropisms = tropisms[!is.na(tropisms)]

# Find groups of proteins expanding various tropisms
vogdata = codus[!is.na(codus$VOG),]

# Add tropism
vogdata$tropism = sapply(as.character(vogdata$Species), function(x) metadata[x,"tropism"])

# Find VOG expanding several tropisms
allvogs = unique(unlist(sapply(as.character(vogdata$VOG),function(x) strsplit(x,";")[[1]])))
vogsummary = data.frame(row.names = allvogs)
for (t in tropisms){
  tempdata = vogdata[vogdata$tropism %in% t,]
  tempvog = unlist(sapply(as.character(tempdata$VOG),function(x) strsplit(x,";")[[1]]))
  vogsummary[,t] = sapply(allvogs,function(x) sum(tempvog %in% x))
}
vogsummary$tropism_count = apply(vogsummary,1,function(x) sum(x!=0))
