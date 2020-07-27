
### Upload codon usage data
codus = read.csv("data/refseq_humanvirus_CoCoPUT.tsv", sep="\t", row.names = 1)
species = as.character(unique(codus$Taxid))

# Upload viral orthology from VOGDB (Jun 21, 2019, vog94, NCBI Refseq release 94)
orthology = read.csv("data/vog.members.tsv", sep = "\t", row.names = 1)
orthframe = data.frame(vog = character(), functionalcategory = character(), prot_id = character(), species = character())
for (v in rownames(orthology)){
  extractdata = sapply(strsplit(as.character(orthology[v,"ProteinIDs"]), ",")[[1]], function(y) strsplit(y,"\\.")[[1]][1:2])
  idx=(extractdata[1,] %in% species)
  if (sum(idx)>1){
    extractdata = extractdata[,idx]
    tempframe = data.frame(vog = v, functionalcategory = orthology[v,"FunctionalCategory"], 
                           prot_id = extractdata[2,], tax_id = extractdata[1,],
                           row.names = 1:ncol(extractdata))
    orthframe = rbind(orthframe,tempframe)
  } else if (sum(idx)==1){
    extractdata = extractdata[,idx]
    tempframe = data.frame(vog = v, functionalcategory = orthology[v,"FunctionalCategory"], 
                           prot_id = extractdata[2], tax_id = extractdata[1],
                           row.names = 1)
    orthframe = rbind(orthframe,tempframe)
  }
}
orthframe = apply(orthframe,2,as.character)

# For each protein, find related functions and assign
to_add = t(sapply(rownames(codus), function(x) if (sum(orthframe[,"prot_id"] %in% x)>1){
  apply(orthframe[orthframe[,"prot_id"] %in% x,c("vog","functionalcategory")],2,paste,collapse =";")}
  else if (sum(orthframe[,"prot_id"] %in% x)==1){orthframe[orthframe[,"prot_id"] %in% x,c("vog","functionalcategory")]}
  else{c("NA","NA")}))

# Create GMT file for GSEA
annotations = read.csv("data/vog.annotations.tsv", sep = "\t", row.names = 1)
fileConn<-file("results/VOGDB_orthology_vog99.gmt")
# Write Vog by Vog
to_write=c()
for (v in rownames(orthology)){
  extractdata = sapply(strsplit(as.character(orthology[v,"ProteinIDs"]), ",")[[1]], function(y) strsplit(y,"\\.")[[1]][2])
  to_write = append(to_write, paste(c(sprintf("%s:%s",v,annotations[v,"ConsensusFunctionalDescription"]),
                   "http://vogdb.org/",
                   extractdata),collapse="\t"))
  
}
writeLines(to_write, fileConn)
close(fileConn)


