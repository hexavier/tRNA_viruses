# Tropism info
#ICTV Master Species List
hostmap = read.csv("data/ICTV Master Species List 2018b.v2.tsv", sep="\t", row.names="Species")

# Human virus
#https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239
viruses = read.csv("data/virus_list.tsv",sep="\t")

# Map virus to virus genus
viruses$family = sapply(as.character(viruses$Species), function(x) if (any(rownames(hostmap) %in% x)){as.character(hostmap[x,"Family"])}else{NA})
viruses$genus = sapply(as.character(viruses$Species), function(x) if (any(rownames(hostmap) %in% x)){as.character(hostmap[x,"Genus"])}else{NA})
viruses$type = sapply(as.character(viruses$Species), function(x) if (any(rownames(hostmap) %in% x)){as.character(hostmap[x,"Genome.Composition"])}else{NA})

# Aiewsakun_2018
seqclass = read.csv("data/Aiewsakun_2018.tsv", sep="\t")
viruses$Aiewsakun2018 = sapply(as.character(viruses$Species), function(x) if (any(seqclass$virus %in% x)){as.character(seqclass[(seqclass$virus %in% x),"classification"])}else{NA})

# Save output (but need to manually curate output)
write.table(viruses,"data/virus_list_toadd.tsv",sep = "\t")
