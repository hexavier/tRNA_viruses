library(ggplot2)

extract_cod <- function (trnas, anticod){
  output = data.frame(row.names = anticod)
  trnas_acod = sapply(rownames(trnas), function(x) substr(x,nchar(x)-2,nchar(x)))
  for (s in colnames(trnas)){
    output[,s] = sapply(anticod, function(x) if(any(trnas_acod==x)){mean(trnas[trnas_acod==x,s])}else{0})
  }
  return(output)
}

transformdata <- function(data,transf){
  aa_idx = regexpr("i?[A-Z][a-z]{2}[A-Z]{3}",rownames(data))==1
  data = data[aa_idx,]
  if (transf=="log"){
    outdata = sapply(data,log)
    # Remove inf values
    outdata[outdata==-Inf] = NaN
    rownames(outdata)=rownames(data)
  }else if (transf=="arcsinh"){
    outdata = sapply(data,asinh)
    rownames(outdata)=rownames(data)
  }else if (transf=="sqrt"){
    outdata = sapply(data,sqrt)
    rownames(outdata)=rownames(data)
  }else if (transf=="rel"){
    # Compute relative data
    outdata = data.frame(matrix(ncol = ncol(data), nrow = nrow(data)),row.names = rownames(data))
    colnames(outdata)= colnames(data)
    aa = sapply(rownames(outdata),function(x) substr(x,1,nchar(x)-3))
    uniqueaa = unique(aa)
    for (n in uniqueaa){
      idx = (aa %in% n)
      idx_data = matrix(as.matrix(data[idx,]), ncol = ncol(data), nrow = sum(idx))
      total = colSums(idx_data)
      outdata[idx,] = t(apply(idx_data,1,function(x) x/total))
      iszero = (total %in% 0)
      if (any(iszero)){
        outdata[idx,iszero] = 1.0/sum(idx)
      }
    }
  }else{
    outdata=data
  }
  return(outdata)
}

# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)

## CALCULATE CU OF STRUCTURAL PROTEINS
# Human CU
codus = read.csv("data/refseq_humanvirus_CoCoPUT.tsv",sep="\t")
is_str = grep("Xs",codus$annotation)
codus=codus[is_str,]
codus_clean = t(codus[,14:ncol(codus)])

# Compute the RCU
rownames(codus_clean) = sapply(rownames(codus_clean),function(x) paste(codons[x,"AA"],x,sep=""))
codon = transformdata(codus_clean,"rel")

# Compute average of each pathway
species = data.frame(row.names=as.character(unique(codus$Species)))
species[,rownames(codon)] = t(sapply(rownames(species), 
                                    function(x) if (sum(codus$Species %in% x)>1){rowMeans(codon[,codus$Species %in% x],na.rm=T)}
                                    else if (sum(codus$Species %in% x)==1){codon[,codus$Species %in% x]}))
# Save output
write.csv(species,"results/virus_Xs_RCUs.csv")
