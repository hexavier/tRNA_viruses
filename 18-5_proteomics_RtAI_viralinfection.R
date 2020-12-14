library(tAI)
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

AAnormalize <- function(data,codons){
  # Compute relative data
  aa = sapply(codons,function(x) substr(x,1,nchar(x)-3))
  uniqueaa = unique(aa)
  outdata = numeric(length=length(data))
  for (n in uniqueaa){
    idx = (aa %in% n)
    #total = sum(data[idx],na.rm=T)
    total = max(data[idx],na.rm=T)
    outdata[idx] = data[idx]/total
    if (total %in% 0){
      outdata[idx] = 1.0/sum(idx)
    }
  }
  return(outdata)
}

## Load trna and weighted CU
# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)

# tRNAs
trna = read.csv("data/viralinfection_tRNAs.csv",row.names = 1)
trna = extract_cod(transformdata(trna,"sqrt"),codons$ANTICODON)

# Codon usage
codus = read.csv("data/proteomics/refseq_proteomics_codoncount.tsv",sep="\t", row.names = 1)
# Keep only columns with codon info
colnames(codus) = sapply(colnames(codus),function(x) paste(codons[x,"AA"],x,sep=""))

## Calculate tAI for genomic CU
codon = extract_cod(transformdata(t(codus),""),rownames(codons)[!(codons$AA %in% c("Stop","Met"))])

# Calculate tAI
initial_s = c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5)
anticodon = apply(trna,2, get.ws, s=initial_s, sking=0)
anticodon = apply(anticodon,2,AAnormalize,paste0(codons[rownames(codon),"AA"],rownames(codon)))# normalize by AA
rownames(anticodon)= rownames(codon)

TAIs = data.frame(matrix(ncol = ncol(anticodon), nrow = ncol(codon)),row.names = colnames(codon)); colnames(TAIs) = colnames(anticodon)
# Calculate tAI
for (sample in colnames(anticodon)){
  # Calculate relative adaptiveness values (ws)
  sample.ws = as.numeric(anticodon[,sample])
  # Calculate tAI for all CUs
  sample.tai <- get.tai(t(codon), sample.ws)
  TAIs[,sample] = sample.tai
}

write.csv(TAIs,"results/RtAI_proteomics_viralinfection.csv")
