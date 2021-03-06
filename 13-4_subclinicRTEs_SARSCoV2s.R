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

AAnormalize <- function(data,codons){
  # Compute relative data
  aa = sapply(codons,function(x) substr(x,1,nchar(x)-3))
  uniqueaa = unique(aa)
  outdata = numeric(length=length(data))
  for (n in uniqueaa){
    idx = (aa %in% n)
    total = max(data[idx],na.rm=T)
    outdata[idx] = data[idx]/total
    if (total %in% 0){
      outdata[idx] = 1.0/sum(idx)
    }
  }
  return(outdata)
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

## Load trna and weighted CU
# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)

# tRNAs
trna = t(read.csv("data/matched_AAcTE.csv",row.names = 1))

# Keep only LUAD, LUSC and HNSC
tissues = sapply(colnames(trna), function(x) strsplit(x, "\\.")[[1]][1])
trna = trna[,tissues %in% c("LUAD","LUSC","HNSC")]

# Add subclinical data
tissues = sapply(colnames(trna), function(x) strsplit(x, "\\.")[[1]][1])
for (type in c("LUAD","LUSC","HNSC")){
  # Get clinical data
  clinic = t(read.csv(sprintf("data/%s.clin.merged.txt",type), row.names = 1))
  # Map patients
  SDA_patients = sapply(colnames(trna)[tissues %in% type], function(x) substr(x,6,17))
  clin_patients=gsub("\\.","-",toupper(rownames(clinic)))
  classification = as.character(sapply(SDA_patients,function(x) as.character(clinic[grep(substr(x,1,12),clin_patients),"patient.anatomic_neoplasm_subdivision"])))
  colnames(trna)[tissues %in% type] = sprintf("%s-%s_%s",rep(type,length(classification)),classification,SDA_patients)
}

# Compute mean of tissue
tissues = sapply(colnames(trna), function(x) strsplit(x, "_")[[1]][1])
trna_mean = data.frame(row.names = rownames(trna))
for (t in unique(tissues)){
  trna_mean[,t] = rowMeans(trna[,tissues %in% t], na.rm = T)
}
anticodon = trna_mean

# Genomic codon usage
codus = read.csv("data/CU_SARScoronaviruses2019ALL.tsv",sep="\t", row.names = 1)
# Keep only columns with codon info
codus_clean = data.frame(sapply(rownames(codus),function(x) as.numeric(codus[x,3:ncol(codus)])), row.names = colnames(codus)[3:ncol(codus)])
rownames(codus_clean) = sapply(rownames(codus_clean),function(x) paste(codons[x,"AA"],x,sep=""))

## Calculate tAI for genomic CU
codon = extract_cod(transformdata(codus_clean,""),rownames(codons)[!(codons$AA %in% c("Stop","Met"))])

TAIs = data.frame(matrix(ncol = ncol(anticodon), nrow = ncol(codon)),row.names = colnames(codon)); colnames(TAIs) = colnames(anticodon)
initial_s = c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5)
# Calculate tAI
for (sample in colnames(anticodon)){
  # Calculate relative adaptiveness values (ws)
  sample.ws = as.numeric(anticodon[,sample])
  # Calculate tAI for all CUs
  sample.tai <- get.tai(t(codon), sample.ws)
  TAIs[,sample] = sample.tai
}
TAIs[,c("description","Accession")] = codus[,c("description","Accession")]
write.csv(TAIs,"results/subclinicRTE_SARSCoV2s.csv")
