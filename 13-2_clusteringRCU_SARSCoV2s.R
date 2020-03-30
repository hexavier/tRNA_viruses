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

## Load trna and weighted CU
# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)

# Genomic codon usage
codus = read.csv("data/CU_SARScoronaviruses2019ALL.tsv",sep="\t", row.names = 1)
# Keep only columns with codon info
codus_clean = data.frame(sapply(rownames(codus),function(x) as.numeric(codus[x,3:ncol(codus)])), row.names = colnames(codus)[3:ncol(codus)])
rownames(codus_clean) = sapply(rownames(codus_clean),function(x) paste(codons[x,"AA"],x,sep=""))

## Calculate tAI for genomic CU
codon = extract_cod(transformdata(codus_clean,"rel"),rownames(codons)[!(codons$AA %in% c("Stop","Met"))])

## Anlyze differences in codus
metadata = read.csv("data/SARScoronaviruses2019ALL.csv", row.names = 1)
mainprot = c("orf1ab polyprotein","spike glycoprotein","nucleocapsid phosphoprotein","membrane glycoprotein")
# Create structure containing codon usage of protein by protein
protcodus = list()
for (p in mainprot){
  codustemp = as.matrix(t(codon[,codus$description %in% p]))
  date = sapply(rownames(codustemp), function(x) as.character(metadata[as.character(codus[x,"Accession"]),"Collection.Date"]))
  place = sapply(rownames(codustemp), function(x) as.character(metadata[as.character(codus[x,"Accession"]),"Locality"]))
  rownames(codustemp) = sprintf("%s-%s",date,place)
  protcodus[[p]] = codustemp
}

plot(hclust(dist(protcodus[[1]])))
