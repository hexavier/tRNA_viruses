library(tAI)
library(ggplot2)
library(ggpubr)

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
anticodon=trna

# Detect tissues
tissues = sapply(colnames(trna), function(x) strsplit(x, "\\.")[[1]][1])

# Genomic codon usage
codus = read.csv("data/refseq_humanvirus_CoCoPUT.tsv",sep="\t", row.names = 1)
# Keep only columns with codon info
codus_clean = data.frame(sapply(rownames(codus),function(x) as.numeric(codus[x,15:ncol(codus)])), row.names = colnames(codus)[15:ncol(codus)])
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
TAIs[,c("annotation","Accession","Species")] = codus[,c("annotation","Accession","Species")]

### Analyze delta RTE ###
Xs_prot = TAIs[grep("Xs",TAIs$annotation),]
Xs = t(sapply(unique(as.character(Xs_prot$Species)), 
              function(x) colMeans(Xs_prot[(Xs_prot$Species %in% x),colnames(Xs_prot)[1:ncol(anticodon)]],na.rm=T)))
Xr_prot = TAIs[grep("Xr",TAIs$annotation),]
Xr = t(sapply(unique(as.character(Xr_prot$Species)), 
              function(x) colMeans(Xr_prot[(Xr_prot$Species %in% x),colnames(Xr_prot)[1:ncol(anticodon)]],na.rm=T)))

# Map species in Xr and Xs
idx = as.character(rownames(Xs)[rownames(Xs) %in% rownames(Xr)])
# Compute delta
Xr = Xr[idx,]
Xs = Xs[idx,]
metadata = read.csv("data/virus_list.tsv",sep="\t",row.names = 1)

# Analyze matching
mapping = c(COAD="Intestine",READ="Intestine",GBM="Neurons",LIHC="Hepatocytes",CHOL="Hepatocytes",
            LUAD="Respiratory tract",LUSC="Respiratory tract",THYM="Immune cells",UCEC="Epithelial cells",
            CESC="Epithelial cells",HNSC="Epithelial cells",SKCM="Epithelial cells")

# Reestructure data in rows
dataset = c()

for (v in rownames(Xr)){
  if (!is.na(metadata[v,"tropism"])){
    colmatch = (tissues %in% names(mapping)[(mapping %in% metadata[v,"tropism"])])
    dataset_temp = data.frame(row.names=1:(sum(colmatch)*2))
    dataset_temp$virus = v
    dataset_temp$sample = rep(colnames(Xr)[colmatch],2)
    dataset_temp$RTE = c(as.numeric(Xr[v,colmatch]),as.numeric(Xs[v,colmatch]))
    dataset_temp$subset = c(rep("Xr",sum(colmatch)),rep("Xs",sum(colmatch)))
    dataset_temp$tropism = metadata[v,"tropism"]
    dataset = rbind(dataset, dataset_temp)
  }
}

# Plot delta by tropism
ggplot(dataset, aes(x=tropism, y=RTE, fill=subset)) +  
  geom_jitter(na.rm=T , size=1, alpha = 0.15 ,aes(color = subset), 
              position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(lwd=1, notch=F, outlier.shape = NA , na.rm=T,alpha=0.2, 
               position = position_dodge(width=0.9)) +
  scale_fill_manual(values=c("#1f439aff","#cc0000ff")) +
  scale_color_manual(values=c("#1f439aff","#cc0000ff")) +
  labs(title="Matching RTE",y = "RTE") + 
  stat_compare_means(aes(label = ..p.format..),method="wilcox.test",paired = T) +
  theme_classic()

# Plot density
plot_hist <- function(dataset, title)
{
  plot(density(dataset[(dataset$subset=="Xs"),"RTE"],na.rm=T),lwd=2,main=title,xlab="RTE",col = rgb(1,0,0,1))
  lines(density(dataset[(dataset$subset=="Xr"),"RTE"],na.rm=T),lwd=2,col = rgb(0,0,1,1))
  col = append(col,rgb(1,0,0,0.5))
  pval = wilcox.test(dataset[(dataset$subset=="Xs"),"RTE"],dataset[(dataset$subset=="Xr"),"RTE"], paired = T)[[3]]
  mtext(sprintf("pv = %.3f",pval),line=-10,adj=0.95)
  
  legend("topright", col=c(rgb(1,0,0,1),rgb(0,0,1,1)),c("Xs","Xr"), lwd=10)
}

plot_hist(dataset, "Gene Subsets")
