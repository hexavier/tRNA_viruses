library(RTCGAToolbox)

get_expression <- function(abbr){
  if (sum(grep(";",abbr))>0){
    names = unlist(strsplit(abbr,";"))
    for (n in names){
      dwload = getFirehoseData(dataset=n, RNASeq2GeneNorm=TRUE, destdir="/home/xhernandez/Downloads/TCGA-mRNAseq")
      if (!exists("output")){
        output = getData(dwload,"RNASeq2GeneNorm")[[1]]@DataMatrix
      }else{
        toadd = getData(dwload,"RNASeq2GeneNorm")[[1]]@DataMatrix
        output = cbind(output, toadd)
      }
    }
  }else{
    dwload = getFirehoseData(dataset=abbr, RNASeq2GeneNorm=TRUE, destdir="/home/xhernandez/Downloads/TCGA-mRNAseq")
    output = getData(dwload,"RNASeq2GeneNorm")[[1]]@DataMatrix
  }
  return(output)
}

## Analyze subclinical subsets for certain types
cancer_types = c("COAD","READ","GBM","BRCA","STAD","KIRP","KICH", "KIRC", "PRAD", "THCA","BLCA", "LIHC", "PAAD", "UCEC", "PCPG", "SKCM", "CHOL", "ESCA", "THYM", "CESC","LUAD","LUSC","HNSC")
genes = c("ACE2","TMPRSS2","BSG","DPP4","ST6GAL1","ST3GAL4")

dataset = c()

for (type in cancer_types){
  # Upload gene expression data
  rnaInfo = get_expression(type)
  
  ## Analyze expression status of genes
  expr_samples = substr(colnames(rnaInfo),1,15)
  mapped_normal = (regexpr("TCGA-[A-Z,0-9]{2}-[A-Z,0-9]{4}-11",expr_samples))==1
  expr_samples = expr_samples[mapped_normal]
  rnaInfo = rnaInfo[,mapped_normal, drop=FALSE]
  
  # Select genes
  genes_temp = genes[(genes %in% rownames(rnaInfo))]
  
  # Save data
  dataset_temp = data.frame(row.names = expr_samples)
  dataset_temp[,genes_temp] = t(rnaInfo[genes_temp,])
  
  if (type=="HNSC"){
    # Subclinical classification
    clinic = t(read.csv(sprintf("data/%s.clin.merged.txt",type), row.names = 1))
    patients=gsub("\\.","-",toupper(rownames(clinic)))
    classification = as.character(sapply(rownames(dataset_temp),function(x) as.character(clinic[grep(substr(x,1,12),patients),"patient.anatomic_neoplasm_subdivision"])))
    dataset_temp$catype = sprintf("%s-%s",rep(type,length(classification)),classification)
  }else{
    dataset_temp$catype = type
  }
  dataset = rbind(dataset, dataset_temp)
}

outtable = data.frame(row.names=unique(dataset$catype))
outtable[,genes] = NA

for (t in rownames(outtable)){
  outtable[t,genes] = colMeans(dataset[dataset$catype==t,genes],na.rm=T)
}

write.csv(outtable, "results/receptors_expression_coronavirus.csv")
