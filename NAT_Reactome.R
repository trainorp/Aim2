############ Import real data ############
load(paste(stem,"BearOmics2/data/transcripts.RData",sep="/"))
load(paste(stem,"BearOmics2/data/genesKey.RData",sep="/"))
load(paste(stem,"BearOmics2/data/metabs.RData",sep="/"))
load(paste(stem,"BearOmics2/data/metabKey.RData",sep="/"))

# Check for metabolites in edge set
metabKey$ChEBI<-""
for(i in 1:nrow(metabKey))
{
  if(metabKey$hmdb[i]!="")
  {
    look1<-ChEBINames$COMPOUND_ID[match(metabKey$hmdb[i],ChEBINames$ACCESSION_NUMBER)]
    metabKey$ChEBI[i]<-ifelse(!is.na(look1),look1,"")
  }
  if(metabKey$ChEBI[i]=="" & metabKey$kegg[i]!="")
  {
    look1<-ChEBINames$COMPOUND_ID[match(metabKey$kegg[i],ChEBINames$ACCESSION_NUMBER)]
    metabKey$ChEBI[i]<-ifelse(!is.na(look1),look1,"")
  }
  if(metabKey$ChEBI[i]=="" & metabKey$cas[i]!="")
  {
    look1<-ChEBINames$COMPOUND_ID[match(metabKey$cas[i],ChEBINames$ACCESSION_NUMBER)]
    metabKey$ChEBI[i]<-ifelse(!is.na(look1),look1,"")
  }
}

# Check for genes in edge set
genesKey$found<-0
for(i in 1:nrow(genesKey))
{
  if(genesKey$uniprotswissprot[i]!="")
  {
    look1<-match(genesKey$uniprotswissprot[i],edges$entityAid) |
      match(genesKey$uniprotswissprot[i],edges$entityBid)
    genesKey$found[i]<-ifelse(!is.na(look1),1,0)
  }
  if(genesKey$found[i]==0 & genesKey$uniprotsptrembl[i]!="")
  {
    look1<-match(genesKey$uniprotsptrembl[i],edges$entityAid) |
      match(genesKey$uniprotsptrembl[i],edges$entityBid)
    genesKey$found[i]<-ifelse(!is.na(look1),1,0)
  }
  if(i%%1000==0) print(i)
}
genesKey<-genesKey %>% left_join(genesKey %>% group_by(gene_id) %>% summarise(foundInd=ifelse(sum(found)>0,1,0)))

# Collapse gene key for unique ENSG:
genesKey$found<-NULL
colFun<-function(x)
{
  x<-unique(x[x!=""])
  return(paste(x,collapse=";"))
}
colFun2<-function(x)
{
  return(as.data.frame(t(apply(x,2,colFun))))
}
byCol<-by(genesKey[,names(genesKey)!="gene_id"],genesKey$gene_id,colFun2)
byCol<-do.call("rbind",byCol)

genesNotFound<-genesKey %>% dplyr::select(gene_id,foundInd,gene_biotype) %>% unique() 
table(genesNotFound$gene_biotype,genesNotFound$foundInd)
genesNotFound<-genesNotFound %>% filter(gene_biotype=="protein_coding")

############ Abundance & Expression processing ############
library(KernSmooth)

fun1<-function(x)
{
  return(x[x>0])
}

hist(log2(fun1(transcripts[1,])))

bkde1<-bkde(log2(fun1(transcripts[1,])))
plot(bkde1$x,bkde1$y,type="l",ylim=c(0,.10))
for(i in 2:nrow(transcripts))
{
  bkde1<-bkde(log2(fun1(transcripts[i,])))
  points(bkde1$x,bkde1$y,type="l")
}

############ Indicator matrix ############
c(metabKey$id,colnames(transcripts))
