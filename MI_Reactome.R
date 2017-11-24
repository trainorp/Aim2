library(entropy)
library(glasso)
library(BayesianGLasso)
library(dplyr)
library(tidyr)

load(file="~/gdrive/Dissertation/Aim2/data/Reactome/reactomeModel.RData")

#Metab data
setwd("~/gdrive/AthroMetab")
key<-read.csv("Data/metabolite_key2.csv")
key[key$id=="M801",]$super<-"Xenobiotics"
key$biochemical<-as.character(key$biochemical)
key$Unknown.Name<-as.character(key$Unknown.Name)
key[grepl("X - ",key$biochemical),]$biochemical<-key[grepl("X - ",key$biochemical),]$Unknown.Name
key$CAS<-as.character(key$CAS)

# Metabolomics key to ChEBI
key$ChEBI<-""
for(i in 1:nrow(key))
{
  found<-ChEBINames[ChEBINames$type=="CAS",]$COMPOUND_ID[match(key$CAS[i],
                                                               ChEBINames[ChEBINames$type=="CAS",]$ACCESSION_NUMBER)]
  if(!is.na(found)) key$ChEBI[i]<-found
  
  found<-ChEBINames[ChEBINames$type=="KEGG",]$COMPOUND_ID[match(key$KEGG[i],
                                                                ChEBINames[ChEBINames$type=="KEGG",]$ACCESSION_NUMBER)]
  if(!is.na(found)) key$ChEBI[i]<-found
  
  found<-ChEBINames[ChEBINames$type=="HMDB",]$COMPOUND_ID[match(key$HMDB[i],
                                                                ChEBINames[ChEBINames$type=="HMDB",]$ACCESSION_NUMBER)]
  if(!is.na(found)) key$ChEBI[i]<-found
}

df1<-read.csv("Data/scaled.csv")
rownames(df1)<-paste(gsub(" ","",df1$group),df1$timepoint,df1$ptid,sep="_")
df1$group<-df1$timepoint<-df1$ptid<-NULL

# Entropy:
m1<-as.matrix(df1)
m1<-apply(m1,2,log2)

entropDf<-data.frame(id=colnames(m1),entropy=0)
for(i in 1:ncol(m1))
{
  tryCatch(
    {
      entropDf$entropy[i]<-entropy.empirical(discretize(m1[,i],numBins=10),unit="log2")
    },
    error=function(e){print(i)}
  )
}

# Join entropy and key:
key<-key %>% left_join(entropDf)
key$include<-(key$entropy>2**.5 & key$ChEBI!="")
key<-key[key$include,]
key$ChEBI<-paste("ChEBI",key$ChEBI,sep="_")

# Subset m1 for those included:
m1<-m1[,colnames(m1) %in% key$id]

# Make distance matrix:
dists1<-dists[grepl("ChEBI",rownames(dists)),grepl("ChEBI",colnames(dists))]
dists1<-dists1[rownames(dists1)%in%key$ChEBI,colnames(dists1)%in%key$ChEBI]

########### gLasso analysis ############
glasso1<-glasso(cov(m1),rho=.5)
glasso1Omega<-glasso1$wi

bglasso1<-blockGLasso(X=m1,iterations=10,burnIn=10)
bglasso1PI<-posteriorInference(bglasso1)
idk<-bglasso1PI$posteriorMedian

colnames(idk)<-rownames(idk)<-key$biochemical[match(colnames(m1),key$id)]
