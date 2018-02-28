options(stringsAsFactors=FALSE)
library(MASS)
library(BayesianGLasso)
library(tidyverse)

setwd("~/gdrive/Dissertation/Aim2")

############ Untargeted data ############
setwd("~/gdrive/AthroMetab/")
spec<-read.csv("Data/AthroACSRawSpectra.csv")
key<-read.csv("Data/metabolite_key2.csv")
key[key$id=="M801",]$super<-"Xenobiotics"
key$biochemical<-as.character(key$biochemical)
key$Unknown.Name<-as.character(key$Unknown.Name)
key[grepl("X - ",key$biochemical),]$biochemical<-key[grepl("X - ",key$biochemical),]$Unknown.Name
key<-spec %>% left_join(key,c("comp_id"="MEBID"))
key$pubchem<-as.character(key$pubchem)

for(i in 1:nrow(key)){
  if(grepl(";",key$pubchem[i])){
    key$pubchem[i]<-unlist(strsplit(key$pubchem[i],";",fixed=TRUE))[1]
  }
}

include<-key$id[!grepl("Unknown",key$biochemical)]

setwd("~/gdrive/Dissertation/Aim2")

############ Structural similarity ############
sims<-read.table("athroMetab_color_tanimotos.txt")
sims<-sims[!(sims$V1 %in% key$id[grepl("Tentative",key$biochemical)] | sims$V2 %in% key$id[grepl("Tentative",key$biochemical)]),]
simMat<-matrix(NA,nrow=length(unique(sims$V1)),ncol=length(unique(sims$V1)))
rownames(simMat)<-colnames(simMat)<-unique(sims$V1)
for(i in 1:nrow(sims)){
  simMat[match(sims$V1[i],rownames(simMat)),match(sims$V2[i],colnames(simMat))]<-sims$V3[i]
  print(i)
}
simMat2<-simMat
rownames(simMat2)<-key$biochemical[match(rownames(simMat2),key$id)]
colnames(simMat2)<-key$biochemical[match(colnames(simMat2),key$id)]

source('~/gdrive/Dissertation/Aim2/heatmap3.R')
# heatmap3(1-simMat2) # May need to go back and change the resolution

############ Abundance data ############
df1<-read.csv("~/gdrive/AthroMetab/Data/scaled.csv")
rownames(df1)<-paste(df1$group,df1$timepoint,df1$ptid,sep="_")

# Follow-up and with annotation only:
df1<-df1[df1$timepoint=="TF-U",]
m1<-as.matrix(df1[,!names(df1) %in% c("group","timepoint","ptid")])
m1<-scale(m1,center=TRUE,scale=FALSE)

# Entropy filter:
m1<-m1[,apply(m1,2,function(x) length(unique(x))>19)]

# Filter for those without structural information:
m1<-m1[,colnames(m1) %in% include]

# Random sample:
set.seed(3)
idk<-blockGLasso(m1[,sample(1:ncol(m1),size=200)],iterations=25,burnIn=0)

