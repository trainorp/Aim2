########### Prereqs ###########
options(stringsAsFactors=FALSE)
library(MASS)
library(BayesianGLasso)
library(tidyverse)

setwd("~/gdrive/Dissertation/Aim2")

## Begin don't run:
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

############ Abundance data ############
df1<-read.csv("~/gdrive/AthroMetab/Data/scaled.csv")
rownames(df1)<-paste(df1$group,df1$timepoint,df1$ptid,sep="_")
rm(spec,i)
save.image(file="atheroExampleV3Data.RData")
## End don't run:

# Import data:
load(file="atheroExampleV3Data.RData")

# Follow-up and with annotation only:
df1<-df1[df1$timepoint=="TF-U",]
m1<-as.matrix(df1[,!names(df1) %in% c("group","timepoint","ptid")])
m1<-scale(m1,center=TRUE,scale=FALSE)

# Entropy filter:
m1<-m1[,apply(m1,2,function(x) length(unique(x))>19)]

# Filter for those without structural information:
m1<-m1[,colnames(m1) %in% colnames(simMat)]

# Make sure column order / names match:
simMat<-simMat[rownames(simMat) %in% colnames(m1),colnames(simMat) %in% colnames(m1)]

########### Heatmap ###########
simMat2<-simMat
rownames(simMat2)<-key$biochemical[match(rownames(simMat2),key$id)]
colnames(simMat2)<-key$biochemical[match(colnames(simMat2),key$id)]

source('~/gdrive/Dissertation/Aim2/heatmap3.R')
dev.new()
heatmap3(1-simMat2) # May need to go back and change the resolution

########### Random sample: ###########
# True partial correlations:
pCorFun<-function(x){
  pcors<-matrix(0,nrow=nrow(x),ncol=ncol(x))
  for(j in 1:ncol(pcors)){
    for(i in 1:nrow(pcors)){
      pcors[i,j]<-(-x[i,j]/sqrt(x[i,i]*x[j,j]))
    }
  }
  return(pcors)
}

set.seed(3)
samp<-sample(1:ncol(m1),size=200)
idk<-blockGLasso(m1[,samp],iterations=200,burnIn=0)
idk2<-idk$Omegas[[200]]
colnames(idk2)<-rownames(idk2)<-key$biochemical[match(colnames(m1)[samp],key$id)]
idk3<-pCorFun(idk2)
rownames(idk3)<-colnames(idk3)<-colnames(idk2)

priorHyper<-simMat+.1
colnames(m1) == colnames(priorHyper) # Yep
aiBGL1<-blockGLasso(m1[,samp],iterations=10,burnIn=0,adaptive=TRUE,
                    adaptiveType="priorHyper",priorHyper=priorHyper[samp,samp],
                    gammaPriors=10,gammaPriort=.001)

library(RhpcBLASctl)
ptm<-proc.time()
aiBGL1<-blockGLasso(m1,iterations=5,burnIn=0,adaptive=TRUE,
                    adaptiveType="priorHyper",priorHyper=priorHyper,
                    gammaPriors=40,gammaPriort=.001)
proc.time()-ptm

aiBGL2<-aiBGL1$Omegas[[10]]
aiBGL3<-pCorFun(aiBGL2)
colnames(aiBGL3)<-rownames(aiBGL3)<-key$biochemical[match(colnames(m1),key$id)]
