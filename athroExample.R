options(stringsAsFactors=FALSE)
library(ChemmineR)
library(glasso)
library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(session)

setwd("~/gdrive/Dissertation/Aim2")

#Metab data
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

for(i in 1:nrow(key))
{
  if(grepl(";",key$pubchem[i]))
  {
    key$pubchem[i]<-unlist(strsplit(key$pubchem[i],";",fixed=TRUE))[1]
  }
}

include<-key$id[!grepl("Unknown",key$biochemical)]

setwd("~/gdrive/Dissertation/Aim2")

############ Structural fingerprint ############
pubChems<-as.numeric(key$pubchem)
pubChems<-pubChems[!is.na(pubChems)]
metabIds<-key$id[!is.na(key$pubchem)]
# sdfs<-list()
# for(i in 1:length(pubChems))
# {
#   tryCatch({
#     sdfs[[i]]<-getIds(pubChems[i])},
#     error=function(e) sdfs[[i]]<-"Not Here")
#   print(i)
# }
# save(sdfs,file="sdfs.RData")
load(file="sdfs.RData")

Sdfs<-sdfs[[1]]
for(i in 2:length(sdfs))
{
  Sdfs<-c(Sdfs,sdfs[[i]])
}
Sdfs@ID<-paste0("CMP",1:length(sdfs))
apSet<-sdf2ap(Sdfs)

# Annotation data:
apAnnoTemp<-sapply(datablock(Sdfs),function(x) x[c('PUBCHEM_COMPOUND_CID','PUBCHEM_IUPAC_OPENEYE_NAME','PUBCHEM_MOLECULAR_FORMULA')])
apAnnoDf<-as.data.frame(t(apAnnoTemp))
apAnnoDf$id<-metabIds
apAnnoDf$cmp<-rownames(apAnnoDf)
apAnnoDf<-apAnnoDf %>% left_join(key %>% dplyr::select(id,biochemical,HMDB))

# Pairwise similarity:
apSim<-matrix(NA,ncol=length(cid(apSet)),nrow=length(cid(apSet)))
rownames(apSim)<-colnames(apSim)<-cid(apSet)
for(i in 1:nrow(apSim))
{
  for(j in 1:ncol(apSim))
  {
    apSim[i,j]<-cmp.similarity(apSet[i],apSet[j]) 
  }
}

distForPlot<-1-apSim
rownames(distForPlot)<-apAnnoDf$biochemical[match(rownames(distForPlot),apAnnoDf$cmp)]
colnames(distForPlot)<-apAnnoDf$biochemical[match(colnames(distForPlot),apAnnoDf$cmp)]

# See heatmap3.R for making heatmap with distForPlot

rm(i,j,include,metabIds,Sdfs,sdfs,apAnnoTemp,spec,pubChems)

# Distribution of distances:
dists<-1-apSim
diag(dists)<-NA
distsV<-c(dists)
distsV<-distsV[!is.na(distsV)]

# Match row/column names:
rownames(dists)<-apAnnoDf$id[match(rownames(dists),apAnnoDf$cmp)]
colnames(dists)<-apAnnoDf$id[match(colnames(dists),apAnnoDf$cmp)]

############ Abundance data ############
df1<-read.csv("~/gdrive/AthroMetab/Data/scaled.csv")
rownames(df1)<-paste(df1$group,df1$timepoint,df1$ptid,sep="_")

# Follow-up and with annotation only:
df1<-df1[df1$timepoint=="TF-U",]
m1<-as.matrix(df1[,!names(df1) %in% c("group","timepoint","ptid")])
m1<-scale(m1,center=TRUE,scale=FALSE)

# Entropy filter:
m1<-m1[,apply(m1,2,function(x) length(unique(x))>19)]

#save(m1,file="~/gdrive/Dissertation/Aim2/m1.RData")

# Make prior info matrix:
priorHyper=matrix(1e-6,ncol=ncol(dists),nrow=nrow(dists))
colnames(priorHyper)<-rownames(priorHyper)<-colnames(dists)
for(i in 1:nrow(dists))
{
  for(j in 1:ncol(dists))
  {
    priorHyper[i,j]<-priorHyper[i,j]*dists[i,j]
  }
}

priorHyper2<-matrix(NA,ncol=ncol(m1),nrow=ncol(m1))
colnames(priorHyper2)<-rownames(priorHyper2)<-colnames(m1)
# Fix prior info matrix to have same cols as m1
for(j in 1:ncol(priorHyper2))
{
  for(i in 1:nrow(priorHyper2))
  {
    if(rownames(priorHyper2)[i] %in% rownames(dists) & colnames(priorHyper2)[j] %in% colnames(dists))
    {
      priorHyper2[i,j]<-dists[rownames(priorHyper2)[i],colnames(priorHyper2)[j]]
    }
  }
}
rm(priorHyper)
priorHyper<-priorHyper2
# png(file="PriorMatrixImg.png",height=5,width=6,units="in",res=600)
# image(priorHyper2,col=rev(blue2green2red(100)))
# dev.off()

# Include only those with structure:
m2<-m1
m2<-m2[,colnames(m2) %in% 
 colnames(priorHyper2)[!apply(priorHyper2,2,FUN=function(x) sum(as.integer(is.na(x)))==length(x))]]
priorHyper3<-priorHyper2[match(colnames(m2),rownames(priorHyper2)),
                         match(colnames(m2),colnames(priorHyper2))]
save.session(file="athroExWorkspace_20171011.RData")

############ Actually running the sampler ############
library(session)
restore.session(file="~/gdrive/Dissertation/Aim2/athroExWorkspace_20171011.RData")
# Other arguments for adaptive lasso:
source('~/gdrive/Dissertation/Aim2/BayesianGLasso/R/blockGLasso.default.R')
source('~/gdrive/Dissertation/Aim2/BayesianGLasso/R/blockAdGLasso.default.R')

# ptm<-proc.time()
# idk<-blockAdGLasso.default(m2,iterations=1000,burnIn=100,gammaPriort=1e-6,gammaPriors=10**.75,lambdaii=10**.75,
#                            adaptive=TRUE,adaptiveType="priorHyper",
#                            priorHyper=2*(1-priorHyper3),verbose=TRUE)
# save(idk,file="~/gdrive/Dissertation/Aim2/idk1000.RData")
# proc.time()-ptm

# ptm<-proc.time()
# idk<-blockAdGLasso.default(m2,iterations=1000,burnIn=100,gammaPriort=1e-6,gammaPriors=1e3,lambdaii=1e3,
#                            adaptive=TRUE,adaptiveType="priorHyper",
#                            priorHyper=2*(1-priorHyper3),verbose=TRUE)
# save(idk,file="~/gdrive/Dissertation/Aim2/idk2.RData")
# proc.time()-ptm

#load("~/gdrive/Dissertation/Aim2/idk.RData")
load("~/gdrive/Dissertation/Aim2/out_i1000_b100_s5p6_l5p6_t1e-06.RData")

idk2<-idk$Omegas[[1000]]
idk3<-idk2
for(j in 1:ncol(idk3))
{
  for(i in 1:nrow(idk3))
  {
    idk3[i,j]<-(-idk2[i,j] / sqrt(idk2[i,i]*idk2[j,j]))
  }
}
rownames(idk3)<-colnames(idk3)<-colnames(idk$Sigmas[[2]])

which.max(idk3["M231",])

caffeine<-sapply(idk$Omegas,FUN=function(x) x["M231","M231"])
pzanth<-sapply(idk$Omegas,FUN=function(x) x["M526","M526"])
arach<-sapply(idk$Omegas,FUN=function(x) x["M22","M22"])
ca<-sapply(idk$Omegas,FUN=function(x) x["M231","M22"])
cp<-sapply(idk$Omegas,FUN=function(x) x["M231","M526"])

par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(4,3,2,3))
plot(1:1100,caffeine,type="l",xaxt="n",las=1)
plot(1:1100,pzanth,type="l",xaxt="n",axes=FALSE,frame.plot=TRUE)
axis(4,las=1)
plot(1:1100,cp,type="l",xaxt="n",las=1)
plot(1:1100,-cp/sqrt(caffeine*pzanth),type="l",xaxt="n",axes=FALSE,frame.plot=TRUE)
axis(4,las=1)

plot(1:1100,caffeine,type="l")
plot(1:1100,arach,type="l")
plot(1:1100,ca,type="l")
plot(1:1100,-ca/sqrt(caffeine*arach),type="l")

############ Posterior Inference ############
source('~/gdrive/Dissertation/Aim2/BayesianGLasso/R/postInference.R')
pI100<-posteriorInference(idk,alpha=.5)
pI100Median<-pI100$posteriorMedian
rownames(pI100Median)<-colnames(pI100Median)<-key$biochemical[match(colnames(m2),key$id)]

pI100MedianRho<-matrix(NA,nrow=nrow(pI100Median),ncol=ncol(pI100Median))
rownames(pI100MedianRho)<-colnames(pI100MedianRho)<-key$biochemical[match(colnames(m2),key$id)]
for(j in 1:ncol(pI100MedianRho))
{
  for(i in 1:nrow(pI100MedianRho))
  {
    pI100MedianRho[i,j]<-(-pI100Median[i,j] / sqrt(pI100Median[i,i]*pI100Median[j,j]))
  }
}
diag(pI100MedianRho)<-NA

CIzeros<-matrix(0,nrow=nrow(pI100Median),ncol=ncol(pI100Median))
rownames(CIzeros)<-colnames(CIzeros)<-colnames(pI100Median)
for(j in 1:ncol(CIzeros))
{
  for(i in 1:nrow(CIzeros))
  {
    CIzeros[i,j]<-!(pI100$lowerCI[i,j]<0 & 0<pI100$upperCI[i,j])
  }
}
diag(CIzeros)<-0
CIzeros2<-CIzeros[rowSums(CIzeros)>0,colSums(CIzeros)>0]

############ Graphical model ############
library(igraph)
diag(pI100MedianRho)<-0
g1<-graph.adjacency(pI100MedianRho,mode="undirected",weighted=TRUE)
g2<-graph.adjacency(CIzeros,mode="undirected")

g1b<-delete.edges(g1,which(abs(E(g1)$weight)<.05))

plot(g1b,vertex.size=2,vertex.label=NA)
png(file="g1b.png",height=6,width=6,units="in",res=600)
plot(g1b,vertex.size=2,vertex.label.cex=.1)
dev.off()

myRed<-rgb(255,0,0,max=255,alpha=200,names="myRed")
myBlue<-rgb(0,0,200,max=255,alpha=200,names="myBlue")
E(g1b)$color<-ifelse(E(g1b)$weight<0,"#FF0000C8","#0000C8C8")

png(file="g1b.png",height=12,width=12,units="in",res=600)
plot(g1b,vertex.size=2,vertex.label.cex=.1,layout=layout_with_graphopt,
     edge.width=abs(E(g1b)$weight)*5,vertex.color="grey")
dev.off()

plot(g1b,vertex.size=2,vertex.label=NA,layout=layout_with_fr)
plot(g1b,vertex.size=2,vertex.label=NA,layout=layout_with_graphopt)
plot(g1b,vertex.size=2,vertex.label=NA,layout=layout_with_dh)
plot(g1b,vertex.size=2,vertex.label=NA,layout=layout_nicely)

gl1<-glasso(s=cov(m1),rho=.5)$wi
gl1Rho<-matrix(NA,nrow=nrow(gl1),ncol=ncol(gl1))
rownames(gl1Rho)<-colnames(gl1Rho)<-colnames(pI100Median)
for(j in 1:ncol(gl1Rho))
{
  for(i in 1:nrow(gl1Rho))
  {
    gl1Rho[i,j]<-(-gl1[i,j] / sqrt(gl1[i,i]*gl1[j,j]))
  }
}
gl1Rho2<-gl1Rho[rowSums(gl1Rho)>0,colSums(gl1Rho)>0]
