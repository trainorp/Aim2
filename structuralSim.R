############ Prereqs ############
options(stringsAsFactors=FALSE)
library(glasso)
library(colorRamps)
library(igraph)
library(dplyr)
library(tidyr)
oldPar<-par()

############ Untargeted data ############
setwd("~/gdrive/AthroMetab/Data")
spec<-read.csv("AthroACSRawSpectra.csv")
key<-read.csv("metabolite_key2.csv")
key[grepl("X - ",key$biochemical),]$biochemical<-key[grepl("X - ",key$biochemical),]$Unknown.Name
key<-spec %>% left_join(key,c("comp_id"="MEBID"))
metab<-read.csv("scaled.csv")
metab$timepoint[metab$timepoint=="TF-U"]<-"TF/U"
rownames(metab)<-paste(metab$ptid,metab$timepoint)
metabPheno<-metab %>% select(group,ptid,timepoint)
metab<-metab %>% select(-group,-ptid,-timepoint)
metab<-log2(metab) 
metab<-metab[grepl("TF/U",rownames(metab)),]
for(i in 1:nrow(key))
{
  if(grepl(";",key$pubchem[i]))
  {
    key$pubchem[i]<-unlist(strsplit(key$pubchem[i],";",fixed=TRUE))[1]
  }
}

############ Max Imp filter #############
metab<-metab[,apply(metab,2,function(x) table(x)[1]/32<.6)]

############ PubChem #############
setwd("~/gdrive/Dissertation/KBRIN_UT_2017")

# Get Fingerprints:
# pubChems<-as.numeric(key$pubchem)
# pubChems<-pubChems[!is.na(pubChems)]
# sdfs<-list()
# for(i in 1:length(pubChems))
# {
#   tryCatch({
#     sdfs[[i]]<-getIds(pubChems[i])},
#     error=function(e) sdfs[[i]]<-"fuck")
#   print(i)
# }
# save(sdfs,file="sdfs.RData")
load("sdfs.RData")

# Concatenate sdfs:
sdfSet<-sdfs[[1]]
for(i in 2:length(sdfs))
{
  if(!is.null(sdfs[[i]]))
  {
    sdfSet<-c(sdfSet,sdfs[[i]])
  }
}

# Get names:
pcNames<-c()
for(i in 1:length(sdfSet))
{
  pcNames<-c(pcNames,sdfSet[[i]]@header['Molecule_Name'])
}
sdfSet@ID<-as.character(pcNames)

# Make atom pair set:
apSet<-sdf2ap(sdfSet)
fpSet<-desc2fp(apSet,descnames=4096)

# Similarities
simMat<-sapply(cid(fpSet), function(x) fpSim(x=fpSet[x], fpSet, sorted=FALSE)) 

# Distance matrix:
distMat<-1-simMat

# Names:
keyWpubchem<-key[!is.na(key$pubchem),]
rownames(distMat)<-colnames(distMat)<-
  keyWpubchem$biochemical[match(rownames(distMat),keyWpubchem$pubchem)]

# Clustering: 
hc<-hclust(as.dist(distMat),method="ward.D2")
clust<-cutreeHybrid(hc,distMat,minClusterSize=4,
                    deepSplit=TRUE,cutHeight=1,pamStage=FALSE)

# Dendrogram
#png("chemClust.png",height=4,width=12,units="in",res=600)
par(mar=c(0,0,0,0))
plotDendroAndColors(hc,labels2colors(clust$labels),cex.dendroLabels=.15,
                    groupLabels="",main="",marAll=c(1,4,1,0))
#dev.off()
par(oldPar)

# Heatmap
#png("chemHeat.png",height=12,width=12,units="in",res=600)
#par(bg="transparent")
hm<-heatmap(distMat,cexRow=.175,cexCol=.175,col=rev(matlab.like2(400)[c(seq(1,300,5),301:400)]),
        hclustfun=function(x) hclust(as.dist(x),method="ward.D2"))
#dev.off()
par(oldPar)

# png("chemHeatBig.png",height=15,width=15,units="in",res=1500)
# heatmap(distMat,cexRow=.15,cexCol=.15,col=rev(matlab.like2(400)[c(seq(1,300,5),301:400)]),
#             hclustfun=function(x) hclust(as.dist(x),method="ward.D2"))
# dev.off()

############ Make penalty matrix ############
rownames(distMat)<-colnames(distMat)<-key$id[match(rownames(distMat),keyWpubchem$biochemical)]
covM<-cov(metab)
keepNames<-intersect(colnames(metab),colnames(distMat))
distMat<-distMat[keepNames,keepNames]
covM<-covM[keepNames,keepNames]

############ gLasso ############
source('~/gdrive/Dissertation/KBRIN_UT_2017/glasso2.R')

# First example:
gLassoDist<-glasso2(covM,rho=2*distMat,penalize.diagonal=FALSE)
wDist<-gLassoDist$w
wiDist<-gLassoDist$wi
rownames(wDist)<-colnames(wDist)<-colnames(covM)

# png(file="graphLayout1.png",height=5,width=7,units="in",res=600)
# par(mfrow=c(2,3),mar=c(1,1,2,1))
# plot(graph1,layout=layout_as_tree, vertex.size=1, vertex.label=NA,main="Tree")
# plot(graph1,layout=layout_in_circle,vertex.size=1,vertex.label=NA,main="Circle")
# plot(graph1,layout=layout_nicely, vertex.size=1, vertex.label=NA,main="Nicely")
# plot(graph1,layout=layout_with_dh, vertex.size=1, vertex.label=NA,main="DH")
# plot(graph1,layout=layout_with_gem, vertex.size=1, vertex.label=NA,main="Gem")
# plot(graph1,layout=layout_with_graphopt, vertex.size=1, vertex.label=NA,main="GraphOpt")
# dev.off()
# 
# png(file="graphLayout2.png",height=2.5,width=4.67,units="in",res=600)
# par(mfrow=c(2,2),mar=c(1,1,2,1))
# plot(graph1,layout=layout_with_mds,vertex.size=1,vertex.label=NA,main="MDS")
# plot(graph1,layout=layout_with_fr, vertex.size=1, vertex.label=NA,main="FR")
# plot(graph1,layout=layout_with_kk,vertex.size=1,vertex.label=NA,main="KK")
# plot(graph1,layout=layout_with_lgl,vertex.size=1,vertex.label=NA,main="LGL")
# dev.off()

############ Adaptive gLasso ############
BICdf<-data.frame(rhoPen=seq(.5,5,.01),logLik=NA,nEdges=NA)
for(i in 1:nrow(BICdf))
{
  gLasso<-glasso2(covM,rho=BICdf$rhoPen[i]*distMat,penalize.diagonal=FALSE)
  BICdf$logLik[i]<-gLasso$loglikNoPen
  diag(gLasso$wi)<-0
  BICdf$nEdges[i]<-sum(gLasso$wi>0)
  print(i)
}
BICdf$BIC<-(-2)*BICdf$logLik+BICdf$nEdges*log(32)
BICdf$eBIC<-(-2)*BICdf$logLik+BICdf$nEdges*log(32)+.06*BICdf$nEdges*1032
plot(BICdf$rhoPen,BICdf$eBIC)
# pen<-1.87
pen<-BICdf$rhoPen[which.min(BICdf$eBIC)]
# number of edges: 662
BICdf %>% filter(rhoPen==pen)

# Final model
gLassoDist<-glasso2(covM,rho=pen*distMat,penalize.diagonal=FALSE)
wDist<-gLassoDist$w
wiDist<-gLassoDist$wi
rownames(wDist)<-colnames(wDist)<-colnames(covM)
rownames(wiDist)<-colnames(wiDist)<-colnames(covM)

############ Regular gLasso ############
BICRegdf<-data.frame(rhoPen=seq(.5,5,.01),logLik=NA,nEdges=NA)
for(i in 1:nrow(BICRegdf))
{
  gLasso<-glasso2(covM,rho=BICRegdf$rhoPen[i],penalize.diagonal=FALSE)
  BICRegdf$logLik[i]<-gLasso$loglikNoPen
  diag(gLasso$wi)<-0
  BICRegdf$nEdges[i]<-sum(gLasso$wi>0)
  print(i)
}
BICRegdf$BIC<-(-2)*BICRegdf$logLik+BICRegdf$nEdges*log(32)
BICRegdf$eBIC<-(-2)*BICRegdf$logLik+BICRegdf$nEdges*log(32)+.06*BICRegdf$nEdges*1032
BICRegdf$rhoPen[which.min(BICRegdf$eBIC)]

# Final model
gLasso<-glasso2(covM,rho=1.268,penalize.diagonal=FALSE)
w<-gLasso$w
wi<-gLasso$wi
rownames(w)<-colnames(w)<-colnames(covM)
rownames(wi)<-colnames(wi)<-colnames(covM)

############ Wi versus structural similarity ############
png(file="dotMat.png",width=8,height=4.25,units="in",res=600)
par(mfrow=c(1,2),mar=c(1,1,2,0))
wiDistOrd<-wiDist[rownames(distMat),colnames(distMat)]>0
wiDistOrd<-wiDistOrd[hm$rowInd,hm$colInd]
image(1L:ncol(wiDistOrd),1L:nrow(wiDistOrd),wiDistOrd,col=c("white","black"),axes=FALSE,
      xlab="",ylab="")
box()
mtext("Structure-Adaptive gLasso")
par(mar=c(1,0,2,1))
wiOrd<-wi[rownames(distMat),colnames(distMat)]>0
wiOrd<-wiOrd[hm$rowInd,hm$colInd]
image(1L:ncol(wiOrd),1L:nrow(wiOrd),wiOrd,col=c("white","black"),axes=FALSE,
      xlab="",ylab="")
mtext("gLasso")
box()
dev.off()
par(oldPar)

############ Make some graphs ############
wiDist2<-wiDist
diag(wiDist2)<-0
graph2<-graph.adjacency(wiDist2,mode="undirected",weighted=TRUE)

# Edge widths:
ew<-abs(E(graph2)$weight)*5
eCDF<-ecdf(ew)
widths<-qunif(eCDF(ew),min=.25,max=1.5)**1.6
myRed<-rgb(255,0,0,max=255,alpha=200,names="myRed")
myBlue<-rgb(0,0,200,max=255,alpha=200,names="myBlue")
E(graph2)$color<-ifelse(E(graph2)$weight>0,"#FF0000C8","#0000C8C8")

# Regular Size
png(file="plasmaInteractome.png",res=600,height=6,width=6,units="in")
par(mar=c(0,0,1.5,0))
set.seed(21)
plot(graph2,layout=layout_nicely,vertex.size=2,vertex.label=NA,main="Plasma Interactome",
     vertex.color="grey",edge.width=widths)
dev.off()

# Large!
vNames<-key$biochemical[match(colnames(covM),key$id)]
png(file="plasmaInteractomeLarge.png",res=600,height=20,width=20,units="in")
par(mar=c(0,0,1.5,0))
set.seed(33)
plot(graph2,layout=layout_nicely,vertex.size=0,vertex.label=vNames,main="Plasma Interactome",
     vertex.color="grey",edge.width=widths,vertex.label.cex=.4)
dev.off()

############ Cotinine ############
cotWi<-wiDist
diag(cotWi)<-0
rownames(cotWi)<-colnames(cotWi)<-colnames(covM)
cotID<-which(rownames(covM)=="M254")
cotMask<-abs(wiDist[cotID,])>1e-16
cotID2<-names(which(cotMask))
cotWi<-cotWi[cotID2,cotID2]

cotVNames<-key$biochemical[match(colnames(cotWi),key$id)]
cotGraph<-graph.adjacency(cotWi,mode="undirected",weighted=TRUE)
cotEW<-abs(E(cotGraph)$weight)*5
coteCDF<-ecdf(cotEW)
cotWidths<-qunif(coteCDF(cotEW),min=.2,max=1.2)**7
E(cotGraph)$color<-ifelse(E(cotGraph)$weight>0,"#FF000064","#0000C864")

png(file="Cotinine.png",height=6,width=6,units="in",res=600)
par(mar=c(1,1,2,1))
set.seed(3)
plot(cotGraph,layout=layout_nicely,vertex.size=2,vertex.label=cotVNames,
     main="Cotinine Interactome",
     vertex.color=rgb(190,190,190,max=255,alpha=100),
     edge.width=cotWidths,vertex.label.cex=.6,vertex.label.color="black",
     vertex.frame.color=rgb(190,190,190,max=255,alpha=100))
dev.off()

table(apply(cotWi,2,function(x) sum(x))==0)

############ Cytoscape ############
library(NetPathMiner)
cotGraph2<-cotGraph
E(cotGraph2)$weight<-abs(E(cotGraph2)$weight)
plotCytoscapeGML(cotGraph2,file="idk.gml",vertex.label=pregVNames)
