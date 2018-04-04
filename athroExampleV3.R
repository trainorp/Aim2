########### Prereqs ###########
options(stringsAsFactors=FALSE)
library(MASS)
library(BayesianGLasso)
library(tidyverse)
library(igraph)
library(RCy3)

setwd("~/gdrive/Dissertation/Aim2")

## Begin don't run:
########### Partial correlation Function ###########
pCorFun<-function(x){
  pcors<-matrix(0,nrow=nrow(x),ncol=ncol(x))
  for(j in 1:ncol(pcors)){
    for(i in 1:nrow(pcors)){
      pcors[i,j]<-(-x[i,j]/sqrt(x[i,i]*x[j,j]))
    }
  }
  return(pcors)
}

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
sims<-sims[!(sims$V1=="M681" | sims$V2=="M681"),]
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

# Follow-up and with annotation only:
df1<-df1[df1$timepoint=="TF-U" | (df1$group=="sCAD" & df1$timepoint=="T0"),]
m1<-as.matrix(df1[,!names(df1) %in% c("group","timepoint","ptid")])
m1<-scale(m1,center=TRUE,scale=FALSE)

# Entropy filter:
m1<-m1[,apply(m1,2,function(x) length(unique(x))>14)]

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
heatmap3(1-simMat2,labRow="",labCol="",nam="NoLab")
heatmap3(1-simMat2,resMult=1.25,nam="Lab")

# Cut tree for smaller cluster picture:
hc<-hclust(as.dist(1-simMat2),method="ward.D2")
plot(hc,cex=.2)
ct<-data.frame(clust=cutree(hc,h=1.01))

# Small dataset:
simMat3<-simMat2
simMat3<-simMat3[rownames(simMat3) %in% rownames(ct)[ct$clust==9],colnames(simMat3) %in% rownames(ct)[ct$clust==9]]
colnames(simMat3)<-rownames(simMat3)<-gsub("(alpha or beta)","",colnames(simMat3))
heatmap3(1-simMat3,nam="Sub",lwdPar=2,cFac=6)

########### Run the sampler: ###########
# Structure Adaptive:
priorHyper<-2*(simMat**2)+.1
save.image(file="atheroExampleV3Data.RData")
## End don't run:

ptm<-proc.time()
aiBGL1<-blockGLasso(m1,iterations=10,burnIn=0,adaptive=TRUE,
                    adaptiveType="priorHyper",priorHyper=priorHyper,
                    lambdaii=62,gammaPriors=8,gammaPriort=.001)
proc.time()-ptm
summary(c(aiBGL1$lambdas[[2]])[c(aiBGL1$lambdas[[2]])>0])
save(aiBGL1,file="aiBGL1.RData")

########### Posterior Inference ###########
# Import data:
load(file="atheroExampleV3Data.RData")
load(file="aiBGL1.RData")

# Concentration matrix:
pIaiBGL1<-posteriorInference(aiBGL1)
aiBGL1Con<-pIaiBGL1$posteriorMedian
colnames(aiBGL1Con)<-rownames(aiBGL1Con)<-key$biochemical[match(colnames(m1),key$id)]

# Partial correlation matrix:
aiBGL1$Omegas<-lapply(aiBGL1$Omegas,pCorFun)
pIaiBGL1<-posteriorInference(aiBGL1)
aiBGL1Cor<-pIaiBGL1$posteriorMedian
colnames(aiBGL1Cor)<-rownames(aiBGL1Cor)<-key$biochemical[match(colnames(m1),key$id)]
rm(aiBGL1)

# Save image:
save.image(file="atheroExampleV3DataPart2.RData")

########### Cytoscape Graph ###########
load(file="atheroExampleV3DataPart2.RData")

# Weird name issue:
colnames(aiBGL1Con)[colnames(aiBGL1Con)=="tryptophan betaine "]<-"tryptophan betaine"
colnames(aiBGL1Cor)[colnames(aiBGL1Cor)=="tryptophan betaine "]<-"tryptophan betaine"
rownames(aiBGL1Con)[rownames(aiBGL1Con)=="tryptophan betaine "]<-"tryptophan betaine"
rownames(aiBGL1Cor)[rownames(aiBGL1Cor)=="tryptophan betaine "]<-"tryptophan betaine"

# Cytoscape graph function:
graphFun<-function(mat){
  # Graph from adjacency 
  g<-graph_from_adjacency_matrix(abs(get(mat)),mode="undirected",diag=FALSE,weighted=TRUE)
  e<-as.data.frame(get.edgelist(g))
  E(g)$color<-c("darkred","navyblue")[as.integer(get(mat)[lower.tri(get(mat))]>0)+1L]
  e$col<-get.edge.attribute(g,"color")
  
  # Delete edges:
  g<-delete_edges(g,which(E(g)$weight<.01))

  # Graph as graphNEL
  gNel<-as_graphnel(g)
  gNel<-initEdgeAttribute(gNel,'weight','numeric',0)
  gNel<-initEdgeAttribute(gNel,"color","char","black")
  
  # Export Weights attribute:
  w<-data.frame(weights=unlist(edgeData(gNel,attr="weight")))
  w$name<-rownames(w)
  w$name<-gsub("\\|"," (unspecified) ",w$name)
  write.csv(w,file=paste0(mat,"Weights",".csv"),row.names=FALSE)
  
  # Export Color attribute:
  co<-data.frame(colors=unlist(edgeData(gNel,attr="color")))
  co$name<-rownames(co)
  co$name<-gsub("\\|"," (unspecified) ",co$name)
  write.csv(co,file=paste0(mat,"Colors",".csv"),row.names=FALSE)
  
  return(gNel)
}

# Cytoscape Windows:
deleteAllWindows(CytoscapeConnection())
cw<-CytoscapeWindow('aiBGL1Cor',graph=graphFun("aiBGL1Cor"),overwrite=TRUE)
displayGraph(cw)
cw2<-CytoscapeWindow('aiBGL1Con',graph=graphFun("aiBGL1Con"),overwrite=TRUE)
displayGraph(cw2)

########### igraph Graph ###########
load(file="atheroExampleV3DataPart2.RData")

g<-graph_from_adjacency_matrix(abs(aiBGL1Cor),mode="undirected",diag=FALSE,weighted=TRUE)
