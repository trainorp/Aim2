########### Prereqs ###########
options(stringsAsFactors=FALSE,scipen=600)
library(tidyverse)
library(clusterGeneration)
library(igraph)
library(BayesianGLasso)

########### Partial correlation function ###########
pCorFun<-function(x){
  pcors<-matrix(0,nrow=nrow(x),ncol=ncol(x))
  for(j in 1:ncol(pcors)){
    for(i in 1:nrow(pcors)){
      pcors[i,j]<-(-x[i,j]/sqrt(x[i,i]*x[j,j]))
    }
  }
  return(pcors)
}

########### Simulate data ###########
nClust<-2
nInClust<-20
nSamp<-1000000
sigmas<-list()
simData<-list()
for(i in 1:nClust){
  set.seed(i+33)
  sigma<-genPositiveDefMat(nInClust,covMethod="c-vine",eta=1)$Sigma
  sigmas[[i]]<-sigma
  set.seed(i+3333)
  simData[[i]]<-mvrnorm(n=nSamp,mu=rep(0,nInClust),Sigma=sigma)
}
simData<-do.call("cbind",simData)
Sigmas<-matrix(0,nrow=ncol(simData),ncol=ncol(simData))
for(i in 1:length(sigmas)){
  start<-(i-1)*nInClust+1
  end<-(i)*nInClust
  Sigmas[start:end,start:end]<-sigmas[[i]]
}
Conc<-solve(Sigmas)
# heatmap(simData)

########### Smaller training sample ###########
set.seed(333)
Samp<-sample(1L:nrow(simData),size=50)
simDataSamp<-simData[Samp,]

########### Simulated prior information ###########
goodPrior<-Conc
ConcP<-pCorFun(Conc)

########### BGL ###########
simGrid<-expand.grid(adaptive=c(TRUE,FALSE),adaptiveType=c("norm","priorHyper"),
                     stringsAsFactors=FALSE)
simGrid<-simGrid[!(simGrid$adaptive==FALSE & simGrid$adaptiveType=="priorHyper"),]
simGrid$err<-NA
simGrid$gammaPriors<-1e4
simGrid$gammaPriors[simGrid$adaptive & simGrid$adaptiveType=="norm"]<-1e0
simGrid$gammaPriort<-1e-1
simGrid$gammaPriors[simGrid$adaptive & simGrid$adaptiveType=="norm"]
simGrid$lambdaii<-50
simGrid$lambdaii[simGrid$adaptive & simGrid$adaptiveType=="norm"]

for(i in 1:nrow(simGrid)){
  bgl1<-blockGLasso(simDataSamp,iterations=10000,burnIn=0,adaptive=simGrid$adaptive[i],
                    adaptiveType=simGrid$adaptiveType[i],
                    priorHyper=1000*(abs(goodPrior)**1.22)+0.1,
                    gammaPriors=10,gammaPriort=.1,lambdaii=10)
  bgl1PI<-posteriorInference(bgl1)
  bgl1PM<-bgl1PI$posteriorMean
  bgl1PMp<-pCorFun(bgl1PM)
  bgl1Lambda<-bgl1$lambdas[[10000]]
  simGrid$err[i]<-mean((bgl1PM-Conc)**2)
}

########### Evaluation ###########
simDataCov<-cov(simData)
simDataCor<-cor(simData)
simDataConc<-solve(cov(simData))

simDataCorG<-graph_from_adjacency_matrix(abs(simDataCor),mode="undirected",
                                        diag=FALSE,weighted=TRUE)
plot(simDataCorG)

simDataConcG<-graph_from_adjacency_matrix(abs(simDataConc),mode="undirected",
                                         diag=FALSE,weighted=TRUE)
plot(simDataConcG)

