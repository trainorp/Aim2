########### Prereqs ###########
options(stringsAsFactors=FALSE,scipen=600)
library(tidyverse)
library(clusterGeneration)
library(igraph)
library(BayesianGLasso)

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

########### BGL ###########
bgl1<-blockGLasso(simData,iterations=1000,burnIn=0,adaptive=FALSE)
bgl1Omegas<-bgl1$Omegas[[1000]]

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

