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
simDataTest<-simData[!(1L:nrow(simData) %in% Samp),]

########### Simulated prior information ###########
goodPrior<-Conc
ConcP<-pCorFun(Conc)

########### BGL ###########
simGridPriorAd<-data.frame(adaptive=TRUE,adaptiveType="priorHyper",
                           gammaPriors=10**seq(0,4,1),gammaPriort=1e-1,lambdaii=10)
simGridPriorNorm<-data.frame(adaptive=TRUE,adaptiveType="norm",
                           gammaPriors=10**seq(0,4,1),gammaPriort=1e-1,lambdaii=10)
simGridReg<-data.frame(adaptive=FALSE,adaptiveType="norm",gammaPriors=NA,
                       gammaPriort=NA,lambdaii=NA)

simGrid<-rbind(simGridPriorAd,simGridPriorNorm,simGridReg)
simGrid$likelihood<-NA
for(i in 1:nrow(simGrid)){
  bgl1<-blockGLasso(simDataSamp,iterations=1000,burnIn=0,adaptive=simGrid$adaptive[i],
                    adaptiveType=simGrid$adaptiveType[i],
                    priorHyper=100*(abs(goodPrior))+0.1,
                    gammaPriors=simGrid$gammaPriors[i],
                    gammaPriort=simGrid$gammaPriort[i],
                    lambdaii=simGrid$lambdaii[i])
  bgl1PI<-posteriorInference(bgl1)
  bgl1PM<-bgl1PI$posteriorMean
  bgl1PMp<-pCorFun(bgl1PM)
  bgl1Lambda<-bgl1$lambdas[[1000]]
  simGrid$likelihood[i]<-sum(mvtnorm::dmvnorm(x=simDataTest,
              mean=rep(0,ncol(simDataTest)),sigma=solve(bgl1PM),log=TRUE))
}

########### Data likelihood function ###########
sum(mvtnorm::dmvnorm(x=simDataTest,mean=rep(0,ncol(simDataTest)),
                 sigma=solve(bgl1PM),log=TRUE))

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

