########### Prereqs ###########
library(MASS)
library(dplyr)
library(tidyr)
library(BayesianGLasso)
library(igraph)

setwd("~/gdrive/Dissertation/Aim2")

# AR(1) example
sig<-toeplitz(.9**(0:9)) #toeplitz(.9**(0:109))
#sig[1,10]<-.9
omega<-solve(sig)
pcors<-matrix(0,nrow=nrow(omega),ncol=ncol(omega))
for(j in 1:ncol(pcors)){
  for(i in 1:nrow(pcors)){
    pcors[i,j]<-(-omega[i,j]/sqrt(omega[i,i]*omega[j,j]))
  }
}
set.seed(3)
x1<-mvrnorm(n=20,mu=rep(0,ncol(sig)),Sigma=sig)

########### Regular BGL Test ############
BGL1<-blockGLasso(x1,iterations=1000,burnIn=500,adaptive=FALSE,lambdaPriora=1,
            lambdaPriorb=1/10)

pIBGL1<-posteriorInference(BGL1)
medBGL1<-pIBGL1$posteriorMedian
gBGL1<-graph_from_adjacency_matrix(abs(medBGL1),mode="undirected",diag=FALSE,weighted=TRUE)
E(gBGL1)$width<-(E(gBGL1)$weight**2)/4
medBGL1[lower.tri(medBGL1,diag=TRUE)]<-NA
na.omit(c(medBGL1))
E(gBGL1)$color<-c("darkred","navyblue")[as.integer(na.omit(c(t(medBGL1)))>0)+1L]
plot(gBGL1)

aBGL1<-blockGLasso(x1,iterations=1000,burnIn=500,adaptive=TRUE,adaptiveType="norm",
                   gammaPriors=10**(-2),gammaPriort=10**(-6))
pIaBGL1<-posteriorInference(aBGL1)
medaBGL1<-pIaBGL1$posteriorMedian
gaBGL1<-graph_from_adjacency_matrix(abs(medaBGL1),mode="undirected",diag=FALSE,weighted=TRUE)
E(gaBGL1)$width<-(E(gaBGL1)$weight**2)/4
medaBGL1[lower.tri(medaBGL1,diag=TRUE)]<-NA
na.omit(c(medaBGL1))
E(gaBGL1)$color<-c("darkred","navyblue")[as.integer(na.omit(c(t(medaBGL1)))>0)+1L]
plot(gaBGL1)
