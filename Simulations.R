########### Prereqs ###########
library(MASS)
library(dplyr)
library(tidyr)
library(BayesianGLasso)
library(igraph)
library(ggplot2)

setwd("~/gdrive/Dissertation/Aim2")

########### Simulated data ###########
nRV<-15L
nObs<-10L
topParam<-.9
topParamSim<-.975

# AR(1) example
sig<-toeplitz(topParam**(0:(nRV-1L))) # True covariance matrix
omega<-solve(sig) # True concentration matrix

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
pCorFun(omega)

# Simulated structural similarity:
sim<-toeplitz(topParamSim**(0:(nRV-1L))) # True similarity matrix

# Simulated mv random normal:
set.seed(3)
x1<-mvrnorm(n=nObs,mu=rep(0,ncol(sig)),Sigma=sig)

########### Regular BGL Test ############
BGL1<-blockGLasso(x1,iterations=1000,burnIn=500,adaptive=FALSE,lambdaPriora=1,
            lambdaPriorb=1/10)

# Posterior inference object:
pIBGL1<-posteriorInference(BGL1)

# Posterior median:
medBGL1<-pIBGL1$posteriorMedian

# Signed absolute error:
errBGL1<-omega-medBGL1

# Graph:
gBGL1<-graph_from_adjacency_matrix(abs(medBGL1),mode="undirected",diag=FALSE,weighted=TRUE)
E(gBGL1)$width<-(E(gBGL1)$weight**2)/4
medBGL1[lower.tri(medBGL1,diag=TRUE)]<-NA
na.omit(c(medBGL1))
E(gBGL1)$color<-c("darkred","navyblue")[as.integer(na.omit(c(t(medBGL1)))>0)+1L]
plot(gBGL1)

aBGL1<-blockGLasso(x1,iterations=1000,burnIn=500,adaptive=TRUE,adaptiveType="norm",
                   gammaPriors=10**(-2),gammaPriort=10**(-1))
pIaBGL1<-posteriorInference(aBGL1)
medaBGL1<-pIaBGL1$posteriorMedian
gaBGL1<-graph_from_adjacency_matrix(abs(medaBGL1),mode="undirected",diag=FALSE,weighted=TRUE)
E(gaBGL1)$width<-(E(gaBGL1)$weight**2)/4
medaBGL1[lower.tri(medaBGL1,diag=TRUE)]<-NA
na.omit(c(medaBGL1))
E(gaBGL1)$color<-c("darkred","navyblue")[as.integer(na.omit(c(t(medaBGL1)))>0)+1L]
plot(gaBGL1)

aBGL2<-blockGLasso(x1,iterations=1000,burnIn=500,adaptive=TRUE,adaptiveType="norm",
                   gammaPriors=10**(-2),gammaPriort=10**(1))
pIaBGL2<-posteriorInference(aBGL2)
medaBGL2<-pIaBGL2$posteriorMedian
gaBGL2<-graph_from_adjacency_matrix(abs(medaBGL2),mode="undirected",diag=FALSE,weighted=TRUE)
E(gaBGL2)$width<-(E(gaBGL2)$weight**2)/4
medaBGL2[lower.tri(medaBGL2,diag=TRUE)]<-NA
na.omit(c(medaBGL2))
E(gaBGL2)$color<-c("darkred","navyblue")[as.integer(na.omit(c(t(medaBGL2)))>0)+1L]
plot(gaBGL2)

df1<-rbind(data.frame(pen="large",value=c(medaBGL2)),
           data.frame(pen="small",value=c(medaBGL1)))
ggplot(df1,aes(x=pen,y=value,fill=pen))+geom_boxplot()+theme_bw()

########### Informative adaptive simulations ############
