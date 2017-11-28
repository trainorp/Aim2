########### Prereqs ###########
options(stringsAsFactors=FALSE)
oldPar<-par()

library(MASS)
library(dplyr)
library(tidyr)
library(igraph)
library(ggplot2)

source("~/gdrive/Dissertation/Aim2/BayesianGLasso/R/blockGLasso.default.R")
source("~/gdrive/Dissertation/Aim2/BayesianGLasso/R/blockAdGLasso.default.R")

setwd("~/gdrive/Dissertation/Aim2")

########### Lambda analysis ############
s<-.1

fun1<-function(dis,omega,x){
  return((1+dis)/(s+omega))
}
df1<-expand.grid(dis=seq(0,1,.2),omega=seq(0,1,.01))
df1$lambda<-NA

for(i in 1:nrow(df1)){
  df1$lambda[i]<-fun1(dis=df1$dis[i],omega=df1$omega[i])
}

df1$dis<-factor(df1$dis)

# png(filename="lambdasVsSim.png",width=4.5,height=3,units="in",res=600)
ggplot(data=df1,aes(x=omega,y=lambda,color=dis,group=dis))+geom_line()+theme_bw()+
  xlab(expression(paste("|",tilde(omega)[ij],"|")))+
  ylab(expression(paste(E(lambda["ij"]))))+
  scale_color_discrete(name="Tanimoto\nDissimilarity")
# dev.off()

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

# Concentration matrix graph:
Om1<-omega
gOm1<-graph_from_adjacency_matrix(abs(Om1),mode="undirected",diag=FALSE,weighted=TRUE)
E(gOm1)$width<-(E(gOm1)$weight**2)/4
Om1[lower.tri(Om1,diag=TRUE)]<-NA
E(gOm1)$color<-c("darkred","navyblue")[as.integer(na.omit(c(t(Om1)))>0)+1L]

# png(file="Om1.png",height=5,width=5,units="in",res=300)
par(mar=c(1,1,1,1))
set.seed(2)
plot(gOm1)
# dev.off()

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
E(gBGL1)$color<-c("darkred","navyblue")[as.integer(na.omit(c(t(medBGL1)))>0)+1L]
plot(gBGL1)

# Export plot
# png(file="gBGL1.png",height=5,width=5,units="in",res=300)
par(mar=c(1,1,1,1))
set.seed(3)
plot(gBGL1)
# dev.off()

# Adaptive with small gamma t
aBGL1<-blockGLasso(x1,iterations=1000,burnIn=500,adaptive=TRUE,adaptiveType="norm",
                   gammaPriors=10**(-2),gammaPriort=10**(-1))
pIaBGL1<-posteriorInference(aBGL1)
medaBGL1<-pIaBGL1$posteriorMedian
gaBGL1<-graph_from_adjacency_matrix(abs(medaBGL1),mode="undirected",diag=FALSE,weighted=TRUE)
E(gaBGL1)$width<-(E(gaBGL1)$weight**2)/4
medaBGL1[lower.tri(medaBGL1,diag=TRUE)]<-NA
E(gaBGL1)$color<-c("darkred","navyblue")[as.integer(na.omit(c(t(medaBGL1)))>0)+1L]
plot(gaBGL1)

# Adaptive with large gamma t
aBGL2<-blockGLasso(x1,iterations=1000,burnIn=500,adaptive=TRUE,adaptiveType="norm",
                   gammaPriors=10**(-2),gammaPriort=10**(1))
pIaBGL2<-posteriorInference(aBGL2)
medaBGL2<-pIaBGL2$posteriorMedian
gaBGL2<-graph_from_adjacency_matrix(abs(medaBGL2),mode="undirected",diag=FALSE,weighted=TRUE)
E(gaBGL2)$width<-(E(gaBGL2)$weight**2)/4
medaBGL2[lower.tri(medaBGL2,diag=TRUE)]<-NA
E(gaBGL2)$color<-c("darkred","navyblue")[as.integer(na.omit(c(t(medaBGL2)))>0)+1L]
plot(gaBGL2)

df1<-rbind(data.frame(pen="large",value=c(medaBGL2)),
           data.frame(pen="small",value=c(medaBGL1)))
ggplot(df1,aes(x=pen,y=value,fill=pen))+geom_boxplot()+theme_bw()

########### Informative adaptive simulations ############
aiBGL1<-blockGLasso(x1,iterations=1000,burnIn=500,adaptive=TRUE,adaptiveType="priorHyper",
                   priorHyper=(sim)**30,gammaPriors=10**(-2),gammaPriort=10**(-1))

pIaiBGL1<-posteriorInference(aiBGL1)
medaiBGL1<-pIaiBGL1$posteriorMedian
gaiBGL1<-graph_from_adjacency_matrix(abs(medaiBGL1),mode="undirected",diag=FALSE,weighted=TRUE)
E(gaiBGL1)$width<-(E(gaiBGL1)$weight**2)/4
medaiBGL1[lower.tri(medaiBGL1,diag=TRUE)]<-NA
E(gaiBGL1)$color<-c("darkred","navyblue")[as.integer(na.omit(c(t(medaiBGL1)))>0)+1L]

# png(file="gaiBGL1.png",height=5,width=5,units="in",res=300)
par(mar=c(1,1,1,1))
set.seed(3)
plot(gaiBGL1)
# dev.off()

# png(file="OmAll.png",height=4,width=12,units="in",res=300)
par(oma=c(1,1,1,1),mar=c(0,0,0,0),mfrow=c(1,3))
set.seed(2)
plot(gOm1)
set.seed(3)
plot(gBGL1)
set.seed(3)
plot(gaiBGL1)
# dev.off()
