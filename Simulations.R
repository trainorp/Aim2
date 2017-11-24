########### Prereqs ###########
library(MASS)
library(dplyr)
library(tidyr)
setwd("~/gdrive/Dissertation/Aim2")

# AR(1) example
sig<-toeplitz(.9**(0:9)) #toeplitz(.9**(0:109))
#sig[1,10]<-.9
omega<-solve(sig)
pcors<-matrix(0,nrow=nrow(omega),ncol=ncol(omega))
for(j in 1:ncol(pcors))
{
  for(i in 1:nrow(pcors))
  {
    pcors[i,j]<-(-omega[i,j]/sqrt(omega[i,i]*omega[j,j]))
  }
}
x1<-mvrnorm(n=50,mu=rep(0,ncol(sig)),Sigma=sig)

########### Regular BGL Test ############
X=x1;iterations=2000;burnIn=1000;
lambdaPriora=1;lambdaPriorb=1/10;
illStart=c("identity","glasso");rho=.1;
verbose=TRUE
idk1<-blockGLasso.default(X=x1,iterations=100,burnIn = 0)
idk2<-blockAdGLasso.default(X=x1,iterations=1)
plot(1:1001,sapply(idk2$Omegas,function(x) x[1,1]),type="l")

########### Adaptive BGL Test ############
# Simulated a priori knowledge of ar(1) structure:
prior<-1e-6/toeplitz(rev(1:ncol(sig)))
prior<-1e-6/toeplitz(1:ncol(sig))
prior[1,10]<-1/2
idk<-blockAdGLasso.default(X=x1,iterations=500,burnIn=1,adaptiveType="priorHyper",
                           priorHyper=prior,gammaPriors=1e-2,gammaPriort=1e-6)
idk2<-blockAdGLasso.default(X=x1,iterations=500,burnIn=1,gammaPriors=1e-2,gammaPriort=1e-6)

plot(1:501,sapply(idk$Omegas,function(x) x[1,2]),type="l")
points(1:501,sapply(idk2$Omegas,function(x) x[1,2]),type="l",col="red")
summary(sapply(idk$Omegas,function(x) x[1,2]))
summary(sapply(idk2$Omegas,function(x) x[1,2]))
plot(1:501,sapply(idk$Omegas,function(x) x[1,3]),type="l")
points(1:501,sapply(idk2$Omegas,function(x) x[1,3]),type="l",col="red")
summary(sapply(idk$Omegas,function(x) x[1,3]))
summary(sapply(idk2$Omegas,function(x) x[1,3]))
plot(1:501,sapply(idk$Omegas,function(x) x[2,3]),type="l")
points(1:501,sapply(idk2$Omegas,function(x) x[2,3]),type="l",col="red")
summary(sapply(idk$Omegas,function(x) x[2,3]))
summary(sapply(idk2$Omegas,function(x) x[2,3]))
plot(1:501,sapply(idk$Omegas,function(x) x[1,10]),type="l")
points(1:501,sapply(idk2$Omegas,function(x) x[1,10]),type="l",col="red")
summary(sapply(idk$Omegas,function(x) x[1,10]))
summary(sapply(idk2$Omegas,function(x) x[1,10]))

sapply(idk$Sigmas,function(x) x[1,2])
posteriorInference(idk)

set.seed(3)
s1<-genPositiveDefMat(11,covMethod="c-vine")$Sigma
#write.table(s1,file="cov1.csv",row.names=FALSE,col.names=FALSE,sep=",")
set.seed(33)
s2<-.99**toeplitz(0:9)
#write.table(s2,file="cov2.csv",row.names=FALSE,col.names=FALSE,sep=",")

set.seed(4)
x1<-mvrnorm(n=10,mu=rep(0,11),Sigma=s1)
#write.table(x1,file="x1.csv",row.names=FALSE,col.names=FALSE,sep=",")
set.seed(5)
x2<-mvrnorm(n=10,mu=rep(0,10),Sigma=s2)
#write.table(x2,file="x2.csv",row.names=FALSE,col.names=FALSE,sep=",")

omegaT1<-solve(s1)
omegaT2<-solve(s2)

save(x1,file="x1.RData")
load(file="x1.RData")

source('~/gdrive/Dissertation/Aim2/BayesianGLasso/R/blockGLasso.R')
idk<-blockGLasso(X=x1)
