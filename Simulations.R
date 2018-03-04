########### Prereqs ###########
library(MASS)
library(tidyverse)
library(igraph)
library(clusterGeneration)
library(BayesianGLasso)

options(stringsAsFactors=FALSE)
oldPar<-par()
setwd("~/gdrive/Dissertation/Aim2")

########### Lambda analysis ############
fun1<-function(r,s,dis,omega,x){
  return((1+r)/(s+omega+dis))
}
df1<-expand.grid(dis=seq(0.1,1,.1),omega=seq(0,1,.01))
df1$lambda<-NA

for(i in 1:nrow(df1)){
  df1$lambda[i]<-fun1(r=1,s=.1,dis=df1$dis[i],omega=df1$omega[i])
}

df1$dis<-factor(df1$dis)

png(filename="Plots/lambdasVsSim.png",width=4.5,height=3.5,units="in",res=600)
ggplot(data=df1,aes(x=omega,y=lambda,color=dis,group=dis))+geom_line()+theme_bw()+
  xlab(expression(paste("|",tilde(omega)[ij],"|")))+
  ylab(expression(paste(E(lambda["ij"]))))+
  scale_color_discrete(name="Similarity")
dev.off()

########### Simulated AR(1) data ###########
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

png(file="Plots/AR1_Om.png",height=5,width=5,units="in",res=300)
par(mar=c(1,1,1,1))
set.seed(2)
plot(gOm1)
dev.off()

# Correlation adjacency plot:
gCor<-graph_from_adjacency_matrix(cor(x1),mode="undirected",diag=FALSE,weighted=TRUE)
E(gCor)$width<-(E(gCor)$weight**2)
E(gCor)$color<-"darkred"

png(file="Plots/AR1_Cor.png",height=5,width=5,units="in",res=300)
par(mar=c(1,1,1,1))
set.seed(32)
plot(gCor)
dev.off()

########### Simulate Random multivariate Gaussian (c-vines) ############
set.seed(33)
sigR<-genPositiveDefMat(15,covMethod="c-vine")$Sigma
omegaR<-solve(sigR)
pCorFun(omegaR)

# Plot as graph:
Om1R<-omegaR
gOm1R<-graph_from_adjacency_matrix(abs(Om1R),mode="undirected",diag=FALSE,weighted=TRUE)
E(gOm1R)$width<-(E(gOm1R)$weight**2)/4
Om1R[lower.tri(Om1R,diag=TRUE)]<-NA
E(gOm1R)$color<-c("darkred","navyblue")[as.integer(na.omit(c(t(Om1R)))>0)+1L]
set.seed(333)
plot(gOm1R)

########### Regular BGL Test (AR1) ############
BGLres<-data.frame()
BGLgrid<-expand.grid(gammaPriorr=c(1,2,4,8,16),gammaPriors=10**(seq(-2,2,by=1)))
iterations<-10000
burnIn<-1000
for(i in 1:nrow(BGLgrid)){
  BGL1<-blockGLasso(x1,iterations=iterations,burnIn=burnIn,adaptive=FALSE,
          lambdaPriora=BGLgrid$gammaPriorr[i],lambdaPriorb=BGLgrid$gammaPriors[i])
  
  # Posterior inference object:
  pIBGL1<-posteriorInference(BGL1)
  
  # Posterior median:
  medBGL1<-pIBGL1$posteriorMedian
  medBGL1Sigma<-solve(medBGL1)
  
  BGL1Errs<-data.frame(var="err",
                       val=sapply(BGL1$Omegas,function(x) mean(abs(omega-x)))[(burnIn+1):(burnIn+iterations)])
  BGL1Lambdas<-data.frame(var="lambdas",val=BGL1$lambdas[(burnIn+1):(burnIn+iterations)])
  BGL1res<-rbind(BGL1Errs,BGL1Lambdas)
  BGL1res$gammaPriorr<-BGLgrid$gammaPriorr[i]
  BGL1res$gammaPriors<-BGLgrid$gammaPriors[i]
  BGLres<-rbind(BGLres,BGL1res)
}

# Boxplot of error as a function of lambda priors
BGLres$gammaPriorr<-factor(BGLres$gammaPriorr)
BGLres$gammaPriors<-factor(BGLres$gammaPriors)
BGLres$Var<-factor(BGLres$var,levels=c("lambdas","err"),labels=c("lambda","`Mean Absolute Error`"))

# Boxplots of lambdas and Mean Absolute Error
png(filename="Plots/AR1_LambdaError.png",height=4,width=5.5,units="in",res=300)
ggplot(BGLres,aes(x=gammaPriors,y=val,fill=gammaPriorr))+
  geom_boxplot()+xlab(expression(paste("Gamma hyperparameter ",italic(s))))+ylab("Value")+
  scale_fill_discrete(name=expression(italic(r)))+theme_bw()+
  facet_wrap(~Var,nrow=2,scale="free_y",labeller=label_parsed)
dev.off()

# Graph:
BGL1<-blockGLasso(x1,iterations=iterations,burnIn=burnIn,adaptive=FALSE,
                  lambdaPriora=16,lambdaPriorb=.01)
pIBGL1<-posteriorInference(BGL1)
medBGL1<-pIBGL1$posteriorMedian
medBGL1Sigma<-solve(medBGL1)

gBGL1<-graph_from_adjacency_matrix(abs(medBGL1),mode="undirected",diag=FALSE,weighted=TRUE)
E(gBGL1)$width<-(E(gBGL1)$weight**2)
medBGL1[lower.tri(medBGL1,diag=TRUE)]<-NA
E(gBGL1)$color<-c("darkred","navyblue")[as.integer(na.omit(c(t(medBGL1)))>0)+1L]
plot(gBGL1)

# Export plot
png(file="Plots/gBGL1.png",height=5,width=5,units="in",res=300)
par(mar=c(1,1,1,1))
set.seed(3)
plot(gBGL1)
dev.off()

########### Adaptive BGL Test (AR1) Noninformative ############
aBGLres<-data.frame()
BGLgrid<-expand.grid(gammaPriorr=10**seq(-2,1.5,.5),gammaPriors=10**(seq(-3,0,by=.5)))
iterations<-10000
burnIn<-1000
for(i in 1:nrow(BGLgrid)){
  aBGL1<-blockGLasso(x1,iterations=iterations,burnIn=burnIn,adaptive=TRUE,adaptiveType="norm",
        gammaPriorr=BGLgrid$gammaPriorr[i],gammaPriors=BGLgrid$gammaPriors[i])
  
  # Posterior inference object:
  pIaBGL1<-posteriorInference(aBGL1)
  
  # Posterior median:
  medaBGL1<-pIaBGL1$posteriorMedian
  medaBGL1Sigma<-solve(medaBGL1)
  
  aBGL1Errs<-data.frame(var="err",
                       val=sapply(aBGL1$Omegas,function(x) mean(abs(omega-x)))[(burnIn+1):(burnIn+iterations)])
  lambdaMatList<-aBGL1$lambdas[(burnIn+1):(burnIn+iterations)]
  aBGL1Lambdas<-data.frame(var="lambdas",val=sapply(lambdaMatList,function(x) median(x[x>0])))
  aBGL1res<-rbind(aBGL1Errs,aBGL1Lambdas)
  aBGL1res$gammaPriorr<-BGLgrid$gammaPriorr[i]
  aBGL1res$gammaPriors<-BGLgrid$gammaPriors[i]
  aBGLres<-rbind(aBGLres,aBGL1res)
}

# Boxplot of error as a function of priors
aBGLres$gammaPriorr<-factor(aBGLres$gammaPriorr)
aBGLres$gammaPriors<-factor(aBGLres$gammaPriors)
ggplot(aBGLres %>% filter(var=="err"),aes(x=gammaPriorr,y=val,fill=gammaPriors))+
  geom_boxplot()+ylab("Error")+theme_bw()

# Graph:
aBGL1<-blockGLasso(x1,iterations=iterations,burnIn=burnIn,adaptive=TRUE,adaptiveType="norm",
                   gammaPriorr=.316,gammaPriors=.0316)

pIaBGL1<-posteriorInference(aBGL1)
medaBGL1<-pIaBGL1$posteriorMedian
medaBGL1Sigma<-solve(medaBGL1)

gaBGL1<-graph_from_adjacency_matrix(abs(medaBGL1),mode="undirected",diag=FALSE,weighted=TRUE)
E(gaBGL1)$width<-(E(gaBGL1)$weight**2)
medaBGL1[lower.tri(medaBGL1,diag=TRUE)]<-NA
E(gaBGL1)$color<-c("darkred","navyblue")[as.integer(na.omit(c(t(medaBGL1)))>0)+1L]
plot(gaBGL1)

########### Adaptive BGL Test (AR1) Informative ############
BGLres<-data.frame()
BGLgrid<-expand.grid(gammaPriorr=10**seq(-2,1.5,.5),gammaPriors=10**(seq(-3,0,by=.5)))
iterations<-1000
burnIn<-100
priorHyper<-abs(solve(sim))+.1
for(i in 1:nrow(BGLgrid)){
  aiBGL1<-blockGLasso(x1,iterations=iterations,burnIn=burnIn,adaptive=TRUE,
                     adaptiveType="priorHyper",priorHyper=abs(solve(sim)),
                     gammaPriorr=BGLgrid$gammaPriorr[i],gammaPriors=BGLgrid$gammaPriors[i])
  
  # Posterior inference object:
  pIaiBGL1<-posteriorInference(aiBGL1)
  
  # Posterior median:
  medaiBGL1<-pIaiBGL1$posteriorMedian
  medaiBGL1Sigma<-solve(medaiBGL1)
  
  aiBGL1Errs<-data.frame(var="err",
                        val=sapply(aiBGL1$Omegas,function(x) mean(abs(omega-x)))[(burnIn+1):(burnIn+iterations)])
  lambdaMatList<-aiBGL1$lambdas[(burnIn+1):(burnIn+iterations)]
  aiBGL1Lambdas<-data.frame(var="lambdas",val=sapply(lambdaMatList,function(x) median(x[x>0])))
  aiBGL1res<-rbind(aiBGL1Errs,aiBGL1Lambdas)
  aiBGL1res$gammaPriorr<-BGLgrid$gammaPriorr[i]
  aiBGL1res$gammaPriors<-BGLgrid$gammaPriors[i]
  BGLres<-rbind(BGLres,aiBGL1res)
}

# Boxplot of error as a function of priors
BGLres$gammaPriorr<-factor(BGLres$gammaPriorr)
BGLres$gammaPriors<-factor(BGLres$gammaPriors)
ggplot(BGLres %>% filter(var=="err"),aes(x=gammaPriorr,y=val,fill=gammaPriors))+
  geom_boxplot()+ylab("Error")+theme_bw()

# Lambda as a function of priors:
ggplot(BGLres %>% filter(var=="lambdas"),aes(x=gammaPriorr,y=log(val),fill=gammaPriors))+
  geom_boxplot()+ylab("Lambda")+theme_bw()

# Graph:
aiBGL1<-blockGLasso(x1,iterations=iterations,burnIn=burnIn,adaptive=TRUE,
                   adaptiveType="priorHyper",priorHyper=priorHyper,
                   gammaPriorr=1,gammaPriors=.001)

pIaiBGL1<-posteriorInference(aiBGL1)
medaiBGL1<-pIaiBGL1$posteriorMedian
medaiBGL1Sigma<-solve(medaiBGL1)

gaiBGL1<-graph_from_adjacency_matrix(abs(medaiBGL1),mode="undirected",diag=FALSE,weighted=TRUE)
E(gaiBGL1)$width<-(E(gaiBGL1)$weight**2)
medaiBGL1[lower.tri(medaiBGL1,diag=TRUE)]<-NA
E(gaiBGL1)$color<-c("darkred","navyblue")[as.integer(na.omit(c(t(medaiBGL1)))>0)+1L]
plot(gaiBGL1)
