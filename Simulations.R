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

png(filename="lambdasVsSim.png",width=4.5,height=3.5,units="in",res=600)
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

png(file="Om1.png",height=5,width=5,units="in",res=300)
par(mar=c(1,1,1,1))
set.seed(2)
plot(gOm1)
dev.off()

# Correlation adjacency plot:
gCor<-graph_from_adjacency_matrix(cor(x1),mode="undirected",diag=FALSE,weighted=TRUE)
E(gCor)$width<-(E(gCor)$weight**2)
E(gCor)$color<-"darkred"

png(file="gCor.png",height=5,width=5,units="in",res=300)
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
set.seed(36)
plot(gOm1R)

########### Regular BGL Test (AR1) ############
BGLres<-data.frame()
BGLgrid<-expand.grid(lambdaPriora=c(1,2,4,8,16),lambdaPriorb=10**(seq(-1,2,by=.5)))
for(i in 1:nrow(BGLgrid)){
  iterations<-10000
  burnIn<-1000
  lambdaPriora<-BGLgrid$lambdaPriora[i]
  lambdaPriorb<-BGLgrid$lambdaPriorb[i]
  BGL1<-blockGLasso(x1,iterations=iterations,burnIn=burnIn,adaptive=FALSE,
                    lambdaPriora=lambdaPriora,lambdaPriorb=lambdaPriorb)
  
  # Posterior inference object:
  pIBGL1<-posteriorInference(BGL1)
  
  # Posterior median:
  medBGL1<-pIBGL1$posteriorMedian
  medBGL1Sigma<-solve(medBGL1)
  
  BGL1Errs<-data.frame(var="err",
                       val=sapply(BGL1$Omegas,function(x) mean(abs(omega-x)))[(burnIn+1):(burnIn+iterations)])
  BGL1Lambdas<-data.frame(var="lambdas",val=BGL1$lambdas[(burnIn+1):(burnIn+iterations)])
  BGL1res<-rbind(BGL1Errs,BGL1Lambdas)
  BGL1res$lambdaPriora<-lambdaPriora
  BGL1res$lambdaPriorb<-lambdaPriorb
  BGLres<-rbind(BGLres,BGL1res)
}

# Boxplot of error as a function of lambda priors
BGLres$lambdaPriora<-factor(BGLres$lambdaPriora)
BGLres$lambdaPriorb<-factor(BGLres$lambdaPriorb)
ggplot(BGLres %>% filter(var=="err"),aes(x=lambdaPriorb,y=val,fill=lambdaPriora))+
  geom_boxplot()+ylab("Error")+theme_bw()

# Lambda as a function of priors:
ggplot(BGLres %>% filter(var=="lambdas"),aes(x=lambdaPriorb,y=val,fill=lambdaPriora))+
  geom_boxplot()+ylab("Lambda")+theme_bw()

# Lambda versus error:
BGLresSum<-BGLres %>% group_by(var,lambdaPriora,lambdaPriorb) %>% summarise(val=median(val))
BGLresSum<-BGLresSum %>% spread(key=var,value=val)
ggplot(BGLresSum,aes(x=lambdas,y=err))+geom_point()+theme_bw()

# Graph:
gBGL1<-graph_from_adjacency_matrix(abs(medBGL1),mode="undirected",diag=FALSE,weighted=TRUE)
E(gBGL1)$width<-(E(gBGL1)$weight**2)/4
medBGL1[lower.tri(medBGL1,diag=TRUE)]<-NA
E(gBGL1)$color<-c("darkred","navyblue")[as.integer(na.omit(c(t(medBGL1)))>0)+1L]
plot(gBGL1)

# Export plot
png(file="gBGL1.png",height=5,width=5,units="in",res=300)
par(mar=c(1,1,1,1))
set.seed(3)
plot(gBGL1)
dev.off()

########### Regular BGL Test (D-vines) ############

########### Adaptive BGL Test (AR1) ############
# Adaptive with small gamma t
aBGL1<-blockGLasso(x1,iterations=1000,burnIn=500,adaptive=TRUE,adaptiveType="norm",
                   gammaPriors=10**(-2),gammaPriort=10**(-1))
pIaBGL1<-posteriorInference(aBGL1)
medaBGL1<-pIaBGL1$posteriorMedian
hist(sapply(aBGL1$Omegas,function(x) mean(abs(omega-x))))

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

########### Informative adaptive simulations ############
aiBGL1<-blockGLasso(x1,iterations=1000,burnIn=500,adaptive=TRUE,adaptiveType="priorHyper",
                   priorHyper=20*(sim)**10,gammaPriors=10**(1),gammaPriort=10**(-1))

# Analysis of distribution for lambda
exOmegas<-aiBGL1$Omegas[[999]]
exLambdas<-aiBGL1$lambdas[[999]]
plot(c(exLambdas)~c(sim))
mean(exLambdas[exLambdas>0])
abline(lm(c(exLambdas)~c(sim)))

# Make plot:
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
