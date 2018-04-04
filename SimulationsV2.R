########### Prereqs ###########
library(MASS)
library(tidyverse)
library(igraph)
library(clusterGeneration)
library(BayesianGLasso)
library(igraph)

options(stringsAsFactors=FALSE)
oldPar<-par()
setwd("~/gdrive/Dissertation/Aim2")

########### Lambda analysis ############
fun1<-function(r,s,dis,omega,x){
  return((1+r)/(s+omega+dis))
}
df1<-expand.grid(dis=seq(0.1,1,.2),omega=seq(0,1,.01))
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
nRV<-20L
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
pCors<-pCorFun(omega)
pCorsInd<-abs(pCors)>.1

# Simulated structural similarity:
sim<-toeplitz(topParamSim**(0:(nRV-1L))) # True similarity matrix

# Simulated mv random normal:
set.seed(33)
x1<-mvrnorm(n=nObs,mu=rep(0,ncol(sig)),Sigma=4*sig)

# Concentration matrix graph:
Om1<-omega
# Om1<-pCors
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

########### Simulation function ############
simGrid<-expand.grid(gammaPriorr=10**seq(-2,1.5,.5),gammaPriors=10**(seq(-3,0,by=.5)),
            adaptive=c(FALSE,TRUE),adaptiveType=c("norm","priorHyper"),
            stringsAsFactors = FALSE)
simGrid<-simGrid %>% filter(!(adaptive==FALSE & adaptiveType=="priorHyper"))
simGrid$iterations<-1000
simGrid$burnIn<-100
simGrid$f1<-simGrid$ppv<-simGrid$spec<-simGrid$sens<-NA
for(i in 1:nrow(simGrid)){
  # Gibbs sampler:
  bgl<-NULL
  attempt<-0
  while(is.null(bgl) & attempt<3){
    attempt<-attempt+1
    try(
      if(simGrid$adaptive[i]){
        bgl<-blockGLasso(x1,iterations=simGrid$iterations[i],burnIn=simGrid$burnIn[i],
                         adaptive=simGrid$adaptive[i],adaptiveType=simGrid$adaptiveType[i],
                         priorHyper=abs(solve(sim)),gammaPriorr=simGrid$gammaPriorr[i],
                         gammaPriors=simGrid$gammaPriors[i])
      }else{
        bgl<-blockGLasso(x1,iterations=simGrid$iterations[i],burnIn=simGrid$burnIn[i],
                         adaptive=simGrid$adaptive[i],lambdaPriora=simGrid$gammaPriorr[i],
                         lambdaPriorb=simGrid$gammaPriors[i])
      }
    )
  }
  if(attempt>=3 & is.null(bgl)){
    next
  }
  
  # Posterior inference object:
  bgl$Omegas<-lapply(bgl$Omegas,pCorFun)
  pIBgl<-posteriorInference(bgl)
  
  # Posterior median:
  medBgl<-pIBgl$posteriorMedian

  # Topological error analysis:
  pCorsMedBgl<-pCorFun(medBgl)
  pCorsIndMedBgl<-abs(pCorsMedBgl)>.2
  tabBgl<-xtabs(~true+pred,data=data.frame(true=c(pCorsInd),pred=c(pCorsIndMedBgl)))
  simGrid$sens[i]<-tabBgl['TRUE','TRUE']/sum(tabBgl['TRUE',])
  simGrid$spec[i]<-tabBgl['FALSE','FALSE']/sum(tabBgl['FALSE',])
  simGrid$ppv[i]<-tabBgl['TRUE','TRUE']/sum(tabBgl[,'TRUE'])
  simGrid$f1[i]<-(2*simGrid$sens[i]*simGrid$ppv[i])/(simGrid$sens[i]+simGrid$ppv[i])
}

########### Output graphs ###########
graphFun<-function(adaptive,adaptiveType,gammaPriorr,gammaPriors,
                   lambdaPriora=NULL,lambdaPriorb=NULL){
  # Sampler:
  bgl<-blockGLasso(x1,iterations=10000,burnIn=1000,adaptive=adaptive,
                   adaptiveType=adaptiveType,priorHyper=(abs(solve(sim))**2),
                   gammaPriorr=gammaPriorr,gammaPriors=gammaPriors,
                   lambdaPriora=lambdaPriora,lambdaPriorb=lambdaPriorb)
  pIBgl<-posteriorInference(bgl)
  bglMed<-pIBgl$posteriorMedian
  bglCor<-pCorFun(bglMed)
  
  # Graph
  bglG<-graph_from_adjacency_matrix(abs(bglMed),mode="undirected",diag=FALSE,weighted=TRUE)
  E(bglG)$width<-(E(bglG)$weight**2)
  bglMed[lower.tri(bglMed,diag=TRUE)]<-NA
  E(bglG)$color<-c("darkred","navyblue")[as.integer(na.omit(c(t(bglMed)))>0)+1L]
  return(bglG)
}
aIG<-graphFun(adaptive=TRUE,adaptiveType="priorHyper",gammaPriorr=1,gammaPriors=10)
aNG<-graphFun(adaptive=TRUE,adaptiveType="norm",gammaPriorr=1,gammaPriors=5)
rG<-graphFun(adaptive=FALSE,adaptiveType=NULL,gammaPriorr=NULL,gammaPriors=NULL,
             lambdaPriora=1,lambdaPriorb=10)

plot(aIG)
plot(aNG)
plot(rG)
