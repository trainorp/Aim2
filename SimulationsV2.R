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

png(filename="Plots/lambdasVsSim.png",width=5.5,height=3.5,units="in",res=600)
ggplot(data=df1,aes(x=omega,y=lambda,color=dis,group=dis))+geom_line()+theme_bw()+
  xlab(expression(paste("|",tilde(omega)[ij],"|")))+
  ylab(expression(paste(E(lambda["ij"]))))+
  scale_color_discrete(name="Similarity")
dev.off()

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

########### One Round Simulated AR(1) data ###########
nRV<-20L
nObs<-10L
topParam<-.9
topParamSim<-.975

# AR(1) example
sig<-toeplitz(topParam**(0:(nRV-1L))) # True covariance matrix
omega<-solve(sig) # True concentration matrix

# True partial correlations:
pCors<-pCorFun(omega)
pCorsInd<-abs(pCors)>.1

# Simulated structural similarity:
sim<-toeplitz(topParamSim**(0:(nRV-1L))) # True similarity matrix

# Simulated mv random normal:
set.seed(333)
x1<-mvrnorm(n=nObs,mu=rep(0,ncol(sig)),Sigma=4*sig)

# Concentration matrix graph:
Om1<-omega
# Om1<-pCors
gOm1<-graph_from_adjacency_matrix(abs(Om1),mode="undirected",diag=FALSE,weighted=TRUE)
E(gOm1)$width<-(E(gOm1)$weight**2)/4
Om1[lower.tri(Om1,diag=TRUE)]<-NA
E(gOm1)$color<-c("darkred","navyblue")[as.integer(na.omit(c(t(Om1)))>0)+1L]

png(file="Plots/AR1_Om.png",height=5,width=5,units="in",res=600)
par(bg=NA,mar=c(0,0,0,0))
set.seed(3)
plot(gOm1)
dev.off()

# Correlation adjacency plot:
gCor<-graph_from_adjacency_matrix(cor(x1),mode="undirected",diag=FALSE,weighted=TRUE)
E(gCor)$width<-(E(gCor)$weight**2)
E(gCor)$color<-"darkred"

png(file="Plots/AR1_Cor.png",height=5,width=5,units="in",res=300)
par(bg=NA,mar=c(0,0,0,0))
set.seed(6)
plot(gCor)
dev.off()

########### Simulation function ############
simGridFun<-function(x1,simGrid){
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
                           priorHyper=eval(parse(text=simGrid$priorHyper[i])),
                           gammaPriors=simGrid$gammaPriorr[i],gammaPriort=simGrid$gammaPriors[i])
        }else{
          bgl<-blockGLasso(x1,iterations=simGrid$iterations[i],burnIn=simGrid$burnIn[i],
                           adaptive=simGrid$adaptive[i],lambdaPriora=simGrid$gammaPriorr[i],
                           lambdaPriorb=simGrid$gammaPriors[i])
        }
      )
    }
    if(attempt>=3 & is.null(bgl)) next
    
    # Likelihood:
    sigmas<-lapply(bgl$Omegas[(bgl$burnIn+1:length(bgl$Omegas)-bgl$burnIn)],solve)
    simGrid$ll[i]<-mean(sapply(sigmas,function(x) 
      -sum(log(mvtnorm::dmvnorm(x=x1,mean=rep(0,ncol(x1)),sigma=x)))))
    
    # Posterior inference object:
    bgl$Omegas<-lapply(bgl$Omegas,pCorFun)
    pIBgl<-posteriorInference(bgl)
    
    # Posterior median:
    medBgl<-pIBgl$posteriorMedian
    
    # Average lambdas:
    if(simGrid$adaptive[i]){
      simGrid$medLambda[i]<-median(sapply(bgl$lambdas,function(x) median(x[upper.tri(x)])))
      simGrid$meanLambda[i]<-mean(sapply(bgl$lambdas,function(x) mean(x[upper.tri(x)])))
    }else{
      simGrid$medLambda[i]<-median(bgl$lambdas)
      simGrid$meanLambda[i]<-mean(bgl$lambdas)
    }
    
    # Topological error analysis:
    pCorsMedBgl<-pCorFun(medBgl)
    pCorsIndMedBgl<-abs(pCorsMedBgl)>.2
    tabBgl<-xtabs(~true+pred,data=data.frame(true=c(pCorsInd),pred=c(pCorsIndMedBgl)))
    simGrid$sens[i]<-tabBgl['TRUE','TRUE']/sum(tabBgl['TRUE',])
    simGrid$spec[i]<-tabBgl['FALSE','FALSE']/sum(tabBgl['FALSE',])
    simGrid$ppv[i]<-tabBgl['TRUE','TRUE']/sum(tabBgl[,'TRUE'])
    simGrid$f1[i]<-(2*simGrid$sens[i]*simGrid$ppv[i])/(simGrid$sens[i]+simGrid$ppv[i])
    
    simGrid$auc[i]<-as.numeric(pROC::roc(response=c(pCorsInd),predictor=c(abs(pCorsMedBgl)))$auc)
  }
  return(simGrid)
}

# Simulation Grid:
pHGood<-abs(solve(sim))
simGrid<-expand.grid(gammaPriorr=10**seq(-2,1.5,1),gammaPriors=10**(seq(-3,0,by=.5)),
                     adaptive=c(FALSE,TRUE),adaptiveType=c("norm","priorHyper"),
                     priorHyper="pHGood",stringsAsFactors = FALSE)
simGrid<-simGrid %>% filter(!(adaptive==FALSE & adaptiveType=="priorHyper"))
simGrid$iterations<-1000
simGrid$burnIn<-100
simGrid$medLambda<-simGrid$meanLambda<-simGrid$ll<-simGrid$auc<-simGrid$f1<-
  simGrid$ppv<-simGrid$spec<-simGrid$sens<-NA

# One iteration:
ptm<-proc.time()
simGridBig<-simGridFun(x1,simGrid)
proc.time()-ptm

# For hyperparameter optimization:
for(j in 1:20){
  set.seed(j)
  x2<-mvrnorm(n=nObs,mu=rep(0,ncol(sig)),Sigma=4*sig)
  simGridBig<-rbind(simGridBig,simGridFun(x2,simGrid))
}
simGridBigSum<-simGridBig %>% group_by(gammaPriorr,gammaPriors,adaptive,adaptiveType) %>% 
  summarize(auc=mean(auc)) %>% arrange(adaptive,adaptiveType,desc(auc))

# Hyperparameter optimization result:
# r= 0.10 s=0.01 Regular bgl
# r=1.0 s=0.01 structure adaptive
# r=1.0 s=0.0316 norm adaptive

simGrid<-data.frame(gammaPriorr=c(.1,1,1,1),gammaPriors=c(0.01,0.01,0.031622777,0.031622777),
                    adaptive=c(FALSE,TRUE,TRUE,TRUE),adaptiveType=c("norm","norm","priorHyper","priorHyper"),
                    priorHyper=c("pHGood","pHGood","pHGood","pHBad"),iterations=1000,burnIn=100)
simGrid$medLambda<-simGrid$meanLambda<-simGrid$auc<-simGrid$f1<-
  simGrid$ppv<-simGrid$spec<-simGrid$sens<-NA
pHGood<-abs(solve(sim))

# Simulations for performance analysis:
library(doParallel)
cl<-makeCluster(4)
registerDoParallel(cl)
ptm<-proc.time()
simGridBig<-foreach(j=1:2500,.combine="rbind",.packages=c("clusterGeneration","BayesianGLasso","tidyverse","pROC"),
             .export=ls(.GlobalEnv),.errorhandling="pass",.inorder=FALSE) %dopar% {
  set.seed(j+333)
  pHBad<-abs(solve(sim+matrix(rnorm(n=nrow(sim)*ncol(sim),mean=0,sd=.025),nrow=nrow(sim),
                              ncol=ncol(sim))))
  x2<-mvrnorm(n=nObs,mu=rep(0,ncol(sig)),Sigma=4*sig)
  simGridFun(x2,simGrid)
}
proc.time()-ptm
stopCluster(cl)
save(simGridBig,file="simGridBig.RData")
load(file="simGridBig.RData")

########### Simulation analysis ###########
simGridBig<-simGridBig %>% select(-gammaPriorr,-gammaPriors,-iterations,-burnIn)
simGridBig$Technique<-"BGL"
simGridBig$Technique[simGridBig$adaptive & simGridBig$adaptiveType=="norm"]<-"Adaptive BGL"
simGridBig$Technique[simGridBig$adaptive & simGridBig$adaptiveType=="priorHyper" & 
                       simGridBig$priorHyper=="pHGood"]<-"Chem. Structure\nAdaptive BGL (Good prior)"
simGridBig$Technique[simGridBig$adaptive & simGridBig$adaptiveType=="priorHyper" & 
                       simGridBig$priorHyper=="pHBad"]<-"Chem. Structure\nAdaptive BGL (Poor prior)"
simGridBig$Technique<-factor(simGridBig$Technique)

########### Performance plots ###########
png(file="Plots/AR1AUC.png",height=4,width=7,units="in",res=500,bg = "transparent")
ggplot(simGridBig,aes(x=auc,y=..density..,fill=Technique))+
  geom_histogram(binwidth=.005,alpha=.75,position="identity")+
  theme_bw()+xlab("AUC")+ylab("Density")+
  guides(fill=guide_legend(keywidth=1,keyheight=1.5,default.unit="line"))+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))
dev.off()

png(file="Plots/AR1F1.png",height=4,width=7,units="in",res=500,bg = "transparent")
ggplot(simGridBig,aes(x=f1,y=..density..,fill=Technique))+
  geom_histogram(binwidth=.02,alpha=.75,position="identity")+
  theme_bw()+xlab("F1 Measure")+ylab("Density")+
  guides(fill=guide_legend(keywidth=1,keyheight=1.5,default.unit="line"))+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))
dev.off()

########### Performance summary ###########
AR1Summary<-simGridBig %>% dplyr::select(Technique,sens,spec,auc,f1,ll)
sumFun<-function(x){
  xMean<-round(mean(x,na.rm=TRUE),digits=4)
  xSD<-round(sd(x,na.rm=TRUE),digits=3)
  out<-paste0(xMean," + ",xSD)
}
AR1Summary<-AR1Summary %>% group_by(Technique) %>% 
  summarize(sens=sumFun(sens),spec=sumFun(spec),sumFun(auc),sumFun(f1))

########### Sim Grid plots ##########
ggplot(simGrid %>% filter(adaptive==FALSE),
       aes(x=gammaPriors,color=as.factor(gammaPriorr),y=medLambda))+
  geom_point()+geom_line()

ggplot(simGrid %>% filter(adaptive==TRUE,adaptiveType=="norm"),
       aes(x=gammaPriors,color=as.factor(gammaPriorr),y=medLambda))+
  geom_point()+geom_line()

ggplot(simGrid %>% filter(adaptive==TRUE,adaptiveType=="priorHyper"),
       aes(x=gammaPriors,color=as.factor(gammaPriorr),y=medLambda))+
  geom_point()+geom_line()

########### Output graphs ###########
graphFun<-function(adaptive,adaptiveType,gammaPriorr,gammaPriors,
                   lambdaPriora=NULL,lambdaPriorb=NULL){
  # Sampler:
  bgl<-blockGLasso(x1,iterations=10000,burnIn=1000,adaptive=adaptive,
                   adaptiveType=adaptiveType,priorHyper=abs(solve(sim)),
                   gammaPriors=gammaPriorr,gammaPriort=gammaPriors,
                   lambdaPriora=lambdaPriora,lambdaPriorb=lambdaPriorb)
  pIBgl<-posteriorInference(bgl)
  bglMed<-pIBgl$posteriorMedian
  bglCor<-pCorFun(bglMed)
  
  # Graph
  bglG<-graph_from_adjacency_matrix(abs(bglMed),mode="undirected",diag=FALSE,weighted=TRUE)
  E(bglG)$width<-1.5*(E(bglG)$weight**2)
  bglMed[lower.tri(bglMed,diag=TRUE)]<-NA
  E(bglG)$color<-c("darkred","navyblue")[as.integer(na.omit(c(t(bglMed)))>0)+1L]
  return(bglG)
}
aIG<-graphFun(adaptive=TRUE,adaptiveType="priorHyper",gammaPriorr=1,gammaPriors=.01)
aNG<-graphFun(adaptive=TRUE,adaptiveType="norm",gammaPriorr=1,gammaPriors=.3)
rG<-graphFun(adaptive=FALSE,adaptiveType=NULL,gammaPriorr=NULL,gammaPriors=NULL,
             lambdaPriora=.01,lambdaPriorb=3)

png(file="Plots/aIG.png",height=6,width=6,units="in",res=900) # e
par(bg=NA,mar=c(0,0,0,0))
set.seed(3)
plot(aIG)
dev.off()

png(file="Plots/aNG.png",height=6,width=6,units="in",res=900) # d
par(bg=NA,mar=c(0,0,0,0))
set.seed(3)
plot(aNG)
dev.off()

png(file="Plots/rG.png",height=6,width=6,units="in",res=900) # c
par(bg=NA,mar=c(0,0,0,0))
set.seed(3)
plot(rG)
dev.off()
