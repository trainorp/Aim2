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

# png(filename="Plots/lambdasVsSim.png",width=4.5,height=3.5,units="in",res=600)
ggplot(data=df1,aes(x=omega,y=lambda,color=dis,group=dis))+geom_line()+theme_bw()+
  xlab(expression(paste("|",tilde(omega)[ij],"|")))+
  ylab(expression(paste(E(lambda["ij"]))))+
  scale_color_discrete(name="Similarity")
# dev.off()

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
  
  # Posterior inference object:
  pIBgl<-posteriorInference(bgl)
  
  # Posterior median:
  medBgl<-pIBgl$posteriorMedian
  medBglSigma<-solve(medBgl)
  
  # Topological error analysis:
  pCorsMedBgl<-pCorFun(medBgl)
  pCorsIndMedBgl<-abs(pCorsMedBgl)>.1
  tabBgl<-xtabs(~true+pred,data=data.frame(true=c(pCorsInd),pred=c(pCorsIndMedBgl)))
  simGrid$sens[i]<-tabBgl['TRUE','TRUE']/sum(tabBgl['TRUE',])
  simGrid$spec[i]<-tabBgl['FALSE','FALSE']/sum(tabBgl['FALSE',])
  simGrid$ppv[i]<-tabBgl['TRUE','TRUE']/sum(tabBgl[,'TRUE'])
  simGrid$f1[i]<-(2*simGrid$sens[i]*simGrid$ppv[i])/(simGrid$sens[i]+simGrid$ppv[i])
}
