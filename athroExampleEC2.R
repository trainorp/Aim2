########### Prereqs ###########
options(stringsAsFactors=FALSE)
library(methods)
library(MASS)
library(BayesianGLasso)
library(tidyverse)
library(igraph)

setwd("~/gdrive/Dissertation/Aim2")

# Import data:
load(file="atheroExampleV3Data.RData")

# Follow-up and with annotation only:
df1<-df1[df1$timepoint=="TF-U" | (df1$group=="sCAD" & df1$timepoint=="T0"),]
m1<-as.matrix(df1[,!names(df1) %in% c("group","timepoint","ptid")])
m1<-scale(m1,center=TRUE,scale=FALSE)

# Entropy filter:
m1<-m1[,apply(m1,2,function(x) length(unique(x))>14)]

# Filter for those without structural information:
m1<-m1[,colnames(m1) %in% colnames(simMat)]
print(dim(m1))

# Make sure column order / names match:
simMat<-simMat[rownames(simMat) %in% colnames(m1),colnames(simMat) %in% colnames(m1)]
print(dim(simMat))

# Structure Adaptive:
priorHyper<-2*(simMat**2)+.1

ptm<-proc.time()
aiBGL1<-blockGLasso(m1,iterations=1000,burnIn=250,adaptive=TRUE,
                    adaptiveType="priorHyper",priorHyper=priorHyper,
                    lambdaii=62,gammaPriors=8,gammaPriort=.001)
proc.time()-ptm

save(aiBGL1,file="aiBGL1.RData")
