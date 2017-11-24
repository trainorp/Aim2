library(clusterGeneration)

p<-10
n<-10000
mu<-1

set.seed(3)
cMat<-genPositiveDefMat(p,covMethod = "eigen")$Sigma
omegaMat<-solve(cMat)
s21<-solve(-cMat,rep(mu,p))

set.seed(33)
x1<-mvrnorm(n=n,mu=t(-cMat %*% s21),Sigma=cMat)

set.seed(333)
rnorm1<-stats::rnorm(p)
x2<-t(-cMat %*% s21)+solve(cMat,rnorm1)
for(i in 2:n)
{
  set.seed(333+i)
  x2<-rbind(x2,t(-cMat %*% s21)+solve(cMat,stats::rnorm(p)))
}

covX1<-cov(x1)
covX2<-cov(x2)
