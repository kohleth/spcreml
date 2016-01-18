library(geoR)
library(foreach)
library(doParallel)

N=500

sim=grf(N,cov.pars=c(0.5,3000),mean=5,nugget=0.3,xlim=c(0,100000),ylim=c(0,100000))
sim$borders=NULL
# reml <- likfit(sim, ini=c(0.5, 4000), fix.nug = FALSE, lik.met = "REML")
# summary(reml)

# ml <- likfit(sim, ini=c(0.5, 4000), fix.nug = FALSE, lik.met = "ML")

## check recl---------------
trend=lm(data~1,data=sim,y=TRUE,x=TRUE)
Xmat=trend$x
Y=trend$y
Ngrp=4
grp=kmeans(sim[[1]],Ngrp,nstart=100)
plot(sim[[1]],col=grp$cluster)
grplist=combn(1:Ngrp,2,simplify=FALSE)
pointsk=lapply(1:Ngrp,function(x)which(grp$cluster%in%x))
pointsj=lapply(grplist,function(x)which(grp$cluster%in%x))
g2list=c(combn(1:length(grplist),2,simplify = FALSE),combn(1:length(grplist),3,simplify = FALSE))

dSexpds2=function(D,theta0,...)exp(-D/theta0[2])
dSexpdphi=function(D,theta0,...)D*theta0[1]/theta0[2]^2*exp(-D/theta0[2])

dSsphds2=function(D,theta0,...){
  d=D/theta0[2]
  out=1-1.5*d+0.5*d^3
  out[D>theta0[2]]=0
  out
}
dSsphdphi=function(D,theta0,...){
  out=theta0[1]*(1.5*D/theta0[2]^2-1.5*D^3/theta0[2]^4)
  out[D>theta0[2]]=0
  out
}

dSdtausq=function(D,theta0,...)diag(1,nrow=nrow(D))
dSdtausq1=function(D,theta0,coordsM,...)diag(coordsM[,3]==0.025,nrow=nrow(D))
dSdtausq2=function(D,theta0,coordsM,...)diag(coordsM[,3]==0.1,nrow=nrow(D))
dSdtausq3=function(D,theta0,coordsM,...)diag(coordsM[,3]==0.225,nrow=nrow(D))
dSdtausq4=function(D,theta0,coordsM,...)diag(coordsM[,3]==0.45,nrow=nrow(D))
dSdtausq5=function(D,theta0,coordsM,...)diag(coordsM[,3]==0.8,nrow=nrow(D))
dSdtausq6=function(D,theta0,coordsM,...)diag(coordsM[,3]==1.5,nrow=nrow(D))

dSexp=list(dSexpds2,dSexpdphi,dSdtausq)
dSexph=list(dSexpds2,dSexpdphi,dSdtausq1,dSdtausq2,dSdtausq3,dSdtausq4,dSdtausq5,dSdtausq6)
dSsph=list(dSsphds2,dSsphdphi,dSdtausq)
dSsphh=list(dSsphds2,dSsphdphi,dSdtausq1,dSdtausq2,dSdtausq3,dSdtausq4,dSdtausq5,dSdtausq6)


cl=makeCluster(3)
registerDoParallel(cl)

fngr=gr(theta = c("sigma2"=0.5,"phi"=3000,"nugget"=0.3),Y = Y,coordsMat = sim[[1]],Xmat = Xmat,dS = dSexp,pointsj = pointsj,cov.model = "exp")
numgr=foreach(jj=1:length(pointsj))%do%{
  numDeriv::grad(recl,c("sigma2"=0.5,"phi"=3000,"nugget"=0.3),Y = Y,coordsMat = sim[[1]],Xmat = Xmat,dS = dSexp,pointsj = pointsj[jj],cov.model = "exp")
}

stopCluster(cl)
