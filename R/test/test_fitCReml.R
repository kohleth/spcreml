library(geoR)
library(foreach)
library(doParallel)

N=500

sim=grf(N,cov.pars=c(0.5,3000),mean=5,nugget=0.3,xlim=c(0,100000),ylim=c(0,100000))
sim$borders=NULL
reml <- likfit(sim, ini=c(0.5, 4000), fix.nug = FALSE, lik.met = "REML")
summary(reml)

ml <- likfit(sim, ini=c(0.5, 4000), fix.nug = FALSE, lik.met = "ML")

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

cl=makeCluster(3)
registerDoParallel(cl)

simfit1=fitCReml(form = data~1,data=list(data=sim$data),dS = dSexp,coordsMat = sim[[1]],
                init = c("sigma2"=0.5,"phi"=4000,"nugget"=0.1),pointsj=pointsj,pointsjpair=g2list,pointsk=pointsk,cov.model="exp",lower=c(1e-5,1e-5,1e-5))

stopCluster(cl)
