library(geoR)

N=500

sim=grf(N,cov.pars=c(0.5,3000),mean=5,nugget=0.3,xlim=c(0,100000),ylim=c(0,100000))
sim$borders=NULL
reml <- likfit(sim, ini=c(0.5, 4000), fix.nug = FALSE, lik.met = "REML")
summary(reml)

ml <- likfit(sim, ini=c(0.5, 4000), fix.nug = FALSE, lik.met = "ML")

plot(variog(sim))
lines(reml)

Ngrp=4
grp=kmeans(sim[[1]],Ngrp,nstart=100)
plot(sim[[1]],col=grp$cluster)
grplist=combn(1:Ngrp,2,simplify=FALSE)
pointsk=lapply(1:Ngrp,function(x)which(grp$cluster%in%x))
pointsj=lapply(grplist,function(x)which(grp$cluster%in%x))
g2list=list()
for(ii in 1:(Ngrp-1)){
  for(jj in (ii+1):Ngrp){
    if(length(intersect(grplist[[ii]],grplist[[jj]]))>0){
      g2list=c(g2list,list(c(ii,jj)))
    }
  }
}
g2list=c(combn(1:length(grplist),2,simplify = FALSE),combn(1:length(grplist),3,simplify = FALSE))
simfit1=fitCLSM(form = data~1,data=list(data=sim$data),dS = dSexp,coordsMat = sim[[1]],
        init = c(0.5,4000,0.1),pointsj=pointsj,pointsjpair=g2list,pointsk=pointsk,cov.model="exp",heter=FALSE,lower=c(1e-5,1e-5,1e-5))

format(data.frame("ml"=c(ml$cov.pars,ml$nugget,ml$beta,ml$beta.var),
           "reml"=c(reml$cov.pars,reml$nugget,reml$beta,reml$beta.var),
           "recl"=c(simfit1$par,simfit1$beta,simfit1$varbetahat),
           "true"=c(0.5,3000,0.3,5,NA),row.names=c("s2","phi","tau2","beta","var(beta)")),scientific=FALSE)


#### kriging
bbox=apply(sim[[1]],2,range)
predgr=expand.grid(seq(bbox[1,1],bbox[2,1],l=50),seq(bbox[1,2],bbox[2,2],l=50))
Nk=list(c(2,3,4),c(1,3,4),c(1,2,4),c(1,2,3))
newcluster=kmeans(predgr,4,nstart = 100)
newpointsk=lapply(1:Ngrp,function(x)which(newcluster$cluster==x))

mlkgcl=krige.control(obj.model=ml)
mlkg=krige.conv(geodata=sim,locations=predgr,krige = mlkgcl)
remlkgcl=krige.control(obj.model=reml)
remlkg=krige.conv(geodata=sim,locations=predgr,krige = remlkgcl)

colnames(predgr)=colnames(simfit1$coordsMat)
newdf=data.frame(predgr,data=1)

simpred=predict.reclsm(fit = simfit1,newdf = newdf,newpointsk = newpointsk,pointsk = pointsk,
               Nk =Nk,newcoordsMat = predgr)

plot(simpred,remlkg$predict)

plot(predgr,col=newcluster$cluster)
points(predgr[tempid,],col=6,pch=2)
points(sim[[1]],pch=3)

## try to use gstat
library(gstat)
data(meuse)
locs=meuse[rep(1,nrow(sim[[1]])),]
locs[,c(1:2)]=sim[[1]]
locs$z=sim[[2]]
coordinates(locs)=~x+y
data(meuse.grid)
newdat=meuse.grid[rep(1,nrow(predgr)),]
newdat[,c(1:2)]=predgr
gridded(newdat) = ~x+y
gstatkg=gstat::krige(z~1,locs,newdat,beta=simfit1$betahat[1],
      model=vgm(simfit1$par[1],model="Exp",range=simfit1$par[2],nugget=simfit1$par[3]))

### now try large dataset
data(meuse)
locs=meuse[rep(1,nrow(sim[[1]])),]
locs[,c(1:2)]=sim[[1]]
locs$z=sim[[2]]
coordinates(locs)=~x+y
data(meuse.grid)
Bpredgr=expand.grid(seq(bbox[1,1],bbox[2,1],l=115),seq(bbox[1,2],bbox[2,2],l=115))
Bnewdat=meuse.grid[rep(1,nrow(Bpredgr)),]
Bnewdat[,c(1:2)]=Bpredgr
gridded(Bnewdat) = ~x+y
Bgstatkg=gstat::krige(z~1,locs,Bnewdat,
                     model=vgm(simfit1$par[1],model="Exp",range=simfit1$par[2],nugget=simfit1$par[3]))

reclkgcl=krige.control(obj.model=reml)
reclkgcl$beta=simfit1$betahat
reclkgcl$cov.pars=simfit1$par[1:2]
reclkgcl$nugget=simfit1$par[3]
geoRreclpred=krige.conv(geodata=sim,locations=predgr,krige = remlkgcl)

plot(remlkg$pred,geoRreclpred$pred)
