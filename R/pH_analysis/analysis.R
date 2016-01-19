library(latticeExtra)
library(grid)
library(raster)
library(rgdal)
library(plyr)
library(Matrix)
library(doParallel)
library(foreach)
library(itertools)
library(gstat)
library(rgeos)


load("pHdata.rda")
load("blocking.rda")

## find Sp neighbours
neighMat=gTouches(S1sp,byid=T)
colnames(neighMat)=rownames(neighMat)=S1sp@data[,1]

### 
pointsk=lapply(colnames(neighMat),function(x)which(phcDat$sw_units3==x))

grplist=formGrps(neighMat)
pointsj=lapply(grplist,function(x)which(phcDat$sw_units3%in%x))


## forming second level neighbourhood for var(beta)
neigh2=matrix(NA,nrow=length(pointsj),ncol=length(pointsj))
for(i in 1:length(pointsj)){
  for(ii in i:length(pointsj)){
    neigh2[i,ii]=any(duplicated(c(pointsj[[i]],pointsj[[ii]])))
  }
}
g2list=alply(which(neigh2,arr.ind = T),1)


### find trend form  -------------------
trend=lm(ph~(sx+sy)*depth+(veg_fpar_mean_pvv2+slope+pot+etaaan+prescott)*depth,data=phcDat,x=TRUE)

coordsDf=phcDat[,c("sx","sy","depth")]
# coordsDf=phcDat[,c("x","y","depth")]

coordsMat=as.matrix(coordsDf)


### C-REML---------------------
cl=makeCluster(min(c(20,detectCores())))
registerDoParallel(cl)

fit1=fitCLSM(form = formula(trend),data = phcDat,
             dS = dSexp,coordsMat = coordsMat,init = c(0.5,0.01,0.3),
             pointsj = pointsj,pointsjpair = g2list,
             pointsk = NULL,cov.model = "exp",
             heter = FALSE,lower=rep(1e-3,3))

## Loading covariates for fixed effect prediction--------
fixedrt=stack("covar.tif")
fixedDF0=as.data.frame(rasterToPoints(fixedrt))
fixedDF0$sy=(fixedDF0$y-xylim["y","min"])/diff(xylim["y",])
fixedDF0$sx=(fixedDF0$x-xylim["x","min"])/diff(xylim["x",])

## preparing sp object for prediction----
phcDatsp=phcDat
coordinates(phcDatsp)=~sx+sy+depth


## local kriging -----------------------
cl=makeCluster(20)
registerDoParallel(cl)
depths=unique(phcDat$depth)
## SK
krg=foreach(i=1:length(depths))%:%
  foreach(df=isplitRows(fixedDF0,chunks = length(cl)),.packages=c("gstat","sp"),.combine="rbind")%dopar%{
    df$depth=depths[i]
    dfsp=df
    coordinates(dfsp)=~sx+sy+depth
    krige(formula(trend),locations=phcDatsp,newdata=dfsp,
          model=vgm(fit1$par[1],model = "Exp",range = fit1$par[2],nugget = fit1$par[3]),
          beta=fit1$betahat,maxdist=fit1$par[2]*2)
  }


## convert result from sp into raster
nonNAid=which(!is.na(values(fixedrt[[1]])))

rt=foreach(i=1:length(depths),.combine=stack)%do%{
  v=rep(NA,ncell(fixedrt))
  v0=krg[[i]][[1]]
  ext=quantile(v0,probs=c(0.01,0.99),na.rm=TRUE)
  v0[v0<ext[1]]=ext[1]
  v0[v0>ext[2]]=ext[2]
  v[nonNAid]=v0
  setValues(raster(fixedrt),values = v)
}

names(rt)=paste(depths,"m")
splot=spplot(rt,layout=c(1,length(depths)),
             col.regions=heat.colors,
             at=do.breaks(range(phcDat$ph),16))
splot$condlevels$name=paste(depths,"m")

## prediction error
rtse=foreach(i=1:length(depths),.combine=stack)%do%{
  v=rep(NA,ncell(fixedrt))
  v0=sqrt(krg[[i]][[2]])
  v[nonNAid]=v0
  setValues(raster(fixedrt),values = v)
}

names(rtse)=paste(depths,"m")
sesplot=spplot(rtse,layout=c(1,length(depths)))
sesplot$condlevels$name=paste(depths,"m")

dataplot=splot2(phcDat$ph,phcDat)
dataplot=update(dataplot,scale=list(x=list(draw=FALSE),y=list(draw=FALSE)),ylab="",xlab="")
dataplot$condlevels[[1]]=paste(depths,"m")


print(dataplot,split=c(1,1,3,1),more=TRUE)
print(splot,split=c(2,1,3,1),more=TRUE)
print(sesplot,split=c(3,1,3,1))
grid.text("data", x = 0.13, y = 0.99, just = c("centre", "top"),gp=gpar(fontsize=12))
grid.text("prediction", x = 0.48, y = 0.99, just = c("centre", "top"),gp=gpar(fontsize=12))
grid.text("prediction s.e.", x = 0.8, y = 0.99, just = c("centre", "top"),gp=gpar(fontsize=12))
