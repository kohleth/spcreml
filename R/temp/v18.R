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
library(openxlsx)
source("fns5.R")



swunitdata=read.xlsx("Data/sw_pHc_covdata.xlsx",detectDates = TRUE)
data0=read.xlsx("Data/sw_pHc_covariates.xlsx",detectDates=TRUE)

data1=data.frame(data0,swunitdata[,-(1:10)])
data1=convColClass(data1,setdiff(8:79,15:17),as.numeric)
data1=convColClass(data1,15:17,as.factor)
data1=convColClass(data1,76:78,as.factor)
### check if I have .tiff for all covariates
# covnames=names(data0)[11:73]
# tifnames=gsub(".tif","",list.files("P://cmi105264/Working/Nathan/sw_covariates",pattern=glob2rx("*.tif")))
# covmatchup=data.frame(covnames,"match"=tifnames[match(covnames,tifnames)],"pmatch"=tifnames[pmatch(covnames,tifnames)])
# for(nn in 1:length(covnames)){
#   if(is.na(covmatchup[nn,"pmatch"]))data1[,covnames[nn]]=NULL
# }

S1sp=readOGR("C:/Users/kc58/pHmap/Covariates/swunits/swunits3.shp",layer="swunits3")
# S2sp=readOGR("C:/Users/kc58/pHmap/Covariates/swunits/swunits4.shp",layer="swunits4")
# S3sp=readOGR("C:/Users/kc58/pHmap/Covariates/swunits/swunits5.shp",layer="swunits5")

# elevationrt=raster("../Covariates/elevation.tif")
# elevationpt=rasterToPoints(elevationrt,spatial=TRUE)
# elevationsp=readGDAL("../Covariates/elevation.tif")

# nonNAid=which(!is.na(values(elevationrt)))
# if(length(nonNAid)!=nrow(elevationpt))stop("check this line!")

## make easier names
phcDat=data1
names(phcDat)=tolower(names(phcDat))
phcDat=rename(phcDat,c("x_vg94"="x","y_vg94"="y","phc_value"="ph"))
phcDat=phcDat[complete.cases(phcDat),]
phcDat[,c("x","y")]=geoR::jitter2d(phcDat[,c("x","y")],max=1)

## add depth
phcDat$depth=(with(phcDat,(lowerdepth+upperdepth)/2/100))
phcDat=phcDat[order(phcDat$depth,phcDat$x,phcDat$y),]

## add scaled coordinates
xylim=bbox(S1sp)
phcDat$sx=(phcDat$x-xylim["x","min"])/diff(xylim["x",])
phcDat$sy=(phcDat$y-xylim["y","min"])/diff(xylim["y",])
# phcDat$sx=(phcDat$x)/1e6
# phcDat$sy=(phcDat$y)/1e6


## find Sp neighbours
S1sp=subset(S1sp,sapply(S1sp@polygons,slot,"area")>100^2)
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


### detrending -------------------
trend=lm(ph~(sx+sy)*depth+(veg_fpar_mean_pvv2+slope+pot+etaaan+prescott)*depth,data=phcDat,x=TRUE)

coordsDf=phcDat[,c("sx","sy","depth")]
# coordsDf=phcDat[,c("x","y","depth")]

coordsMat=as.matrix(coordsDf)


rg=1e4
Zsp=data.frame(z=resid(trend),phcDat[,c("sx","sy","depth")])
coordinates(Zsp)=~sx+sy+depth
Zvg=variogram(z~1,data=Zsp,cressie=FALSE,cutoff=0.1)
plot(Zvg$dist,Zvg$gamma,ylim=c(0,1),xlab="standardized distance (m)",ylab="semi-variance")
lines(variogramLine(vgm(psill = 0.28222,model = "Exp",range = 0.00438,nugget = 0.26507),maxdist=0.1))

gstatvg=fit.variogram(Zvg,model=vgm(psill = 0.3,model = "Exp",range = 0.005,nugget = 0.25))
lines(gstatvg,col=2)
lines(variogramLine(vgm(psill = 0.57,model = "Exp",range = 0.00394,nugget = 0),maxdist=0.1),col=2)
# 
# for(i in 1:6){
#   lines(variogramLine(vgm(psill = 0.243295583,model = "Exp",range = 0.004988,nugget = fit2$par[i+2]),maxdist=0.1),col=i+7)
#   
# }




for(d in 1:6){
  Zsp=data.frame(z=resid(trend),phcDat[,c("sx","sy","depth")])
  Zsp=subset(Zsp,depth==unique(phcDat$depth)[d])
  coordinates(Zsp)=~sx+sy+depth
  Zvg=variogram(z~1,data=Zsp,cressie=FALSE,cutoff=0.1)
  points(Zvg$dist,Zvg$gamma,ylim=c(0,0.9),pch=3,col=d+7)
}



cl=makeCluster(20)
registerDoParallel(cl)

fit1=fitCLSM(form = formula(trend),data = phcDat,
             dS = dSexp,coordsMat = coordsMat,init = c(0.5,0.01,0.3),
             pointsj = pointsj,pointsjpair = g2list,
             pointsk = NULL,cov.model = "exp",
             heter = FALSE,lower=rep(1e-3,3))
save(fit1,file="./modelfit/fit1.rda")


fit2=fitCLSM(form = formula(trend),data = phcDat,
             dS = dSexph,coordsMat = coordsMat,init = c(0.5,0.01,rep(0.3,6)),
             pointsj = pointsj,pointsjpair = g2list,
             pointsk = NULL,cov.model = "exp",
             heter = FALSE,lower=rep(1e-3,8))
save(fit2,file="./modelfit/fit2.rda")

stopCluster(cl)



### prediction
# using gstat for prediction
fixedrt=stack(paste0("C:/Cubist1/Covariates","/",c("slope","pot","etaaan","prescott","veg_fpar_mean_pvv2"),".tif"))
fixedDF0=as.data.frame(rasterToPoints(fixedrt))
fixedDF0$sy=(fixedDF0$y-xylim["y","min"])/diff(xylim["y",])
fixedDF0$sx=(fixedDF0$x-xylim["x","min"])/diff(xylim["x",])



# fixedDF=do.call(rbind,rep(list(fixedDF0),length(phcDat$depth)))
# fixedDF$depth=unique(phcD)

phcDatsp=phcDat
coordinates(phcDatsp)=~sx+sy+depth

cl=makeCluster(20)
registerDoParallel(cl)
depths=unique(phcDat$depth)
krg=foreach(i=1:length(depths))%:%
  foreach(df=isplitRows(fixedDF0,chunks = length(cl)),.packages=c("gstat","sp"),.combine="rbind")%dopar%{
    df$depth=depths[i]
    dfsp=df
    coordinates(dfsp)=~sx+sy+depth
    krige(formula(trend),locations=phcDatsp,newdata=dfsp,
          model=vgm(fit1$par[1],model = "Exp",range = fit1$par[2],nugget = fit1$par[3]),
          beta=fit1$betahat,maxdist=fit1$par[2]*2)
  }


stopCluster(cl)


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
writeRaster(rt,filename = "prediction.tif")
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
writeRaster(rtse,filename = "predictionse.tif")
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

# newdf=phcDat[1:1000,]
# edges=which(neighMat,arr.ind = TRUE,useNames = T)
# newpointsk=lapply(colnames(neighMat),function(x)which(phcDat$sw_units3[1:1000]==x))
# Nk=apply(neighMat,1,which)
# 
# newpointskl=apply(edges,1,function(x)which(phcDat$sw_units3[1:1000]%in%colnames(neighMat)[x]))
# emptypair=which(sapply(newpointskl,length)==0)
# newpointskl=newpointskl[-emptypair]
# 
# pred=predict.reclsm(fit = fit2,newdf = newdf,newpointsk = newpointsk,
#                     pointsk=pointsk,Nk = Nk,newcoordsMat = coordsMat[1:1000,])
