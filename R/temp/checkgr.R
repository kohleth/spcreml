recltest=function(theta,neighMat,coordsMat,data,Xmat,grplist,...){
  swnames=colnames(neighMat)
  loglik=0
  c2=matrix(NA,nrow=length(grplist),ncol=2)
  i=1
  for(j in grplist){
    pointsj=which(data$sw_units3%in%j)
    
    Slist=geoR::varcov.spatial(dists.lowertri = dist(coordsMat[pointsj,]),
                               cov.model = "exp",cov.pars=theta[1:2],
                               nugget=theta[3],det=TRUE)
    Yj=data$ph[pointsj]
    Xj=Xmat[pointsj,]
    
    logdetS=Slist$log.det.to.half*2
    logdetXSinvX=determinant(crossprod(Xj,solve(Slist$varcov,Xj)),logarithm = TRUE)
    SinvX=solve(Slist$varcov,Xj)
    XSinvX=crossprod(Xj,solve(Slist$varcov,Xj))
    betahat=try(solve(XSinvX,crossprod(SinvX,Yj)),silent=TRUE)
    if(inherits(betahat,"try-error")){
      browser()
      #         ei=eigen(XSinvX,symmetric = TRUE)
      #         if(any(ei$value<0))browser()
      #         betahat=try(crossprod(t(ei$vec)/sqrt(ei$val))%*%crossprod(SinvX,Yj),silent=TRUE)
      betahat=try(qr.solve(qr(XSinvX,LAPACK = TRUE),crossprod(SinvX,Yj)),silent=TRUE)
      if(inherits(betahat,"try-error"))browser()
    }
    res=Yj-Xj%*%betahat
    ssres=crossprod(res,solve(Slist$varcov,res))
    PYj=solve(Slist$varcov,Yj)-SinvX%*%solve(XSinvX,crossprod(SinvX,Yj))
    c2[i,]=c(c(crossprod(Yj,PYj)),ssres)
    i=i+1
    
    loglik=loglik-(logdetS+logdetXSinvX$modulus+ssres)/2
  }
  list(loglik=loglik,c2=c2)
}

# debug(recltest)
temp=recltest(theta=c(0.3,700,0.3),neighMat=neighMat,coordsMat=coordsMat,data=phcDat,Xmat=Xmat,grplist=grplist)





reclone=function(theta,neighMat,coordsMat,data,Xmat,j,...){
  swnames=colnames(neighMat)
  pointsj=which(data$sw_units3%in%j)
  coordsMatj=coordsMat[pointsj,]
  Sigmaj=geoR::varcov.spatial(dists.lowertri = dist(coordsMatj),
                             cov.model = "exp",cov.pars=theta[1:2],
                             nugget=0)$varcov
  Sigmaj=Sigmaj+diag(theta[-(1:2)][as.factor(coordsMatj[,3])])
  Yj=data$ph[pointsj]
  Xj=Xmat[pointsj,]
  
  logdetS=determinant(Sigmaj,log=TRUE)$modulus
  logdetXSinvX=determinant(crossprod(Xj,solve(Sigmaj,Xj)),logarithm = TRUE)
  SinvX=solve(Sigmaj,Xj)
  XSinvX=crossprod(Xj,solve(Sigmaj,Xj))
  betahat=try(solve(XSinvX,crossprod(SinvX,Yj)),silent=TRUE)
  if(inherits(betahat,"try-error")){
    browser()
    #         ei=eigen(XSinvX,symmetric = TRUE)
    #         if(any(ei$value<0))browser()
    #         betahat=try(crossprod(t(ei$vec)/sqrt(ei$val))%*%crossprod(SinvX,Yj),silent=TRUE)
    betahat=try(qr.solve(qr(XSinvX,LAPACK = TRUE),crossprod(SinvX,Yj)),silent=TRUE)
    if(inherits(betahat,"try-error"))browser()
  }
  res=Yj-Xj%*%betahat
  ssres=crossprod(res,solve(Sigmaj,res))
  #   PYj=solve(Sigmaj,Yj)-SinvX%*%solve(XSinvX,crossprod(SinvX,Yj))
  #   c2[i,]=c(c(crossprod(Yj,PYj)),ssres)
  #   i=i+1
  
  -(logdetS+logdetXSinvX$modulus+ssres)/2
  
}

grone=function(theta,neighMat,coordsMat,data,Xmat,dS,j,...){
  swnames=colnames(neighMat)
  R=length(dS)  
  u=rep(0,R)
  H=matrix(0,nrow=R,ncol=R)
  pointsj=which(data$sw_units3%in%j)
  coordsMatj=coordsMat[pointsj,]
  Dj=as.matrix(dist(coordsMatj))
  Sigmaj=geoR::varcov.spatial(dists.lowertri = dist(coordsMatj),
                              cov.model = "exp",cov.pars=theta[1:2],nugget=0)$varcov
  Sigmaj=Sigmaj+diag(theta[-(1:2)][as.factor(coordsMatj[,3])])
  
  Yj=data$ph[pointsj]
  Xj=Xmat[pointsj,]
  SinvX=solve(Sigmaj,Xj)
  XSinvX=crossprod(Xj,solve(Sigmaj,Xj))
  qj=solve(Sigmaj,Yj)-SinvX%*%solve(XSinvX,crossprod(SinvX,Yj))
  for(r in 1:R){
    dSdthetar=dS[[r]](Dj,theta,coordsM=coordsMatj)
    Wjr=solve(Sigmaj,dSdthetar)-SinvX%*%solve(XSinvX,crossprod(SinvX,dSdthetar))
    u[r]=-0.5*sum(diag(Wjr))+0.5*crossprod(qj,dSdthetar)%*%qj
    #       for(s in r:R){
    #         Wjs=Sigmaj_list$inverse%*%dS[[s]](Dj,theta)
    #         H[r,s]=H[s,r]=H[r,s]+0.5*sum(diag(Wjr%*%Wjs))
    #       }
  }
  u
}

hsone=function(theta,neighMat,coordsMat,data,Xmat,dS,j,...){
  swnames=colnames(neighMat)
  R=length(dS)  
  u=rep(0,R)
  H=matrix(0,nrow=R,ncol=R)
  pointsj=which(data$sw_units3%in%j)
  coordsMatj=coordsMat[pointsj,]
  Dj=as.matrix(dist(coordsMatj))
  Sigmaj=geoR::varcov.spatial(dists.lowertri = dist(coordsMatj),
                              cov.model = "exp",cov.pars=theta[1:2],nugget=theta[3])$varcov
  Sigmaj=Sigmaj+diag(theta[-(1:2)][as.factor(coordsMatj[,3])])
  Yj=data$ph[pointsj]
  Xj=Xmat[pointsj,]
  SinvX=solve(Sigmaj,Xj)
  XSinvX=crossprod(Xj,solve(Sigmaj,Xj))
  #qj=solve(Sigmaj,Yj)-SinvX%*%solve(XSinvX,crossprod(SinvX,Yj))
  for(r in 1:R){
    dSdthetar=dS[[r]](Dj,theta,coordsM=coordsMatj)
    Wjr=solve(Sigmaj,dSdthetar)-SinvX%*%solve(XSinvX,crossprod(SinvX,dSdthetar))
    #u[r]=-0.5*sum(diag(Wjr))+0.5*crossprod(qj,dSdthetar,coordsM=coordsMatj)%*%qj
          for(s in r:R){
            Wjs=solve(Sigmaj,dS[[s]](Dj,theta,coordsM=coordsMatj))-SinvX%*%solve(XSinvX,crossprod(SinvX,dS[[s]](Dj,theta,coordsM=coordsMatj)))
            H[r,s]=H[s,r]=H[r,s]-0.5*traceWW(Wjr,Wjs)
          }
  }
  H
}

comparegr2=foreach(jj=grplist,.packages=c("geoR","Matrix"))%dopar%{
  tempgrn1=grone(theta=c(0.3,700,rep(0.3,6)),neighMat=neighMat,coordsMat=coordsMat,data=phcDat,Xmat=Xmat,dS=dSexph,j=jj)
  numgn1=numDeriv::grad(reclone,c(0.3,700,rep(0.3,6)),neighMat=neighMat,coordsMat=coordsMat,data=phcDat,Xmat=Xmat,dS=dSexph,j=jj)
  list(num=numgn1,analytic=tempgrn1)
}



comparehs=foreach(jj=grplist,.packages=c("geoR","Matrix"))%dopar%{
  temphsn1=hsone(theta=c(0.3,700,rep(0.3,6)),neighMat=neighMat,coordsMat=coordsMat,data=phcDat,Xmat=Xmat,dS=dSexph,j=jj)
  numhsn1=numDeriv::hessian(reclone,c(0.3,700,rep(0.3,6)),neighMat=neighMat,coordsMat=coordsMat,data=phcDat,Xmat=Xmat,dS=dSexph,j=jj)
  list(num=numhsn1,analytic=temphsn1)
}

sapply(comparegr,function(x)max(range(abs(x[[1]]-x[[2]])))<1e-5)
sapply(comparehs,function(x)max(range(abs(x[[1]]-x[[2]])))<1e-5)

sapply(1:(ncol(neighMat)-1),function(a){
  k=colnames(neighMat)[a]
  N_k=which(neighMat[k,]==TRUE)
  
  N_k=N_k[as.numeric(names(N_k))>as.numeric(k)]
  length(which(phcDat$sw_units3%in%c(N_k,k)))
})

debug(grone)
tempgrn1=grone(theta=c(0.3,700,0.3),neighMat=neighMat,coordsMat=coordsMat,data=phcDat,swnames=colnames(neighMat),Xmat=model.matrix(trend),dS=dSexp,R=3,jj=7)
