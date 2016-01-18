splot2=function(var,data){
  nvar=length(unique(data$depth))
  Vcol=level.colors(var,at=do.breaks(range(var),16),heat.colors)
  xyplot(y~x|as.factor(depth),data=data,groups=Vcol,
         aspect="iso",key.auto=TRUE,
         panel = function(x, y, groups, ..., subscripts) {
           fill <- groups[subscripts]
           panel.grid(h = -1, v = -1)
           panel.xyplot(x, y, pch = 15,col=fill,cex=0.1,  ...)
         },layout=c(1,nvar),index.cond=list(nvar:1),
         legend=list(right=list(fun=draw.colorkey,
                                args=list(key=list(col=heat.colors,
                                                   at=do.breaks(range(var),16)),
                                          draw=FALSE)))
  )
}

vecToRt=function(Y){
  for(i in 1:ndepths){
    temp=rep(NA,length(elevationrt))
    temp[nonNAid]=Y[(i-1)*newN/ndepths+(1:(newN/ndepths))]
    assign(paste0("depth",i),setValues(elevationrt,temp))
  }
  stack(mget(paste0("depth",1:ndepths)))
}


convColClass=function(dd,cols,fn){
  for(jj in cols)dd[,jj]=do.call(fn,list(dd[,jj]))
  dd
}

## functions to compute derivative of Sigma
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

traceWW=function(A,B){
  stopifnot(ncol(A)==nrow(B))
  sum(sapply(1:nrow(A),function(ii)A[ii,]%*%B[,ii]))
}

# n=100
# mat1=matrix(rnorm(n^2),n)
# mat2=matrix(rnorm(n^2),n)
# 
# microbenchmark::microbenchmark(traceWW(mat1,mat2),sum(diag(mat1%*%mat2)))
#identical(traceWW(mat1,mat2), sum(diag(mat1%*%mat2)))

formGrps=function(neighMat,minS=100){
  swnames=colnames(neighMat)
  grplist0=list()
  ## form initial group
  for(k in swnames[-length(swnames)]){
    N_k=which(neighMat[k,]==TRUE)
    N_k=N_k[as.numeric(names(N_k))>as.numeric(k)]
    if(length(N_k)>0){
      grp=lapply(names(N_k),function(l){
        out=c(l,k)
        attr(out,"N")=length(which(phcDat$sw_units3%in%out))
        out
      })
      Ngrp=sapply(grp,attr,"N")
      if(min(Ngrp)<100){
        grp=unique(unlist(grp))
        attr(grp,"N")=length(which(phcDat$sw_units3%in%grp))
      }
      grplist0=c(grplist0,list(grp))
    }
  }
  out=list()
  for(i in 1:length(grplist0)){
    if(inherits(grplist0[i][[1]],"list")){
      out=c(out,unlist(grplist0[i],FALSE))
    }else{
      out=c(out,grplist0[i])
    }
  }
  ## if small group still exist
  Nout=sapply(out,attr,"N")
  smallgrp=which(Nout<50)
  if(length(smallgrp)>0){
    N100=which(Nout>99&Nout<200)[1]
    combn=Reduce(union,out[c(smallgrp,N100)])
    attr(combn,"N")=length(which(phcDat$sw_units3%in%combn))
    out[smallgrp]=NULL
    out=c(out,list(combn))
  }
  ## check N is correct
  stopifnot(sapply(out,function(x)length(which(phcDat$sw_units3%in%x))==attr(x,"N")))
  return(out)
}

genSigmaj=function(coordsM,theta,cov.model,heter){
  if(heter){
    Sigmaj=geoR::varcov.spatial(dists.lowertri = dist(coordsM),
                                cov.model = cov.model,cov.pars=theta[1:2],
                                nugget=0)$varcov
    Sigmaj=Sigmaj+diag(theta[-(1:2)][as.factor(coordsM[,3])])
  }else{
    Sigmaj=geoR::varcov.spatial(dists.lowertri = dist(coordsM),
                                cov.model = cov.model,cov.pars=theta[1:2],
                                nugget=theta[3])$varcov
  }
  return(Sigmaj)
}

dist2=function(x,s=1e6,...){
  x[,1:2]=x[,1:2]*s
  dist(x,...)
}

recl=function(theta,coordsMat,Y,Xmat,pointsj,heter=TRUE,cov.model="exp",...){
  recloglik=foreach(pointsj=pointsj,.packages=c("geoR","Matrix"),.combine="+",.export=c("dist","genSigmaj"))%dopar%{
    coordsMatj=coordsMat[pointsj,]
    
    Sigmaj=genSigmaj(coordsM=coordsMatj,theta=theta,cov.model=cov.model,heter=heter)
    
    Yj=Y[pointsj]
    Xj=Xmat[pointsj,]
    
    logdetS=determinant(Sigmaj,log=TRUE)$modulus
    SinvX=solve(Sigmaj,Xj)
    XSinvX=crossprod(Xj,solve(Sigmaj,Xj))
    logdetXSinvX=determinant(XSinvX,logarithm = TRUE)$modulus
    PYj=solve(Sigmaj,Yj)-SinvX%*%solve(XSinvX,crossprod(SinvX,Yj))
    
    -(logdetS+logdetXSinvX+c(crossprod(Yj,PYj)))/2
  }
  return(recloglik)
}

grHs=function(theta,coordsMat,Y,Xmat,dS,pointsj,heter=TRUE,cov.model="exp",gr=TRUE,hs=TRUE,...){
  R=length(dS)
  out=foreach(pointsj=pointsj,.packages=c("geoR"),.export=c("dist","traceWW","genSigmaj"))%dopar%{
    u=rep(0,R)
    H=matrix(0,nrow=R,ncol=R)
    coordsMatj=coordsMat[pointsj,]
    Dj=as.matrix(dist(coordsMatj))
    
    Sigmaj=genSigmaj(coordsM=coordsMatj,theta=theta,cov.model=cov.model,heter=heter)
    
    Yj=Y[pointsj]
    Xj=Xmat[pointsj,]
    SinvX=solve(Sigmaj,Xj)
    XSinvX=crossprod(Xj,solve(Sigmaj,Xj))
    
    if(gr){
      qj=solve(Sigmaj,Yj)-SinvX%*%solve(XSinvX,crossprod(SinvX,Yj))
    }
    for(r in 1:R){
      dSdthetar=dS[[r]](Dj,theta,coordsM=coordsMatj)
      Wjr=solve(Sigmaj,dSdthetar)-SinvX%*%solve(XSinvX,crossprod(SinvX,dSdthetar))
      if(gr){
        u[r]=u[r]-0.5*sum(diag(Wjr))+0.5*crossprod(qj,dSdthetar)%*%qj
      }
      if(hs){
        for(s in r:R){
          Wjs=solve(Sigmaj,dS[[s]](Dj,theta,coordsM=coordsMatj))-SinvX%*%solve(XSinvX,crossprod(SinvX,dS[[s]](Dj,theta,coordsM=coordsMatj)))
          H[r,s]=H[s,r]=H[r,s]-0.5*traceWW(Wjr,Wjs)
        }        
      }
    }
    list(H=H,u=u)
  }
  if(gr&hs)return(list(gr=Reduce("+",lapply(out,function(x)x$u)),
                       hs=Reduce("+",lapply(out,function(x)x$H))))
  if(gr)return(Reduce("+",lapply(out,function(x)x$u)))
  if(hs)return(Reduce("+",lapply(out,function(x)x$H)))
}

gr=function(theta,coordsMat,Y,Xmat,dS,pointsj,heter=TRUE,cov.model,...){
  grHs(theta,coordsMat,Y,Xmat,dS,pointsj,heter=heter,cov.model=cov.model,gr=TRUE,hs=FALSE)
}


Hs=function(theta,coordsMat,Y,Xmat,dS,pointsj,heter=TRUE,cov.model,...){
  grHs(theta,coordsMat,Y,Xmat,dS,pointsj,heter=heter,cov.model,gr=FALSE,hs=TRUE)
}


## just getting the score function
calcAB=function(pointsj,Y,coordsMat,theta,Xmat,cov.model,heter,...){
  out=foreach(pointsj=pointsj,.packages="geoR",.export=c("genSigmaj"))%dopar%{
    A=B=0
    Sigmaj=genSigmaj(coordsMat[pointsj,],theta = theta,cov.model = cov.model,heter = heter)
    Yj=Y[pointsj]
    Xj=Xmat[pointsj,]
    A=A+crossprod(Xj,solve(Sigmaj,Xj))
    B=B+crossprod(Xj,solve(Sigmaj,Yj))
    list(A=A,B=B)   
  }
  A=Reduce("+",lapply(out,function(x)x$A))
  B=Reduce("+",lapply(out,function(x)x$B))
  return(list(A=A,B=B))
}

calcJB=function(g2list,pointsj,Y,coordsMat,theta,Xmat,cov.model,heter,...){
  out=foreach(g2=g2list,.packages="geoR",.export=c("genSigmaj"),.combine="+")%dopar%{
#   g2=g2list[[1]]
    j=pointsj[[g2[[1]]]]
    j2=pointsj[[g2[[2]]]]
    Xj=Xmat[j,]
    Xj2=Xmat[j2,]
    Sj=genSigmaj(coordsM = coordsMat[j,],theta = theta,cov.model = cov.model,heter = heter)
    Sj2=genSigmaj(coordsM = coordsMat[j2,],theta = theta,cov.model = cov.model,heter = heter)
    QXj=solve(Sj,Xj)
    QXj2=solve(Sj2,Xj2)
#     dist2set=function(jj,jj2)dist(coordsMat[c(jj,jj2),])
#     Djj2=outer(j,j2,FUN = Vectorize(dist2set))
  D0=as.matrix(dist(coordsMat[c(j,j2),]))
  Djj2=D0[1:length(j),length(j)+(1:length(j2))]
    covYjYj2=theta[1]*exp(-Djj2/theta[2])
    crossprod(QXj,covYjYj2)%*%QXj2
  }
}

fitCLSM=function(form,data,dS,coordsMat,init,pointsj,pointsjpair,pointsk,cov.model=cov.model,heter=FALSE,lower=NULL){
  
  ## pointsj is a list, each component is a vector of indices (an integer from 1 to N) that indicates which row of data belong to sub-region j (which is itself a joint of 2 partitions k and l).
  ## pointsjpair is a list, each component is a pair of indices of pointsj (i.e. an integer from 1 to length(pointsj)) that indicates which 2 components of pointsj share a common point.
  ## pointsk is
  
  trend=lm(form,data=data,y=TRUE,x=TRUE)
  Xmat=trend$x
  Y=trend$y
  
  t1=proc.time()
  opt=optim(par = init,fn = recl,gr = gr,method = "L-BFGS-B",lower = lower,
            control = list(fnscale=-1,trace=3),coordsMat=coordsMat,Y=Y,Xmat=Xmat,pointsj=pointsj,heter=heter,
            cov.model=cov.model,dS=dS)
  t2=proc.time()
  
  AB=calcAB(pointsj = pointsj,Y = Y,coordsMat = coordsMat,theta = opt$par,Xmat = Xmat,cov.model = cov.model,heter = heter)
  t3=proc.time()
  betahat=with(AB,solve(A,B))
  t4=proc.time()
  JB=calcJB(g2list = pointsjpair,pointsj = pointsj,Y = Y,coordsMat = coordsMat,theta = opt$par,Xmat = Xmat,cov.model = cov.model,heter = heter)
  t5=proc.time()
  covbeta=solve(with(AB,A%*%solve(JB,A)))
  t6=proc.time()
  structure(c(opt,list(betahat=betahat,varbetahat=diag(covbeta),
                       form=form,time=t6-t1,
                       cov.model=cov.model,
                       heter=heter,
                       pointsk=pointsk,
                       coordsMat=coordsMat,
                       Xmat=Xmat,
                       Y=Y,                       
                       optTime=t2-t1,ABtime=t3-t2,
                       JBtime=t5-t4)),class="reclsm")
  
}


summary.reclsm=function(obj){
  b=obj$betahat
  seb=sqrt(obj$varbetahat)
  z=b/seb
  p=pnorm(abs(z),0,1,lower.tail = FALSE)*2
  out=data.frame(betahat=b,"se(betahat)"=seb,Z=z,"Pr(Z>|z|)"=p)
  colnames(out)=c("betahat","se(betahat)","Z","Pr(|Z|>=z)")
  out$stars=cut(out[,4],breaks = c(-1e-5,0.001,0.01,0.05,0.1,1),labels = c("***","**","*",".",""))
  print(out,digits=3)
  print("0 *** 0.001 ** 0.01 * 0.05 . 0.1  1")
}


predict.reclsm=function(fit,newdf,newpointsk,pointsk,Nk,newcoordsMat){
  newXmat=model.matrix(fit$form,data=newdf)
  fixedpred=newXmat%*%fit$betahat
  
  Y=fit$Y
  
  cov.model=fit$cov.model
  heter=fit$heter
  
  ## eq(12) of Eidsvik 2014
  
  theta=fit$par
  beta=fit$beta
  coordsMat=fit$coordsMat
  pointsk=fit$pointsk
  Xmat=fit$Xmat
  
  empty=which(sapply(newpointsk,length)==0)
  if(length(empty)>0){newpointsk=newpointsk[-empty]}
  sizeNk=sapply(Nk,length)
  
  ## neighbourhood of k

  
  Ab0=foreach(k=1:length(newpointsk),.packages="geoR",.export="genSigmaj")%:%
    foreach(l=1:sizeNk[k])%dopar%{
      print(paste("kl:",k,l))
      newpk=newpointsk[[k]]
      pk=pointsk[[k]]
      pl=pointsk[[Nk[[k]][l]]]
      S=genSigmaj(rbind(newcoordsMat[newpk,],coordsMat[c(pk,pl),]),theta = theta,cov.model = cov.model,heter = heter)
      Q=solve(S)
      Q00=Q[1:length(newpk),1:length(newpk)]
      
      if(length(pk)>0){
        Q01=Q[1:length(newpk),length(newpk)+(1:length(pk))]
        bk=Q01%*%(Y[pk]-Xmat[pk,]%*%beta)
      }else{
        bk=0
      }
      
      if(length(pl)>0){
        Q02=Q[1:length(newpk),length(newpk)+length(pk)+(1:length(pl))]
        bl=Q02%*%(Y[pl]-Xmat[pl,]%*%beta)      
      }else{
        bl=0
      }
      
      if(!identical(dim(bk),dim(bl)))print("WARNING: dim(bk) is not the same as dim(bl)!")
      b0=-(bk+bl)
      
      list(Q00=Q00,b0=b0)
    }
  
  A=lapply(Ab0,function(x)Reduce("+",lapply(x,"[[",1)))
  b0=lapply(Ab0,function(x)Reduce("+",lapply(x,"[[",2)))
  
  randompred0=unlist((mapply(solve,A,b0,SIMPLIFY=FALSE)))
  randompred=NULL
  randompred[unlist(newpointsk)]=randompred0
  
  pred=fixedpred+randompred
  attr(pred,"fixed")=fixedpred
  attr(pred,"random")=randompred
  
  pred
  
}


rt2DF=function(rtnames,covdir="C:/Cubist1/Covariates"){
  odir=getwd()
  setwd(covdir)
  out=sapply(rtnames,function(x){
    rt=raster(paste0(x,".tif"))
    v0=getValues(rt)
    v0[complete.cases(v0)]
  }
  )
  setwd(odir)
  out
}