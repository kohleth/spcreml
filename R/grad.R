#' Computes the grad and hessian of \code{recl}.
#'
#' @param theta named vector of covariance parameter. Must (at least) contain named elements \emph{sigma2} and \emph{phi}. Can also contain \emph{nugget} and \emph{kappa}.
#' @param dS a list of functions to compute the first derivative of Sigma (dSigma/dtheta).
#' @param coordsMat a N by D coordinate matrix. D being the spatial dimension, typically 2 or 3.
#' @param Y the response vector
#' @param Xmat the design matrix (not the data.frame)
#' @param poinstj A list of row numbers that indicates which row belongs to block j (block j is itself a union of block k and block l). See Details.
#' @param cov.model character string that specifies the spatial covariance model. See \code{?geoR::cov.spatial}
#' @param gr TRUE/FALSE. Whether the grad vector should be returned.
#' @param hs TRUE/FALSE. Whether the Hessian matrix should be returned.

grHs=function(theta,dS,coordsMat,Y,Xmat,pointsj,cov.model="exp",gr=TRUE,hs=TRUE,...){
  R=length(dS)
  out=foreach(pointsj=pointsj,.packages=c("geoR"),.export=c("dist","traceWW","genSigma"))%dopar%{
    u=rep(0,R)
    H=matrix(0,nrow=R,ncol=R)
    coordsMatj=coordsMat[pointsj,]
    Dj=as.matrix(dist(coordsMatj))

    Sigmaj=genSigma(distM = Dj,theta=theta,cov.model=cov.model,addNugget = TRUE)
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

gr=function(theta,coordsMat,Y,Xmat,dS,pointsj,cov.model,...){
  grHs(theta = theta,dS = dS,coordsMat = coordsMat,Y = Y,Xmat = Xmat,pointsj = pointsj,cov.model = cov.model,gr=TRUE,hs=FALSE)
}


Hs=function(theta,coordsMat,Y,Xmat,dS,pointsj,cov.model,...){
  grHs(theta = theta,dS = dS,coordsMat = coordsMat,Y = Y,Xmat = Xmat,pointsj = pointsj,cov.model = cov.model,gr=FALSE,hs=TRUE)
}

traceWW=function(A,B){
  stopifnot(ncol(A)==nrow(B))
  sum(sapply(1:nrow(A),function(ii)A[ii,]%*%B[,ii]))
}
