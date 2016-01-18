#' calcAB, calcJB
#'

calcAB=function(pointsj,Y,coordsMat,theta,Xmat,cov.model,heter,...){
  out=foreach(pointsj=pointsj,.packages="geoR",.export=c("genSigma"))%dopar%{
    A=B=0
    Sigmaj=genSigma(distM=dist(coordsMat[pointsj,]),theta = theta,cov.model = cov.model,addNugget=TRUE)
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
  out=foreach(g2=g2list,.packages="geoR",.export=c("genSigma"),.combine="+")%dopar%{
    #   g2=g2list[[1]]
    j=pointsj[[g2[[1]]]]
    j2=pointsj[[g2[[2]]]]
    Xj=Xmat[j,]
    Xj2=Xmat[j2,]
    Sj=genSigma(distM = dist(coordsMat[j,]),theta = theta,cov.model = cov.model,addNugget = TRUE)
    Sj2=genSigma(distM = dist(coordsMat[j2,]),theta = theta,cov.model = cov.model,addNugget = TRUE)
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
