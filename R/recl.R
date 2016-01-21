#' Calculate the Residual Composite log-Likelihood.
#'
#' \code{recl} calculates the residual composite log-likelihood.
#'
#' @param theta A named numeric vector containing the covariance parameters. It should (at least) have element \emph{sigma2} and \emph{phi}. Other possible elements are \emph{nugget} and \emph{kappa}.
#' @param coordsMat A N by D matrix of coordinates. D is the spatial dimension, typically 2 or 3.
#' @param Y the response vector
#' @param Xmat the design matrix (not the data.frame)
#' @param pointsj A list. Each list contains the row number of the data matrix (Y or Xmat or coordsMat) that correspond to the (k,l) block.
#' @param cov.model This is passed to \code{geoR::cov.spatial}
#'
#' @export

recl=function(theta,coordsMat,Y,Xmat,pointsj,cov.model="exp",...){
  if(length(unique(c(nrow(coordsMat),length(Y),nrow(Xmat))))>1)stop("coordsMat, Y, and X must be of the same length (or have the same number of rows).")


  recloglik=foreach(pointsj=pointsj,.packages=c("geoR","Matrix"),.combine="+",.export=c("dist","genSigma"))%dopar%{

    coordsMatj=coordsMat[pointsj,]
    Sigmaj=genSigma(distM = dist(coordsMatj),theta=theta,cov.model=cov.model,addNugget=TRUE)

    Yj=Y[pointsj]
    Xj=Xmat[pointsj,]

    logdetS=determinant(Sigmaj,log=TRUE)$modulus
    SinvX=solve(Sigmaj,Xj)
    XSinvX=crossprod(Xj,solve(Sigmaj,Xj))
    logdetXSinvX=determinant(XSinvX,logarithm = TRUE)$modulus
    PYj=solve(Sigmaj,Yj)-SinvX%*%solve(XSinvX,crossprod(SinvX,Yj))

    as.numeric(-(logdetS+logdetXSinvX+crossprod(Yj,PYj)))/2
  }
  return(recloglik)
}
