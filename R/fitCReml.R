#' Fit a spatial model by Composite REML
#'
#' @param form a formula that specifies the response and the fixed effect variables.
#' @param data a data.frame that contains the response variable and the fixed effect variables.
#' @param dS a list of functions to compute the first derivative of Sigma (dSigma/dtheta).
#' @param coordsMat a N by D coordinate matrix. D being the spatial dimension, typically 2 or 3.
#' @param init named vector of initial value for optimiation of \code{recl}. Must (at least) contain named elements \emph{sigma2} and \emph{phi}. Can also contain \emph{nugget} and \emph{kappa}.
#' @param poinstj A list of row numbers that indicates which row belongs to block j (block j is itself a union of block k and block l). See Details.
#' @param pointsjpair A list of pair integers from \code{1} to \code{length(pointsj)}. Each pair indicates which pairs of \code{pointjs} components should be bundled together in the variance calculation. See Details.
#' @param pointsk not used.
#' @param cov.model character string that specifies the spatial covariance model. See \code{?geoR::cov.spatial}
#' @param lower Numeric vector of lower bound for the estimating parameters. Default to NULL (no lower bound).
#'
#' @export


fitCReml=function(form,data,dS,coordsMat,init,pointsj,pointsjpair,pointsk,cov.model=cov.model,lower=NULL){

  trend=lm(form,data=data,y=TRUE,x=TRUE)
  Xmat=trend$x
  Y=trend$y

  t1=proc.time()
  opt=optim(par = init,fn = recl,gr = gr,method = "L-BFGS-B",lower = lower,
            control = list(fnscale=-1,trace=3),coordsMat=coordsMat,Y=Y,Xmat=Xmat,pointsj=pointsj,
            cov.model=cov.model,dS=dS)
  t2=proc.time()
  AB=calcAB(pointsj = pointsj,Y = Y,coordsMat = coordsMat,theta = opt$par,Xmat = Xmat,cov.model = cov.model,heter = heter)
  t3=proc.time()
  betahat=with(AB,solve(A,B))
  t4=proc.time()

  list(opt$par,betahat)
  JB=calcJB(g2list = pointsjpair,pointsj = pointsj,Y = Y,coordsMat = coordsMat,theta = opt$par,Xmat = Xmat,cov.model = cov.model)
  t5=proc.time()
  covbeta=solve(with(AB,A%*%solve(JB,A)))
  t6=proc.time()
  structure(c(opt,list(betahat=betahat,varbetahat=diag(covbeta),
                       form=form,time=t6-t1,
                       cov.model=cov.model,
                       pointsk=pointsk,
                       coordsMat=coordsMat,
                       Xmat=Xmat,
                       Y=Y,
                       optTime=t2-t1,ABtime=t3-t2,
                       JBtime=t5-t4)),class="spCReml")

}
