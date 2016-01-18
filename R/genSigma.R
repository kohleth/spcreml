#' Generate the Covariance matrix Sigma.
#'
#' \code{genSigma} generates the covariance matrix.
#'
#' @param distM Either a distance matrix (result from \code{dist}), or a matrix of distances.
#' @param theta A vector of named covariance paramter. Must (at least) contain named element  \emph{sigma2} and \emph{phi}. Other possible elements include \emph{nugget} and \emph{kappa}.
#' @param cov.model Argument that is passed to \code{geoR::cov.spatial}.
#' @param addNugget TRUE/FALSE. If true nugget is added to the covariance matrix.
#'

genSigma=function(distM,theta,cov.model,addNugget=TRUE){

  if(!all(c("sigma2","phi")%in%names(theta)))stop("argument theta must contain named element 'sigma2' and 'phi'.")
  if(inherits(distM,"matrix")){
    if(addNugget&(nrow(distM)!=ncol(distM))){
      stop("argument addNugget should only be TRUE if a symmetric distM is supplied.")
    }
  }
  if(addNugget&is.na(theta["nugget"]))stop("theta does not contain the named element 'nugget', but addNugget=TRUE.")

  # get kappa
  if(cov.model%in%c("matern","powered.exponential","cauchy")){
    kappa=theta["kappa"]
  }else{
    kappa=NA
  }

  Sigma=cov.spatial(obj = distM,
                    cov.model = cov.model,
                    cov.pars=theta[c("sigma2","phi")],
                    kappa=kappa)
  Sigma=as.matrix(Sigma)

  ## add nugget only if computing covariance within 1 set of points (instead of covariance between 2 sets of points)
  if(isSymmetric(Sigma)){
    if(inherits(distM,"dist")){
      Sigma=Sigma+diag(x=theta["nugget"]+theta["sigma2"],nrow=nrow(Sigma))
    }else{
      Sigma=Sigma+diag(x=theta["nugget"],nrow=nrow(Sigma))
    }
  }

  return(Sigma)
}
