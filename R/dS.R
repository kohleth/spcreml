#' Derivative of Distributions functions
#'


# common nugget
dSdtausq=function(D,theta0,...)diag(1,nrow=nrow(D))

# exponential
dSexpds2=function(D,theta0,...)exp(-D/theta0[2])
dSexpdphi=function(D,theta0,...)D*theta0[1]/theta0[2]^2*exp(-D/theta0[2])
dSexp=list(dSexpds2,dSexpdphi,dSdtausq)

# spherical
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
dSsph=list(dSsphds2,dSsphdphi,dSdtausq)

# gaussian
dSgauds2=function(D,theta0,...)exp(-(D/theta0[2])^2)
dSgaudphi=function(D,theta0,...)2*theta0[1]^2*D^2/(theta0[2]^3)*exp(-(D/theta0[2])^2)
dSgau=list(dSgauds2,dSgaudphi,dSdtausq)

#
