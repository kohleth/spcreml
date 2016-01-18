#' Summary for spCReml object

summary.spCReml=function(obj){
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
