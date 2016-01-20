## modified splot-----
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


## form neighbour------------
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
