clustersegments <- function (x,P,bp,vh=TRUE){

  K          = length(which(bp==1))
  rupt       = matrix(ncol=2,nrow=K)
  rupt[,2]   = which(bp==1)
  rupt[,1]   = c(1,rupt[(2:K-1),2]+1)
  phi0       = EMinit(x,rupt,P,vh)
  phi        = EMalgo(x,phi0,rupt,P,vh)
  logdensity = t( apply(rupt,1,FUN=function(y) logdens(x[ y[1]:y[2] ] ,P,phi)))
  tau        = Estep(logdensity,phi)
  pop        = apply(tau,1,which.max)
  cluster    = c()
  tmp        = c()
  tmp2       = c()
  for ( k in (1:K) ) {
    cluster[rupt[k,1]:rupt[k,2]] = rep(pop[k],rupt[k,2]-rupt[k,1]+1)
    tmp[rupt[k,1]:rupt[k,2]]     = rep(phi[pop[k]],rupt[k,2]-rupt[k,1]+1)
    tmp2[rupt[k,1]:rupt[k,2]]    = rep(phi[pop[k]+P],rupt[k,2]-rupt[k,1]+1)
  }
  bp           = rep(0,length(x))
  bp[rupt[,2]] = 1
  output       = data.frame(signal = x, mean = tmp, sd = tmp2 ,cluster = cluster,bp=bp)
  invisible(output)
}

