segclustout <- function (x,param,P,K,draw){

  
  
  phi        = param[[K]]$phi
  rupt       = param[[K]]$rupt
  logdensity = t( apply(rupt,1,FUN=function(y) logdens(   x[ y[1]:y[2] ] ,P,phi)))
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
  bp = rep(0,length(x))
  bp[rupt[,2]] = 1
  output = data.frame(signal = x, mean = tmp, sd = tmp2 ,cluster = cluster,bp=bp)

  if (draw==TRUE){       
    phi = sort(unique(output$mean))
    par(lwd=2)
    plot(x,type="n",xlab = "genomic order",ylab = "signal")    
    for ( k in (1:K) ) {
      nk = rupt[k,2]-rupt[k,1]+1
      lines(rupt[k,1]:rupt[k,2],x[rupt[k,1]:rupt[k,2]],type="p",col=cluster[rupt[k,1]]+1, pch=cluster[rupt[k,1]])
      lines(rupt[k,1]:rupt[k,2],matrix(phi[cluster[rupt[k,1]]],nrow=1,ncol=nk),col = cluster[rupt[k,1]]+1)
    } #end k      
    abline(v= which(output$bp==1),col=1)
  }  
  invisible(output)


}


