segout <- function (x,K,th,draw,vh=TRUE){

  th   =  th[K,1:K]
  rupt =  matrix(ncol=2,c(c(1,th[1:K-1]+1),th))
  tmp  = c()
  tmp2 = c()
  for ( k in (1:K) ) {
    mk = mean(x[rupt[k,1]:rupt[k,2]])
    sk = sd(x[rupt[k,1]:rupt[k,2]])
    tmp[rupt[k,1]:rupt[k,2]]     = rep(mk,rupt[k,2]-rupt[k,1]+1)
    tmp2[rupt[k,1]:rupt[k,2]]    = rep(sk,rupt[k,2]-rupt[k,1]+1)
  }
  bp = rep(0,length(x))
  bp[rupt[,2]] = 1
  if (vh == TRUE){
    sd = mean((x-tmp)^2)
  } else {
    sd = tmp2
  }
  output = data.frame(signal = x , mean = tmp, sd= sd, bp = bp)

  if (draw==TRUE){
    par(lwd=2)
    plot(x,xlab = "genomic order",ylab = "signal")    
    for ( k in (1:K) ) {
      nk = rupt[k,2]-rupt[k,1]+1
      lines(rupt[k,1]:rupt[k,2],rep(mean(x[rupt[k,1]:rupt[k,2]]),nk),col=2)
    } #end k
    abline(v= which(bp==1),col=1)
  }


  invisible(output)
  
}


