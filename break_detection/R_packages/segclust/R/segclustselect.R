segclustselect <- function (x,param,Pmin,Pmax,Kmax,Linc,method,S=0.5,lmin=1,lmax=length(x),vh=TRUE){


  n          = length(x)
  if (method =="sequential"){
    out      = segmean(x,Kmax=Kmax,lmin=lmin,lmax=lmax,vh=vh)
    Kselect  = segselect(out$J.est,Kmax,S=0.75)
    if (Kselect>1){
      Lincmax  = apply(Linc,2,max)
      J        = - Lincmax
      Jtild    = (J[Pmax]-J)/(J[Pmax]-J[Pmin])*(Pmax-Pmin)+1 
      Pselect  = max(which(diff(diff(Jtild))>S))+1
      if (is.finite(Pselect)==F) {Pselect=Pmin}
      Kselect    = which.max(Linc[,Pselect]-0.5*log(n)*c(1:Kmax))
    }
    else P=1
  }

  else if (method == "BIC"){
    
    SSall      = sum((x-mean(x))^2)
    BIC        = matrix(0,ncol=Pmax,nrow=Kmax)

    for (p in c(Pmin:Pmax)){

      if (n>100){
        a = 0.5*(n-p-1)
        b = 0.5*(n-1)       
        gamma.coef = ((a-0.5)*log(a) - a + 0.5* log(2*pi)) - ((b-0.5)*log(b) - b + 0.5* log(2*pi))
      } else{
        gamma.coef = log( gamma(0.5*(n-p-1))) - log( gamma(0.5*(n-1)) )
      }
      
      
      if (p==1){
        for (k in (1:Kmax)){         
          BIC[k,1] = gamma.coef+0.5*log(SSall) - k*log(n)
        }
      } else {

        for (k in (p:Kmax)){

          SSbg       = 0
          np         = c()
          paramtmp   = param[[p]]
          phitmp     = paramtmp[[k]]$phi
          rupt       = paramtmp[[k]]$rupt
          nk         = rupt[,2]-rupt[,1]+1
          logdensity = t( apply(rupt,1,FUN=function(y) logdens(   x[ y[1]:y[2] ] ,p,phitmp)))          
          tau        = Estep(logdensity,phitmp)
          pop        = apply(tau,1,which.max)
          cluster = c()
          for ( j in (1:k) ) {
            cluster[rupt[j,1]:rupt[j,2]] = rep(pop[j],rupt[j,2]-rupt[j,1]+1)
          } 
          for (ell in c(1:p)){
            np[ell] = sum(nk[pop==ell])
            mp = mean(x[cluster==ell])
            SSbg = SSbg+np[ell]*(mp-mean(x))^2        
          }         
          BIC[k,p] = (n-p-1)/2 * log(1+ (SSbg)/(SSall-SSbg)) + gamma.coef + 
              0.5*p*log(SSall)-0.5 * sum(log(np))-(k-0.5)*log(n)                       
        } # end k
      } # end p
      
      BIC[BIC==0] = -Inf
      BIC[is.na(BIC)] = -Inf
      Pselect = which.max(apply(BIC,2,max))
      Kselect = which.max(apply(BIC,1,max))
    }
  }# end else
    
  invisible(list(Pselect = Pselect,Kselect = Kselect))

}
