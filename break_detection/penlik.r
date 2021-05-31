#Caussinus and Mestre method for break point detection

segmean <- function(x,Kmax,lmin=1,lmax=length(x),vh=TRUE) {
  checkoptions = TRUE
  if ((vh==FALSE) & (lmin<=1)){
    checkoptions = FALSE
    cat("Error in segmean : lmin must be greater than 2 when vh=FALSE","\n")
  }
  if (Kmax> floor(length(x)/lmin)){
    checkoptions = FALSE
    cat("Error in segmean : Kmax must be lower than [length(x)/lmin]","\n")
  }
  if (Kmax<floor(length(x)/lmax)){
    checkoptions = FALSE
    cat("Error in segmean : Kmax must be greater than [length(x)/lmax]","\n")
  }
  if (sum(is.na(x))>0){
    checkoptions = FALSE
    cat("Error in segmean : the data must not contain missing values","\n")
  }
  if (checkoptions == TRUE){
    storage.mode(x)<-"double"
    .Call("sc_segmean",x,as.integer(lmin),as.integer(lmax),as.integer(Kmax),as.logical(vh),
          PACKAGE = "segclust")
  }
}

pen.lik <- function(y,kmax) {
           seg.out   = segmean(y,kmax+1,lmin=3)
           pen.lik   = log(seg.out$J.est[2:kmax+1]/seg.out$J.est[1])
           for (k in 1:(kmax-1))
               {pen.lik[k]=pen.lik[k]+2*k*log(length(y))/(length(y)-1)}
           if (min(pen.lik) >= 0) {brk=0}
           else {kopt=which.min(pen.lik)
                 brk=seg.out$t.est[kopt+1,j=1:kopt]}
           return(brk)
}


