DPEM <- function(x,Pmin,Pmax,Kmax,method,draw,lmin=1,lmax=length(x),vh=TRUE,S=0.5){

  if (Pmin==Pmax){
    cat("use hybrid() function instead","\n")
  } else {
    
    Linc       <- matrix(-Inf, ncol=Pmax,nrow= Kmax)
    param.list <- list()
    
    for (P in (Pmin:Pmax)){   
      out.hybrid      <- hybrid(x,P,Kmax,lmin,lmax,vh)
      param.list[[P]] <- out.hybrid$param
      Linc[,P]        <- out.hybrid$Linc    
    }
    
    out.select <- segclustselect(x,param=param.list,Pmin,Pmax,Kmax,Linc,method = method,S,vh)
    output     <- segclustout(x,param.list[[out.select$Pselect]],out.select$Pselect,out.select$Kselect,draw=draw)
    invisible(output)
  }
}
