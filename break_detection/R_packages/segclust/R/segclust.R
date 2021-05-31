
Estep <- function(logdensity,phi) {
	storage.mode(logdensity)<-"double";
	storage.mode(phi)<-"double";
	.Call("sc_estep",logdensity,as.integer(dim(logdensity)[1]),as.integer(dim(logdensity)[2]),phi);
}

hybrid <- function(x,P,Kmax,lmin=1,lmax=length(x),vh=TRUE) {
  
   checkoptions = TRUE
   if ((vh==FALSE) & (lmin<=1)){
     checkoptions = FALSE
     cat("Error in hybrid : lmin must be greater than 2 when vh=FALSE","\n")     
    }

   if (Kmax>floor(length(x)/lmin)){
     checkoptions = FALSE
     cat("Error in hybrid : Kmax must be lower than [length(x)/lmin]","\n")
   }

   if (Kmax<floor(length(x)/lmax)){
     checkoptions = FALSE
     cat("Error in hybrid : Kmax must be greater than [length(x)/lmax]","\n")
   }

   if (P>Kmax){
     checkoptions = FALSE
     cat("Error in hybrid : the number of groups must be lower than the number of segments")
   }
   if (sum(is.na(x))>0){
     checkoptions = FALSE
     cat("Error in hybrid : the data must not contain missing value","\n")
   }

   if (checkoptions==TRUE){
     storage.mode(x)<-"double"
     res         = .Call("sc_hybrid",x,as.integer(P),as.integer(Kmax),as.integer(lmin),as.integer(lmax),as.logical(vh))   
     tmp         = res$Linc
     if (sum(tmp==0)>0){
	cat("The algorithm has faced convergence problems for P =", P,"\n")
	cat("and for configurations with", which(tmp==0)," segments","\n")
	tmp[tmp==0] = -Inf
	res$Linc    = tmp
     }
   }
	invisible(res)	
}

logdens <- function(x,P,phi) {
	storage.mode(x)<-"double";
	storage.mode(phi)<-"double";
	.Call("sc_logdens",x,as.integer(P),phi);
}

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
    .Call("sc_segmean",x,as.integer(lmin),as.integer(lmax),as.integer(Kmax),as.logical(vh))
  }
}

segmixt <- function(x,P,Kmax,phi,lmin=1,lmax=length(x)) {

  checkoptions = TRUE
  if (Kmax> floor(length(x)/lmin)){
    cat("Error in segmixt : Kmax must be lower than [length(x)/lmin]","\n")
    checkoptions = FALSE
  }
   
  if (Kmax<floor(length(x)/lmax)){
     checkoptions = FALSE
     cat("Error in segmixt : Kmax must be greater than [length(x)/lmax]","\n")
   }

  if (P>Kmax){
    checkoptions = FALSE
    cat("Error in segmixt : the number of groups must be lower than the number of segments","\n")
  }
  if (sum(is.na(x))>0){
    checkoptions = FALSE
    cat("Error in segmixt : the data must not contain missing values","\n")
  }
  if (checkoptions == TRUE){
    storage.mode(x)<-"double"
    storage.mode(phi)<-"double"
    .Call("sc_segmixt",x,as.integer(lmin),as.integer(lmax),as.integer(Kmax),phi, as.integer(P))
  }
}

EMalgo <- function(x,phi,rupt,P,vh=TRUE){
  checkoptions = TRUE
  K = dim(rupt)[1]
  if (P>K){
    checkoptions = FALSE
    cat("Error in EMalgo : the number of groups must be lower than the number of segments","\n")
  }
  if (checkoptions == TRUE){
    storage.mode(x)<-"double"
    storage.mode(phi) <-"double"
    storage.mode(rupt) <- "double"
    .Call("sc_EMalgo",x,phi,rupt,as.integer(K),as.integer(P),as.logical(vh))
  }
  
}


EMinit <- function(x,rupt,P,vh=TRUE){
  checkoptions = TRUE
  K = dim(rupt)[1]
  if (P>K){
    checkoptions = FALSE
    cat("Error in EMinit : the number of groups must be lower than the number of segments","\n")
  }
  if (checkoptions == TRUE){
    storage.mode(x)<-"double"
    storage.mode(rupt) <- "double"
    .Call("sc_EMinit",x,rupt,as.integer(K),as.integer(P),as.logical(vh))
  }
  
}
