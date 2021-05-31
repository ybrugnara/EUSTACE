#Call the external Fortran subroutine myeval.f
myeval<-function(string,a,b) {
  mylik <- .Fortran("myeval",string=as.integer(string),a=as.numeric(a),b=as.numeric(b),n1=as.integer(length(b)),n2=as.integer(dim(a)[1]),mylik=0)$mylik
  return(mylik) 
}
