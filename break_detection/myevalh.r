## Code provided by Andrea Toreti

#scrivo la funzione di valutazione del GA vedi Kehagia 2003 (28)
myevalh<-function(string,a,b,media,sigma) {
#browser()
obs<-b
Pi<-a
lung<-length(obs)
K<-length(Pi[,1])
p<-rep(NA,K)
 for(indi in 1:K) {
  p[indi]<-Pi[indi,indi]
  }
ntilde<-rep(NA,K)
lolikd<-rep(NA,K)
    for (ntil in 1:K) {
     nindex<-which(string==ntil)
     ntilde[ntil]<-length(nindex)
     if(ntil!=K){lolikd[ntil]<-(log(p[ntil]))*(ntilde[ntil]-1)}  
    }
lolikd[K]<-0
#controllo che la mia sequenza abbia tutti gli stati
control<-c()
checkpoint<-0
 for (l in 1:K){
  control<-which(string[]==l)
  if(length(control)==0) {checkpoint<-1;break}
 }
 if(checkpoint==0){
 stdobs<-matrix(nrow=1,ncol=lung)
 sigobs<-matrix(nrow=1,ncol=lung)
   for(j in 1:K) {
    index<-which(string[]==j)
     if(length(index)>0){
      stdobs[index]<-((obs[index]-media[j])^2)/(2*sigma[j]*sigma[j])
      sigobs[index]<-log(sigma[j]*sqrt(2*pi))
     }
   }
 lolika<-sum(stdobs)
 lolikb<-sum(sigobs)
 lolikc<-sum(log(1-p[1:(K-1)]))
 lolikdd<-sum(lolikd)
 loglik<-(-lolika)-lolikb+lolikc+lolikdd
 }
 if(checkpoint!=0) {loglik<--1000}
 return(loglik)
}
