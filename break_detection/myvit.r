## Code provided by Andrea Toreti

myViterbi<-function(x,K,Pi,med,sig) {

### code added by YB ###
if (sum(sig<0.001)>0) return(NULL)
########################

#browser()
lung<-length(x)
#la matrice delta ha in ogni riga il vettore (delta_t(1),...,delta_t(K))
delta<-matrix(nrow=lung,ncol=K)
xsi<-matrix(nrow=lung,ncol=K)
cost<-1/(sqrt(2*pi))

#inizializzazione con i log per evitare l'underflow e utilizzo di un fattore di scala 1e+1
#delta[1,1]<-(cost/sig[1])*exp(-((x[1]-med[1])^2)/(2*(sig[1])^2))
delta[1,1]<-(cost/sig[1])*exp(-((x[1]-med[1])^2)/(2*(sig[1])^2))*(1e+1)
xsi[1,1]<-0
 for (i in 2:K) {
  delta[1,i]<-0
  xsi[1,i]<-0
 }

#recursion con il log e con il fattore di scala per evitare l'underflow
 for (t in 2:lung) {
  for (k in 1:K) {
#   delta[t,k]<-max(delta[(t-1),]*Pi[,k])*(cost/sig[k])*exp(-((x[t]-med[k])^2)/(2*(sig[k])^2))
   delta[t,k]<-max(delta[(t-1),]*Pi[,k]*(cost/sig[k])*exp(-((x[t]-med[k])^2)/(2*(sig[k])^2))*(1e+1))
   xsi[t,k]<-which.max(delta[(t-1),]*Pi[,k])
  }
 }

#termination
qstar<-rep(NA,lung)
Pstar<-max(delta[lung,])
qstar[lung]<-which.max(delta[lung,])

#backtracking
 for (t in 1:(lung-1)) {
  qstar[(lung-t)]<-xsi[(lung-t+1),qstar[(lung-t+1)]]
 }

 return(seq=qstar)

}
  

