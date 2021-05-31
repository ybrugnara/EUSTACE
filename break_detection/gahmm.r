## Code provided by Andrea Toreti

gahmm<-function(x,kmax){
#library(HiddenMarkov)
#library(vegan)
#library(genalg)
source('myrbga.r')
source('myevalh.r')
source('myBW.r')
source('myvit.r')

lung<-length(x)
obs<-x
kprob<-rep(NA,kmax)
kprob2<-rep(NA,kmax)
mystate<-matrix(nrow=kmax,ncol=lung)
mestat<-matrix(nrow=kmax,ncol=kmax)
sistat<-matrix(nrow=kmax,ncol=kmax)

for (K in 2:kmax) {
 p<-(lung-K)/(lung-1)

 Pi<-matrix(nrow=K,ncol=K)
 for (j in 1:K){
  for(i in 1:K){
   if(i==j){Pi[i,j]<-p}
   if(i==j-1){Pi[i,j]<-1-p}
   if(i!=j&i!=j-1) {Pi[i,j]<-0}
  }
 }
 Pi[K,K]<-1
 Pi2<-Pi%*%Pi

 delta<-c()
 delta[1]<-1
 for (i in 2:K){delta[i]<-0}


 cromin<-matrix(nrow=1,ncol=lung)
 cromax<-matrix(nrow=1,ncol=lung)
 cromin[]<-1
 cromax[]<-K

 cromax[1]<-1
 cromin[lung]<-K
 source('myeval.r')
 genetic<-myrbga(cromin,cromax,evalFunc=myeval,ydati=obs,ypi=Pi)
 mineval<-min(genetic$evaluations)
 mindex<-which(genetic$evaluations==mineval)
 mybest<-genetic$population[mindex,]


 if(length(mybest)>lung) {
  initstate<-mybest[1,]
 }
 if(length(mybest)==lung) {initstate<-mybest}

 initmedia<-c()
 initsigma<-c()
 ferma<-0
 for(ini in 1:K) {
  initind<-which(initstate==ini)
  if(length(initind)<=3){ferma<-1}
  initmedia[ini]<-mean(obs[initind])
  initsigma[ini]<-sd(obs[initind])
 }

 micontrol<-0

 initlok<-myevalh(initstate,Pi,obs,initmedia,initsigma)
 hd<-dthmm(obs,Pi,delta,pm=list(mean=initmedia,sd=initsigma),distn='norm')
 BW<-myBaumWelch(hd,control=bwcontrol(prt=FALSE,maxiter=20,posdiff=FALSE))
### code added by YB ###
 if(is.null(BW)) return(NULL)
########################
 if(!is.na(BW$LL)&all(diag(BW$Pi>=0.5))) {
  state<-myViterbi(obs,K,BW$Pi,BW$pm$mean,BW$pm$sd)
### code added by YB ###
  if(is.null(state)) return(NULL)
########################
  micontrol<-1
  mynewloko<-myevalh(state,BW$Pi,obs,BW$pm$mean,BW$pm$sd)
 }
 if(micontrol==0) {ferma<-1}
    
 if(ferma==1) {
  mynewloko<-NA
  break
 } 



 iclfactor<-(2*K+(K-1))*(log(lung)/2)



 mynewloko<-mynewloko-iclfactor
 if(is.na(mynewloko)){break}


 mystate[K,]<-state
 mestat[K,1:K]<-BW$pm$mean
 sistat[K,1:K]<-BW$pm$sd

 kprob2[K]<-mynewloko
}

newmedo<-mean(obs)
newsigo<-sd(obs)
mestat[1,1]<-newmedo
sistat[1,1]<-newsigo
omobsodd<-(obs-newmedo)^2
lokomodd<-(lung*log(1/(sqrt(2*pi)*newsigo)))-(1/(2*newsigo^2))*sum(omobsodd)
newlomo<-lokomodd-(log(lung))
kprob2[1]<-newlomo

mystate[1,]<-1

minprobI<-which.max(kprob2)

estimstateI<-mystate[minprobI,]



if(minprobI>1){
 cpI<-rep(NA,(minprobI-1)) 
 for(chap in 2:minprobI) {
  cpindex<-which(estimstateI==chap)
  cpI[(chap-1)]<-cpindex[1]-1
 }
}
else {cpI<-0}



names(kprob2)<-c(1:kmax)
altrost<-kprob2[(minprobI+1)] 
if(!is.na(altrost)){
 altros<-as.integer(names(altrost))
 if(!is.na(kprob2[altros])&(altros>minprobI)) {
  if(abs(kprob2[altros]-kprob2[minprobI])<=10){

    mrpa<-rep(NA,altros)
    mypvalue<-rep(NA,(altros-1))
    for(myrp in 1:altros) {
     mrpa[myrp]<-which(mystate[altros,]==myrp)[1]
    }

    if(altros>2){
     for(myrp in 1:(altros-2)) {
      mmrpp<-mrpp(obs[mrpa[myrp]:(mrpa[(myrp+2)]-1)],mystate[altros,mrpa[myrp]:(mrpa[(myrp+2)]-1)])
      mypvalue[myrp]<-mmrpp$Pvalue
     }
    }

    mmrpp<-mrpp(obs[mrpa[(altros-1)]:lung],mystate[altros,mrpa[(altros-1)]:lung])
    mypvalue[(altros-1)]<-mmrpp$Pvalue   

    if(all(mypvalue[]<=0.05)){
     estimstateI<-mystate[altros,]
     for(chap in 2:altros) {
      cpindex<-which(estimstateI==chap)
      cpI[(chap-1)]<-cpindex[1]-1
     }
    }

  }
 }
}

 
results<-list(all=mystate,ICLstate=estimstateI,ICLLik=kprob2,ICLChan=cpI,mean=mestat,sigma=sistat)
return(results)
}














