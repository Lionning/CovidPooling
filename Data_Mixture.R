########################
# Data and estimation
# MAJ 04/05/2020 vincent.brault@univ-grenoble-alpes.fr
#########################

rm(list=ls())

#### Param
Niter=1000

#### Loading
tab<-read.table("Default Dataset.csv",sep=",")

#### Transformation
#x=log(8*10^(14)*exp(-0.745*Ct),base=10)
Ct=function(x){
  (log(8)+(14-x)*log(10))/(0.745)
}
tab[,2]<-round(tab[,2])
tab[,1]<-Ct(tab[,1])
tab<-tab[nrow(tab):1,]
median(tab[-1,1]-tab[-nrow(tab),1])
ampl<-mean(tab[-1,1]-tab[-nrow(tab),1])
base<-c(tab[1,1]-ampl/2,(tab[-1,1]+tab[-nrow(tab),1])/2,tab[nrow(tab),1]+ampl/2)

#### Section III A

K=numeric(100)
pb<-txtProgressBar(0,100,style=3)
mu3=matrix(0,3,100)
mu4=matrix(0,4,100)
sigma3=matrix(0,3,100)
sigma4=matrix(0,4,100)
### Robustness
library(Rmixmod)
for (iter in 1:100){
  x<-unlist(sapply(1:nrow(tab),function(it){
    runif(tab[it,2],min = base[it],max=base[it+1])
  }))
  xem1<-mixmodCluster(x,2:6)
  K[iter]<-xem1@bestResult@nbCluster
  if (K[iter]==3){
    mu3[,iter]=sort(xem1@bestResult@parameters@mean)
    sigma3[,iter]=sort(unlist(xem1@bestResult@parameters@variance))
    xem3<-xem1
  }else{
    mu4[,iter]=sort(xem1@bestResult@parameters@mean)
    sigma4[,iter]=sort(unlist(xem1@bestResult@parameters@variance))
    xem4<-xem1
  }
  setTxtProgressBar(pb,iter)
}

table(K)
apply(mu3[,!apply(mu3==0,2,all)],1,sd)
apply(mu4[,!apply(mu4==0,2,all)],1,sd)
apply(sigma3[,!apply(mu3==0,2,all)],1,sd)
apply(sigma4[,!apply(mu4==0,2,all)],1,sd)

##### Representation 3 clusters
xem<-xem3
absci<-seq(min(base),max(base),by=0.01)
Mix<-rowSums(sapply(1:3,function(k){
  dnorm(absci,mean = xem@bestResult@parameters@mean[k],
        sd =sqrt(xem@bestResult@parameters@variance[[k]]))*
    xem@bestResult@parameters@proportions[k]}))
hist(x,breaks = base,ylab="Count",xlab="Simulated data",main="",col="grey")
lines(absci,Mix,col="blue",lwd=3)
col=c("orange","green","red")[rank(xem@bestResult@parameters@mean)]
for (k in 1:3){
  lines(absci,dnorm(absci,mean = xem@bestResult@parameters@mean[k],
                    sd =sqrt(xem@bestResult@parameters@variance[[k]]))*
          xem@bestResult@parameters@proportions[k],
        col=col[k],lwd=2)
}
##### Representation 4 clusters
xem<-xem4
absci<-seq(min(base),max(base),by=0.01)
Mix<-rowSums(sapply(1:4,function(k){
  dnorm(absci,mean = xem@bestResult@parameters@mean[k],
        sd =sqrt(xem@bestResult@parameters@variance[[k]]))*
    xem@bestResult@parameters@proportions[k]}))
hist(x,breaks = base,ylab="Count",xlab="Simulated data",main="",col="grey")
lines(absci,Mix,col="blue",lwd=3)
col=c("cyan","orange","green","red")[rank(xem@bestResult@parameters@mean)]
for (k in 1:4){
  lines(absci,dnorm(absci,mean = xem@bestResult@parameters@mean[k],
                    sd =sqrt(xem@bestResult@parameters@variance[[k]]))*
          xem@bestResult@parameters@proportions[k],
        col=col[k],lwd=2)
}