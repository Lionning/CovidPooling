########################
# Censure EM
# MAJ 04/05/2020 vincent.brault@univ-grenoble-alpes.fr
#########################

rm(list=ls())

#### Param
Niter=100

#### Chargement
tab<-read.table("Default Dataset.csv",sep=",")

#### Transfo base
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

#### Init avec Rmixmod
library(Rmixmod)
x<-unlist(sapply(1:nrow(tab),function(it){
  runif(tab[it,2],min = base[it],max=base[it+1])
}))
xem1<-mixmodCluster(x,3)
hist(xem1)


for (Nb in 1:3){
  
  #tab[nrow(tab)-(0:Nb),2]=0
  s=base[nrow(tab)-(Nb)]
  
  #### Sample
  x<-unlist(sapply(1:nrow(tab),function(it){
    runif(tab[it,2],min = base[it],max=base[it+1])
  }))
  hist(x,breaks = base,probability = FALSE)
  points(tab[,1],tab[,2],col="red",pch=19)
  abline(v=s,col="blue",lty=2,lwd=2)
  
  
  #### EM censure
  pik<-xem1@bestResult@parameters@proportions
  muk<-xem1@bestResult@parameters@mean
  sigmak<-sqrt(unlist(xem1@bestResult@parameters@variance))
  K=length(pik)
  pk<-rep(0.1,K)
  
  for (it in 1:Niter){
    ## Step E
    # log
    lsik<-sapply(1:K,function(k){
      sapply(1:length(x),function(i){
        -log(pk[k]+(1-pk[k])*pnorm((s-muk[k])/sigmak[k]))-
          log(sigmak[k]^2)/2-(x[i]-muk[k])^2/(2*sigmak[k]^2)+log(pik[k])
      })
    })
    # remove max + exp
    sik<-exp(lsik-apply(lsik,1,max)%*%matrix(1,nrow=1,ncol=K))
    # scale
    sik<-sik/apply(sik,1,sum)%*%matrix(1,nrow=1,ncol=K)
    
    ## Step M
    pik<-apply(sik,2,mean)
    for (k in 1:K){
      ## Maximisation with nlm
      # Fonction
      LogVraisk<-function(theta){
        if ((theta[2]<0)||(theta[3]<0)||(theta[3]>=1)){
          res<-+Inf
        }else{
        res<-sum(sik[,k])*(log(theta[3]+(1-theta[3])*pnorm((s-theta[1])/abs(theta[2])))+
                             log(2*pi*theta[2]^2)/2)+sum(sik[,k]*(x-theta[1])^2)/(2*theta[2]^2)-
          sum(sik[x>s,k])*log(theta[3])
        }
        res
      }
      # Optim
      nlm<-nlm(LogVraisk,c(muk[k],sigmak[k],pk[k]),hessian = TRUE)
      muk[k]<-nlm$estimate[1]
      sigmak[k]<-abs(nlm$estimate[2])
      pk[k]<-abs(nlm$estimate[3])
    }
  }
  
  
  ### Display
  
  absci<-seq(min(base),max(base),by=0.01)
  hist(x,breaks = base,ylab="Count",xlab="Simulated data",main="")
  Mix<-rowSums(sapply(1:3,function(k){
    dnorm(absci,mean = muk[k],
          sd =sigmak[k])*
      pik[k]}))
  lines(absci,Mix,col="grey",lwd=3)
  col=c("red","green","blue")[rank(muk)]
  for (k in 1:3){
    lines(absci,dnorm(absci,mean = muk[k],
                      sd =sigmak[k])*
            pik[k],
          col=col[k],lty=2,lwd=2)
  }
  
  cat("s=",s,"\n")
  print(muk[order(muk)])
  print(sigmak[order(muk)])
  print(pik[order(muk)])
  print(pk[order(muk)])
}