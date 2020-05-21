########################
# Censure plan ale
# MAJ 10/05/2020 vincent.brault@univ-grenoble-alpes.fr
#########################

rm(list=ls())
setwd("/home/vbrault/Censure/")## Pour les clusters

N=10^(2:5)
S=(-2):3
Essai=10000


######################
# Warning
# Simulations take time
# We recommend to parallelize it
######################

pb<-txtProgressBar(min=0,max=length(N)*length(S)*Essai,style=3)
it=1
for (n in N){
  for (s in S){
    Nom_fichier<-paste0("Plan_ale_n_",n,"_s_",s,".Rdata")
    mu_est=numeric(Essai)
    sigma_est=numeric(Essai)
    for (iter in 1:Essai){
      Y<-sapply(1:n,function(i){
        x<-rnorm(1,0,1)
        while(x>s){
          x<-rnorm(1,0,1)
        }
        x
      })
      LogVrais<-function(theta){
        if (theta[2]<=0){
          res<-+Inf
        }else{
          res<-length(Y)*log(pnorm(s,mean = theta[1],sd=abs(theta[2])))-
            sum(log(dnorm(Y,mean = theta[1],sd=abs(theta[2]))))
        }
        res
      }
      #### Intelligent initialisation
      nlm<-try(nlm(LogVrais,c(0,1),hessian = TRUE),TRUE)
      nlm$minimum<-+Inf
      #### Random initialisation
      for (try in 1:9){
        nlm2<-try(nlm(LogVrais,c(rnorm(1,0,1),rexp(n = 1,rate = 1)),hessian = TRUE),TRUE)
        if (class(nlm2)!="try-error"){
          while(nlm2$code==4){
            nlm2<-try(nlm(LogVrais,nlm2$estimate,hessian = TRUE),TRUE)
            if (class(nlm2)=="try-error"){
              break
            }
          }
          if (class(nlm2)!="try-error"){ 
            if ((nlm2$minimum<nlm$minimum)&(nlm2$code<3)){
              nlm<-nlm2
            }
          }
        }
      }
      mu_est[iter]=nlm$estimate[1]
      sigma_est[iter]=abs(nlm$estimate[2])
      setTxtProgressBar(pb,it)
      it=it+1
    }
    save("mu_est","sigma_est",file=Nom_fichier)
  }
}