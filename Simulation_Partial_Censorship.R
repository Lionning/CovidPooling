########################
# Censure partially
# MAJ 10/05/2020 vincent.brault@univ-grenoble-alpes.fr
#########################

rm(list=ls())
setwd("/home/vbrault/Censure/")## Pour les clusters

N=10^(2:5)
S=(-2):3
Essai=10000
P<-c(0.1,0.5,0.9)

######################
# Warning
# Simulations take time
# We recommend to parallelize it
######################


pb<-txtProgressBar(min=0,max=length(N)*length(S)*length(P)*Essai,style=3)
it=1
for (p in P){
  for (n in N){
    for (s in S){
      Nom_fichier<-paste0("Plan_p_n_",n,"_s_",s,"_p_",10*p,".Rdata")
      mu_est=numeric(Essai)
      sigma_est=numeric(Essai)
      p_est=numeric(Essai)
      for (iter in 1:Essai){
        Y<-sapply(1:n,function(i){
          x<-rnorm(1,0,1)
          while((x>s)&(runif(1)>p)){
            x<-rnorm(1,0,1)
          }
          x
        })
        LogVrais<-function(theta){
          if ((theta[2]<=0)||(theta[3]<0)||(theta[3]>1)){
            res<-+Inf
          }else{
            res<-length(Y)*(log(theta[3]+(1-theta[3])*pnorm((s-theta[1])/abs(theta[2])))+
                              log(2*pi*theta[2]^2)/2)+sum((Y-theta[1])^2)/(2*theta[2]^2)-
              sum(Y>s)*log(theta[3])
          }
          res
        }
        ### Intelligence initialization
        nlm<-try(nlm(LogVrais,c(0,1,p),hessian = TRUE),TRUE)
        if (class(nlm)!="try-error"){
          if (nlm$code>2){
            class(nlm)="try-error"
          }
        }
        while(class(nlm)=="try-error"){
          nlm<-try(nlm(LogVrais,c(rnorm(1,0,1),rexp(n = 1,rate = 1),runif(1)),hessian = TRUE),TRUE)
          if (class(nlm)!="try-error"){
            if (nlm$code>2){
              class(nlm)="try-error"
            }
          }
        }
        mu_est[iter]=nlm$estimate[1]
        sigma_est[iter]=abs(nlm$estimate[2])
        p_est[iter]=nlm$estimate[3]
        setTxtProgressBar(pb,it)
        it=it+1
      }
      save("mu_est","sigma_est","p_est",file=Nom_fichier)
    }
  }
}