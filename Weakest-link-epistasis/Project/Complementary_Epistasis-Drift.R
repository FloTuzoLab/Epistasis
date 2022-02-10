#1.Kimura fixation probability
u_exact<-function(s,Ne,p){
  return((1-exp(-2*Ne*s*p))/(1-exp(-2*Ne*s)))
}

#1.Toy model 1: drawing fitness values directly from distributions
##a.Mutations drawn from beta distributions
b_ds_set=c(-5*10^-1,-10^-1)
beta_set=c(1,2)
bias_leg_set=c("High bias","Low bias")
centered_distribution_set=c("decentered distribution","centered distribution")

par(mfrow=c(2,2),mai=c(0.4,0.8,0.4,0.8))
for(i in 1:2){
  for(j in 1:2){
    bias_leg=bias_leg_set[i]
    centered_distribution=centered_distribution_set[j]
    b_ds=b_ds_set[i]
    beta=beta_set[j]
    curve(dbeta(x,0.5*(1+b_ds)/(1-0.5*(1+b_ds))*beta,beta),col = "blue",main=paste(bias_leg,"&",centered_distribution), xlab = "Fitness value", ylab = "Density", xlim=c(0,1),ylim=c(0,10),cex.main=0.75)
    curve(dbeta(x,1*(0.9+b_ds)/(1-0.9*(1+b_ds))*beta,beta),col = "gold",add=TRUE)
    curve(dbeta(x,1*(1+b_ds)/(1-1*(1+b_ds))*beta,beta),col = "red",add=TRUE)
    legend("top",legend=c("f=0.5","f=0.9","f=1"),lty=1,col=c("blue","gold","red"),cex=0.75)
  }
}

##b.Simulations for small populations10
N_pop_set=c(10)
N_pop=N_pop_set[1]
beta_set=c(1,2)
N_mod_set=c(1,3,9,27)
N_rep=30
bias<-c(rep(rep(b_ds_set[1],length(beta_set)*length(N_mod_set)),N_rep),rep(rep(b_ds_set[2],length(beta_set)*length(N_mod_set)),N_rep))
dist_mut<-rep(c(rep(c(rep(beta_set[1],length(N_mod_set))),N_rep),rep(c(rep(beta_set[2],length(N_mod_set))),N_rep)),length(beta_set))
N_mod_data<-rep(c(rep(N_mod_set[1],N_rep),rep(N_mod_set[2],N_rep),rep(N_mod_set[3],N_rep),rep(N_mod_set[4],N_rep)),length(N_mod_set))
data<-cbind(bias,dist_mut,N_mod_data)
fit_results<-c()
mod_fit_results<-list()
i=0      
for(b in 1:length(b_ds_set)){
  b_ds=b_ds_set[b]
  for(be in 1:length(beta_set)){
    beta=beta_set[be]
    for(m in 1:length(N_mod_set)){
      N_mod=N_mod_set[m]
      N_mut=N_mod*N_pop*10^2
      for (r in 1:N_rep){
        #print(paste(b,be,m,r))
        i=i+1
        f=1-1/N_pop
        #f=0.1
        f_tab<-rep(f,N_mod)
        for (t in 1:N_mut){
          pos_mut<-trunc(runif(1,1,N_mod+1))
          #print(pos_mut)
          f<-min(f_tab)
          f_curr<-rbeta(1,shape1=f_tab[pos_mut]*(1+b_ds)/(1-f_tab[pos_mut]*(1+b_ds))*beta,shape2=beta)
          #print(f_curr)
          if(f_curr==f_tab[pos_mut]){
          }
          else if(f_curr<10^-5){
            p_fix=0
          }
          else if(f_tab[pos_mut]==f){
            sec_min_fit<-min(f_tab[-which(f_tab==min(f_tab))])
            if(f_curr>sec_min_fit){
              p_fix=u_exact((sec_min_fit-f)/f,N_pop,1/N_pop)
            }
            else{
              p_fix=u_exact((f_curr-f)/f,N_pop,1/N_pop)
            }
          }
          else{
            if(f_curr>f){
              p_fix=1/N_pop
            }
            else{
              p_fix=u_exact((f_curr-f)/f,N_pop,1/N_pop)
            }
          }
          #print(p_fix)
          mut_fix=rbinom(1,1,p_fix)
          #print(mut_fix)
          if(mut_fix==1){
            #print("YES")
            f_tab[pos_mut]=f_curr
          }
        }
        fit_results[i]=f
        print(f)
        mod_fit_results[[i]]=f_tab
      }
    }
  }
}
data<-cbind(data,fit_results)
data<-data.frame(data)
#par(mfrow=c(1,1))
#boxplot(data$fit_results)

par(mfrow=c(2,2))
for(i in 1:2){
  for(j in 1:2){
    bias_leg=bias_leg_set[i]
    centered_distribution=centered_distribution_set[j]
    b_ds=b_ds_set[i]
    #print(b_ds)
    beta=beta_set[j]
    #print(beta)
    boxplot(data$fit_results[data$bias==b_ds & data$dist_mut==beta]~data$N_mod_data[data$bias==b_ds & data$dist_mut==beta],ylim=c(0,1), xlab = "Number of complementary modules", ylab = "Fitness",main=paste(bias_leg,"&",centered_distribution),cex.main=0.8)
    abline(h=1-1/N_pop,col="red")
  }
}


##c.Simulations for small populations100
N_pop_set=c(100)
N_pop=N_pop_set[1]
beta_set=c(1,2)
N_mod_set=c(1,3,9,27)
N_rep=30
bias<-c(rep(rep(b_ds_set[1],length(beta_set)*length(N_mod_set)),N_rep),rep(rep(b_ds_set[2],length(beta_set)*length(N_mod_set)),N_rep))
dist_mut<-rep(c(rep(c(rep(beta_set[1],length(N_mod_set))),N_rep),rep(c(rep(beta_set[2],length(N_mod_set))),N_rep)),length(beta_set))
N_mod_data<-rep(c(rep(N_mod_set[1],N_rep),rep(N_mod_set[2],N_rep),rep(N_mod_set[3],N_rep),rep(N_mod_set[4],N_rep)),length(N_mod_set))
data<-cbind(bias,dist_mut,N_mod_data)
fit_results<-c()
mod_fit_results<-list()
i=0      
for(b in 1:length(b_ds_set)){
  b_ds=b_ds_set[b]
  for(be in 1:length(beta_set)){
    beta=beta_set[be]
    for(m in 1:length(N_mod_set)){
      N_mod=N_mod_set[m]
      N_mut=N_mod*N_pop*10^2
      for (r in 1:N_rep){
        print(paste(b,be,m,r))
        i=i+1
        f=1-1/N_pop
        #f=0.1
        f_tab<-rep(f,N_mod)
        for (t in 1:N_mut){
          pos_mut<-trunc(runif(1,1,N_mod+1))
          #print(pos_mut)
          f<-min(f_tab)
          f_curr<-rbeta(1,shape1=f_tab[pos_mut]*(1+b_ds)/(1-f_tab[pos_mut]*(1+b_ds))*beta,shape2=beta)
          #print(f_curr)
          if(f_curr==f_tab[pos_mut]){
          }
          else if(f_curr<10^-5){
            p_fix=0
          }
          else if(f_tab[pos_mut]==f){
            sec_min_fit<-min(f_tab[-which(f_tab==min(f_tab))])
            if(f_curr>sec_min_fit){
              p_fix=u_exact((sec_min_fit-f)/f,N_pop,1/N_pop)
            }
            else{
              p_fix=u_exact((f_curr-f)/f,N_pop,1/N_pop)
            }
          }
          else{
            if(f_curr>f){
              p_fix=1/N_pop
            }
            else{
              p_fix=u_exact((f_curr-f)/f,N_pop,1/N_pop)
            }
          }
          #print(p_fix)
          mut_fix=rbinom(1,1,p_fix)
          #print(mut_fix)
          if(mut_fix==1){
            #print("YES")
            f_tab[pos_mut]=f_curr
          }
        }
        fit_results[i]=f
        print(f)
        mod_fit_results[[i]]=f_tab
      }
    }
  }
}
data<-cbind(data,fit_results)
data<-data.frame(data)
#par(mfrow=c(1,1))
#boxplot(data$fit_results)

par(mfrow=c(2,2))
for(i in 1:2){
  for(j in 1:2){
    bias_leg=bias_leg_set[i]
    centered_distribution=centered_distribution_set[j]
    b_ds=b_ds_set[i]
    #print(b_ds)
    beta=beta_set[j]
    #print(beta)
    boxplot(data$fit_results[data$bias==b_ds & data$dist_mut==beta]~data$N_mod_data[data$bias==b_ds & data$dist_mut==beta],ylim=c(0.7,1), xlab = "Number of complementary modules", ylab = "Fitness",main=paste(bias_leg,"&",centered_distribution),cex.main=0.8)
    abline(h=1-1/N_pop,col="red")
    abline(h=(1-1/(2*N_pop))^3,col="orange")
    abline(h=(1-1/(2*N_pop))^9,col="yellow")
    abline(h=(1-1/(2*N_pop))^27,col="light green")
  }
}

#2.Toy model 2: global epistasis framework
#a.Fitness function
fitness<-function(x,K_x){
  return(x/(x+K_x))
}
#b.Mutational landscape
mu_set=c(1,10)
sigmaX_set=c(0.1,0.2)
par(mfrow=c(2,2),mai=c(0.4,0.8,0.4,0.8))
for(i in 1:2){
  for(j in 1:2){
    mu=mu_set[i]
    sigma_X=mu*sigmaX_set[j]
    curve(dnorm(x,mean=mu,sd=sigma_X),col = "blue",main=paste(bias_leg,"&",centered_distribution), xlab = "Fitness value", ylab = "Density", xlim=c(0,20),ylim=c(0,1),cex.main=0.75)
    curve(dnorm(x,mean=mu-0.1,sd=sigma_X),col = "gold",add=TRUE)
    curve(dnorm(x,mean=mu-0.5,sd=sigma_X),col = "red",add=TRUE)
    legend("top",legend=c("f=0.5","f=0.9","f=1"),lty=1,col=c("blue","gold","red"),cex=0.75)
  }
}

mu_set=c(1,2,3)
sigmaX_set=c(0.5,1)
par(mfrow=c(3,2),mai=c(0.4,0.8,0.4,0.8))
for(i in 1:3){
  for(j in 1:2){
    mu=mu_set[i]
    sigma_X=sigmaX_set[j]
    print(paste(mu,sigma_X))
    curve(dnorm(x,mean=mu,sd=sigma_X),col = "blue",main=paste(bias_leg,"&",centered_distribution), xlab = "Fitness value", ylab = "Density", xlim=c(-3,5),ylim=c(0,2),cex.main=0.75)
    #curve(dnorm(x,mean=mu-0.25,sd=sigma_X),col = "gold",add=TRUE)
    curve(dnorm(x,mean=mu-1*sigma_X,sd=sigma_X),col = "orange",add=TRUE)
    curve(dnorm(x,mean=mu-2*sigma_X,sd=sigma_X),col = "red",add=TRUE)
    legend("top",legend=c("f=0.5","f=0.9","f=1"),lty=1,col=c("blue","gold","red"),cex=0.75)
    abline(v=mu,lty=2)
  }
}

Kx<-1
N_pop_set=c(10,100)
N_pop=N_pop_set[1]
sigmaX_set=c(0.25,0.5)
bias_set=c(0,1,2)
N_mod_set=c(1,2,10,50)
N_rep=30
#bias_tab<-c(rep(rep(bias_set[1],length(sigmaX_set)*length(N_mod_set)),N_rep),rep(rep(b_ds_set[2],length(beta_set)*length(N_mod_set)),N_rep))
bias_tab<-c()
for (i in 1:length(bias_set)){
  bias_tab<-c(bias_tab,rep(rep(bias_set[i],length(sigmaX_set)*length(N_mod_set)),N_rep))
}
sigma_tab<-c()
for (i in 1:length(sigmaX_set)){
  sigma_tab<-c(sigma_tab,rep(c(rep(sigmaX_set[i],length(N_mod_set))),N_rep))
}
sigma_tab<-rep(sigma_tab,length(bias_set))
#dist_mut<-rep(c(rep(c(rep(beta_set[1],length(N_mod_set))),N_rep),rep(c(rep(beta_set[2],length(N_mod_set))),N_rep)),length(beta_set))
N_mod_data<-c()
for (i in 1:length(N_mod_set)){
    N_mod_data<-c(N_mod_data,rep(N_mod_set[i],N_rep))
  }
N_mod_data<-rep(N_mod_data,length(bias_set)*length(sigmaX_set))
data<-cbind(bias_tab,sigma_tab,N_mod_data)
fit_results<-c()
mod_fit_results<-list()
i=0      
for(b in 1:length(bias_set)){
  bias=bias_set[b]
  for(s in 1:length(sigmaX_set)){
    sigma_X=sigmaX_set[s]
    for(m in 1:length(N_mod_set)){
      N_mod=N_mod_set[m]
      N_mut=N_mod*N_pop*10^3
      for (r in 1:N_rep){
        print(paste(b,s,m,r))
        i=i+1
        f0=1-1/N_pop
        x0<-f0/(1-f0)*Kx
        X_tab<-rep(x0,N_mod)
        f_tab<-fitness(X_tab,Kx)
        f<-min(f_tab)
        for (t in 1:N_mut){
          pos_mut<-trunc(runif(1,1,N_mod+1))
          #xmin<-min(X_tab)
          x_curr<-10^(rnorm(1,mean=log10(X_tab[pos_mut])-bias*sigma_X,sd=sigma_X))
          #print(f_curr)
          f_curr<-fitness(x_curr,Kx)
          if(n==1){
            if(f_curr==f_tab[pos_mut]){
              p_fix=1/N_pop
            }
            else if(log10(x_curr)>10){
              p_fix=1/N_pop
            }
            else{
              p_fix=u_exact((f_curr-f)/f,N_pop,1/N_pop)
            }
          }
          else{
            if(f_curr==f_tab[pos_mut]){
              p_fix=1/N_pop
            }
            else if(log10(x_curr)>10){
              p_fix=1/N_pop
            }
            else if(f_tab[pos_mut]==f){
              sec_min_fit<-min(f_tab[-which(f_tab==min(f_tab))])
              if(f_curr>sec_min_fit){
                p_fix=u_exact((sec_min_fit-f)/f,N_pop,1/N_pop)
              }
              else{
                p_fix=u_exact((f_curr-f)/f,N_pop,1/N_pop)
              }
            }
            else{
              if(f_curr>f){
                p_fix=1/N_pop
              }
              else{
                p_fix=u_exact((f_curr-f)/f,N_pop,1/N_pop)
              }
            }
          }
          #print(p_fix)
          mut_fix=rbinom(1,1,p_fix)
          #print(mut_fix)
          if(mut_fix==1){
            #print("YES")
            X_tab[pos_mut]=x_curr
          }
          f_tab<-fitness(X_tab,Kx)
          f<-min(f_tab)
        }
        fit_results[i]=f
        print(f)
        mod_fit_results[[i]]=f_tab
      }
    }
  }
}
data<-cbind(data,fit_results)
data<-data.frame(data)
addtxt<-list(l=4.5,h=1.05,txt=c("A","B","C","D","E","F","G","H","I"),srt = 0,font=2,col="black")
#par(mfrow=c(1,1))
#boxplot(data$fit_results)
setwd(dir="~")
setwd(dir="Desktop/Complementary-epistasis/Figures")
mut_bias_set<-c("No","Weak","High")
sd_set<-c("Weak","High")
par(mfrow=c(3,2),mai=c(0.7,0.8,0.4,0.8))
for(b in 1:length(bias_set)){
  bias=bias_set[b]
  for(s in 1:length(sigmaX_set)){
    sigma_X=sigmaX_set[s]
    #print(beta)
    boxplot(data$fit_results[data$bias==bias & data$sigma_tab==sigma_X]~data$N_mod_data[data$bias==bias & data$sigma_tab==sigma_X],ylim=c(0,1.05), xlab = "Number of complementary modules", ylab = "Fitness",cex.main=1)
    abline(h=1-1/N_pop,col="red")
    abline(h=(1-1/(N_pop))^2,col="orange")
    abline(h=(1-1/(N_pop))^10,col="yellow")
    abline(h=(1-1/(N_pop))^50,col="blue")
    text(addtxt$l,addtxt$h,addtxt$txt[s+(b-1)*length(sigmaX_set)],srt=addtxt$srt,font=addtxt$font,col=addtxt$col,cex=1.25)
  }
}
dev.print(device = jpeg, file = "Evo_Outcomes_Ne10.jpeg", width = 780*3,height=960*3,res=300,type="cairo")
#main=paste(mut_bias_set[b],"mutational bias","&",sd_set[s],"mutational variability")

N_pop=N_pop_set[2]
fit_results_Ne100<-c()
mod_fit_results_Ne100<-list()
i=0      
for(b in 1:length(bias_set)){
  bias=bias_set[b]
  for(s in 1:length(sigmaX_set)){
    sigma_X=sigmaX_set[s]
    for(m in 1:length(N_mod_set)){
      N_mod=N_mod_set[m]
      N_mut=N_mod*N_pop*10^3
      for (r in 1:N_rep){
        print(paste(b,s,m,r))
        i=i+1
        f0=1/2
        x0<-f0/(1-f0)*Kx
        X_tab<-rep(x0,N_mod)
        f_tab<-fitness(X_tab,Kx)
        f<-min(f_tab)
        for (t in 1:N_mut){
          pos_mut<-trunc(runif(1,1,N_mod+1))
          #xmin<-min(X_tab)
          x_curr<-10^(rnorm(1,mean=log10(X_tab[pos_mut])-bias*sigma_X,sd=sigma_X))
          #print(f_curr)
          f_curr<-fitness(x_curr,Kx)
          if(n==1){
            if(f_curr==f_tab[pos_mut]){
              p_fix=1/N_pop
            }
            else if(log10(x_curr)>10){
              p_fix=1/N_pop
            }
            else{
              p_fix=u_exact((f_curr-f)/f,N_pop,1/N_pop)
            }
          }
          else{
            if(f_curr==f_tab[pos_mut]){
              p_fix=1/N_pop
            }
            else if(log10(x_curr)>10){
              p_fix=1/N_pop
            }
            else if(f_tab[pos_mut]==f){
              sec_min_fit<-min(f_tab[-which(f_tab==min(f_tab))])
              if(f_curr>sec_min_fit){
                p_fix=u_exact((sec_min_fit-f)/f,N_pop,1/N_pop)
              }
              else{
                p_fix=u_exact((f_curr-f)/f,N_pop,1/N_pop)
              }
            }
            else{
              if(f_curr>f){
                p_fix=1/N_pop
              }
              else{
                p_fix=u_exact((f_curr-f)/f,N_pop,1/N_pop)
              }
            }
          }
          #print(p_fix)
          mut_fix=rbinom(1,1,p_fix)
          #print(mut_fix)
          if(mut_fix==1){
            #print("YES")
            X_tab[pos_mut]=x_curr
          }
          f_tab<-fitness(X_tab,Kx)
          f<-min(f_tab)
        }
        fit_results_Ne100[i]=f
        print(f)
        mod_fit_results_Ne100[[i]]=f_tab
      }
    }
  }
}
data<-cbind(data,fit_results_Ne100)
data<-data.frame(data)
#par(mfrow=c(1,1))
#boxplot(data$fit_results)

setwd(dir="~")
setwd(dir="Desktop/Complementary-epistasis/Figures")
mut_bias_set<-c("No","Weak","High")
sd_set<-c("Weak","High")
par(mfrow=c(3,2),mai=c(0.7,0.8,0.4,0.8))
for(b in 1:length(bias_set)){
  bias=bias_set[b]
  for(s in 1:length(sigmaX_set)){
    sigma_X=sigmaX_set[s]
    #print(beta)
    boxplot(data$fit_results_Ne100[data$bias==bias & data$sigma_tab==sigma_X]~data$N_mod_data[data$bias==bias & data$sigma_tab==sigma_X],ylim=c(0,1.05), xlab = "Number of complementary modules", ylab = "Fitness",cex.main=1)
    abline(h=1-1/N_pop,col="red")
    abline(h=(1-1/(N_pop))^2,col="orange")
    abline(h=(1-1/(N_pop))^10,col="yellow")
    abline(h=(1-1/(N_pop))^50,col="blue")
    text(addtxt$l,addtxt$h,addtxt$txt[s+(b-1)*length(sigmaX_set)],srt=addtxt$srt,font=addtxt$font,col=addtxt$col,cex=1.25)
  }
}
dev.print(device = jpeg, file = "Evo_Outcomes_Ne100.jpeg", width = 780*3,height=960*3,res=300,type="cairo")
