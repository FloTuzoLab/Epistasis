library(extraDistr)
setwd(dir="~")
setwd(dir="Desktop/Complementary-epistasis/Figures")

##1.Plotting Kumaraswamy distribution, which seems very appropriate for fitness effects
a_set<-c(0.5,1,10)
b_set<-c(1/10,1/5,1)
par(mfrow=c(1,3),mai=c(0.75,0.7,0.6,0.25))
col_set=c(2,1,3)
count=0
addtxt<-list(l=-0.03,h=5,txt=c("A","B","C","D","E","F","G","H","I"),srt = 0,font=2,col="black")
for (a in a_set){
  curve(dkumar(x,a,b_set[1]*a),0,1,col = 2,xlim=c(-0.05,1),ylim=c(0,5),ylab="k(x;a,b)",lwd=2,cex.lab=2,cex.axis=1.5)
  count=count+1
  text_leg=c()
  text(addtxt$l,addtxt$h,addtxt$txt[count],srt=addtxt$srt,font=addtxt$font,col=addtxt$col,cex=1.5)
  for(b in 1:length(b_set)){
    print(paste(a,b))
    curve(dkumar(x,a, b_set[b]*a),0,1 ,col = col_set[b],add=TRUE,lwd=2)
    text_leg=c(text_leg,paste("a=",a,", b=",b_set[b]*a," "))
  }
  legend(x=0.225,y=5,legend=text_leg,lwd=2,col=c(2,1,3),cex=1)
}
dev.print(device = jpeg, file = "Kumar_dist.jpeg", width = 1300*3,height=460*3,res=300,type="cairo")

n=10^6
Ne=10^5
b=2
(4*Ne/(4*Ne+n*b))^2
1-2*n*b/(4*Ne+n*b)
1-n/(2*Ne+n)
sqrt(1-n*b/((Ne+n*b)))

integrate(f, lower, upper, …, subdivisions = 100L,
          rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
          stop.on.error = TRUE, keep.xy = FALSE, aux = NULL)

integrand <- function(x) {1/((x+1)*sqrt(x))}
integrate(integrand, lower = 0, upper = Inf)

integrate(dbeta(shape1=0.1,shape2=0.2), 0.01, 0.05)


n=10^2
Ne=10^3
v=2*Ne-2
a=10
b=2
alpha_K=v/a
beta_K=n*b
beta_func<-function(x){
  x^(alpha_K)*(1-x)^(beta_K-1)
}

beta_M1_func<-function(x){
  x^(alpha_K+1/a)*(1-x)^(beta_K-1)
}
integrate(beta_M1_func,0,1)$value/integrateR(beta_func,0,1)$value
factorial(2*Ne)/(factorial(2*(Ne+n*b)))*2^(2*n*b)*factorial(Ne+n*b-1)/factorial(Ne-1)*factorial(Ne+n*b)/factorial(Ne)
sqrt(Ne/(Ne+n*b))
4*(Ne-2)/(4*Ne+n*b-2)
max(rkumar(10^8,2,10))
sqrt(Ne/(Ne+n*b))/max(rkumar(10^6,2,5))
1-1/2*n*b/(Ne+n*b)
(2*Ne-1)/(n+2*Ne-1)

prod=1
for(i in Ne:(Ne+n*b-1)){
  prod=prod*(i)*(i+1)*2^2/((2*i+1)*(2*i+2))
}
prod

(Ne/(Ne+n*b/4))^2
(Ne/(Ne+n*b))^(1/2)
(2*Ne/(2*Ne+n*b*a))^(1/10)
(2*Ne-1)/(2*Ne-1+n*b)
((2*(Ne+1))/(2*(Ne+1+2*n*b)))^(1/4)
((2*Ne-7/4)/(2*Ne-7/4+n*b/4))^(4)
(1-n*b/(v+n*b*a))

a=5
b=2
n=c(10^-2)
v=10^1
Kumar_func<-function(f){
  f^(v+a-1)*a*b*n*(1-f^a)^(b*n-1)
}

Kumar_func_M1<-function(f){
  f^(v+a)*a*b*n*(1-f^a)^(b*n-1)
}
norm<-integrate(Kumar_func,0,1)$value
M1<-integrate(Kumar_func_M1,0,1)$value
M1/norm

n=c(10^-0)
v=10^2
Kumar_func<-function(f){
  f^(v+a-1)*a*b*n*(1-f^a)^(b*n-1)
}

Kumar_func_M1<-function(f){
  f^(v+a)*a*b*n*(1-f^a)^(b*n-1)
}
norm<-integrate(Kumar_func,0,1)$value
M1<-integrate(Kumar_func_M1,0,1)$value
M1/norm
