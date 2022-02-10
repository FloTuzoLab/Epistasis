library(plot3D)
require(Rmpfr)
library(RColorBrewer)

options(digits=20)#To avoid problems with rounding numbers
#Defining features of plots
ncol=128
ncontour=5
#palet<-colorRampPalette(c("white","black"))(ncol)
jet.colors <- colorRampPalette(c("steelblue1", "yellow", "tomato")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)

Fit_mat<-matrix(nrow=250,ncol=250)
FitList<-list()
Ruggedness_set<-c(1,10,50)
N_reso=250
for(rug in 1:length(Ruggedness_set)){
  Ruggedness<-Ruggedness_set[rug]
  for (k in 1:200){
    optloc<-runif(2,0,2.5)
    print(optloc[1])
    for (i in 1:250){
      for (j in 1:250){
        if(i==1 && j==1){
          Fit_mat[i,j]<-0
        }
        else if(i==1){
          Fit_mat[i,j]<-0
        }
        else if(j==1){
          Fit_mat[i,j]<-0
        }
        else{
          if(k>1){
            Fit_mat[i,j]<-exp(-(i/100-optloc[1])^2*Ruggedness-(j/100-optloc[2])^2*Ruggedness)+Fit_mat[i,j]#exp(-(i/100-optloc[1])^2-(j/100-optloc[2])^2)
          }
          else{
            Fit_mat[i,j]<-exp(-(i/100-optloc[1])^2*Ruggedness-(j/100-optloc[2])^2*Ruggedness)
          }
          #print(Fit_mat[i,j])
        }
      }
    }
  }
  FitList[[rug]]<-Fit_mat
}

image.plot(Fit_mat/max(Fit_mat))

Fitfacet <- (Fit_mat[-1, -1] + Fit_mat[-1, -N_reso] + Fit_mat[-N_reso, -1] + Fit_mat[-N_reso, -N_reso])/4
ncol=256
Fit_cut<-cut(Fitfacet,ncol)
par(mfrow=c(1,3),pin=c(6,5),mai=c(0.5,0,0.5,0.5))
persp3D(c(1:250),c(1:250),xlab="Genotype",ylab="Genotype",zlab="Fitness",zlim=c(0,1.05),FitList[[1]]/max(FitList[[1]]),theta = -20, phi = 20, expand = 1
        ,clab = c("Fitness scale"),shade=0.05,border="black",axes = T, ticktype="simple",ticks = FALSE,#palet[P_cut]
        box=TRUE,cex.lab=1.25,cex.axis=1,r=5,lwd=0.1,lty=1,legend.shrink=0.4,colvar=FitList[[1]]/max(FitList[[1]]),col=palet,lwd=1,labcex=1.5,colkey=FALSE)
text(0.06,0.16,"A",srt=0,font=2,col=1,cex=2,family="serif")
persp3D(c(1:250),c(1:250),xlab="Genotype",ylab="Genotype",zlab="Fitness",zlim=c(0,1.05),FitList[[2]]/max(FitList[[2]]),theta = -20, phi = 20, expand = 1
        ,clab = c("Fitness scale"),shade=0.05,border="black",axes = T, ticktype="simple",ticks = FALSE,#palet[P_cut]
        box=TRUE,cex.lab=1.25,cex.axis=1,r=5,lwd=0.1,lty=1,legend.shrink=0.4,colvar=FitList[[2]]/max(FitList[[2]]),col=palet,lwd=1,labcex=1,colkey=FALSE)
text(0.06,0.16,"B",srt=0,font=2,col=1,cex=2,family="serif")
persp3D(c(1:250),c(1:250),xlab="Genotype",ylab="Genotype",zlab="Fitness",zlim=c(0,1.05),FitList[[3]]/max(FitList[[3]]),theta = -20, phi = 20, expand = 1
        ,clab = c("Fitness scale"),shade=0.05,border="black",axes = T, ticktype="simple",ticks = FALSE,#palet[P_cut]
        box=TRUE,cex.clab=3,cex.lab=1.25,cex.axis=1,r=5,lwd=0.1,lty=1,legend.shrink=1,colvar=FitList[[3]]/max(FitList[[3]]),col=palet,lwd=1,labcex=1.5)
text(0.06,0.16,"C",srt=0,font=2,col=1,cex=2,family="serif")
#colkey(clim=c(0,1),col=c("white"),side=4,length=0.8,width=0.7,cex.axis=0.75,side.clab=2,cex.clab=0.8,dist=-0.105,add=T,tck=2,hadj=2.5)
setwd(dir="~")
setwd(dir="Desktop/")
dev.print(device = jpeg, file = "RuggedLandscapes.jpeg", width = 1350*3,height=460*3,res=300,type="cairo")

Fit_TO_3D<-matrix(nrow=100,ncol=100)
x<-seq(-1,1,length=100)
y<-seq(-1,1,length=100)
for (k in 1:100){
  for (j in 1:100){
    Fit_TO_3D[k,j]<-x[k]^3-3*x[k]*y[j]^2-y[j]*x[k]^2
  }
}
Fitfacet <- (Fit_TO_3D[-1, -1] + Fit_TO_3D[-1, -N_reso] + Fit_TO_3D[-N_reso, -1] + Fit_TO_3D[-N_reso, -N_reso])/4
persp3D(c(1:100),c(1:100),xlab="Genotype",ylab="Genotype",zlab="Fitness",zlim=c(-1,1),Fit_TO_3D/max(Fit_TO_3D),theta = 30, phi = 10, expand =1
        ,clab = c("Fitness scale"),shade=0.05,border="black",axes = T, ticktype="simple",ticks = FALSE,#palet[P_cut]
        box=FALSE,cex.lab=1.25,cex.axis=1,r=5,lwd=0.1,lty=1,legend.shrink=0.4,colvar=Fit_TO_3D/max(Fit_TO_3D),col=palet,lwd=1,labcex=1.5,colkey=FALSE)

