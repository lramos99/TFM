library(robust)
library(MASS)
library(robustbase)
library(modi)
library(xlsx)
library(MASS)
library(corrplot)
library(rcompanion)
library(featureCorMatrix)
library(pracma)
library(rcompanion)


classify=function(pd){
  n=nrow(pd)
  pdq <- data.frame()[1:n, ]#quantitative
  pdb <- data.frame()[1:n, ]#binary
  pdc <- data.frame()[1:n, ]#categorical
  i=2
  while (i<=ncol(pd)){
    if (is.factor(pd[,i])){
      if (length(levels(pd[,i]))==2){pdb=cbind(pdb, pd[,i])}
      else {pdc=cbind(pdc, pd[,i])}
    }
    else {pdq=cbind(pdq, pd[,i])}
    i=i+1
  }
  return(list(pdq, pdb,pdc))}

varg=function(mat){
  n=ncol(mat)
  v=sum(mat)/(2*n^2)
  mat1=mat/v
  return(mat1)
}

mahal1=function(pd){
  cov1=covMcd(pd)$cov
  n=nrow(pd)
  D2=matrix(ncol=n, nrow = n)
  pd1=pd
  for (i in 1: ncol(pd)){
    pd1[,i][is.na(pd[, i])]<-median(pd1[, i], na.rm = T)
  }
  for (i in 1:n){
    D2[i,] <- mahalanobis(pd1, center=as.numeric(pd1[i,]), cov=cov1, pairwise=TRUE)
  }
  return(D2)}

jacc=function(pd){
  n=nrow(pd)
  p=ncol(pd)-1
  s=matrix(ncol=n, nrow=n)
  for (i in 1:n){
    j=1
    while (j<=i){
      a = sum((pd[i,]==1 & pd[j,]==1), na.rm=T)
      d = sum((pd[i,]==0 & pd[j,]==0), na.rm=T)
      if (a==0 & i!= j) s[i,j]=0
      if (i==j) s[i,j]=1
      if (a!= 0) s[i,j]=a/(p-d)
      s[j,i]=s[i,j]
      j=j+1
    }}
  D2=1-s
  return(D2)
}

coin=function(pd){
  p=ncol(pd)
  n=nrow(pd)
  s=matrix(ncol=n, nrow=n)
  for (i in 1:n){
    j=1
    while (j<=i){
      alf = sum(pd[i,]==pd[j,], na.rm=T)
      s[i,j]=alf/p
      s[j,i]=s[i,j]
      j=j+1
    }}
  D2=1-s
  return(D2)
}

####Algorithm time series#####

#step 1
matvar<-function(pd, k){
  X<-matrix(ncol=17,nrow=length(unique(pd[,1])))
  years<-2002:2018
  for (j in 1:17){
    u<-pd[which(pd$year==years[j]),k]
    X[,j]<-as.numeric(u)}
  return(X)}

#step 2
matmean<-function(pd, k){
  X<-matvar(pd,k)
  D<-colMeans(X)
  return(D)}

#step3
mahalstep3=function(pd, k){
  cov1=covMcd(matvar(pd,k))$cov
  n=length(unique(pd[,1]))
  D2=list()
  for (i in 1:n){
    D2[i] <- as.numeric((t(matvar(pd,k)[i,]-as.numeric(matmean(pd,k)))%*%solve(cov1))%*%(matvar(pd,k)[i,]-as.numeric(matmean(pd,k))))
  }
  return(unlist(D2))}

#step 4
signfunct1<-function(pd,k){
  X=mahalstep3(pd,k)
  xk=rep(0,length(X))
  mat=matvar(pd,k)
  n=length(unique(pd[,1]))
  xmean=matmean(pd,k)
  for (i in 1:n){
    s=sum(mat[i,]-xmean)
    if (s>0){sign=1}
    else{sign=-1}
    xk[i]=sign*X[i]
  }
  return(xk)}

signfunct2<-function(pd){
  Xt=matrix(nrow = length(unique(pd[,1])), ncol=13)
  for (k in 3:15){Xt[,k-2]=signfunct1(pd,k)}
  return(Xt[,c(1,3:13)])
}

#5 and part of 6
Dmr<-function(pd){
  Xt<-signfunct2(pd)
  Dmr<-mahal1(Xt)
  Dmr=varg(Dmr)
  return(Dmr)
}

#part of 6 and 7
robmds<-function(pd){
  n<-nrow(pd)
  pdq=as.data.frame(classify(pd)[1])
  pdb=as.data.frame(classify(pd)[2])
  pdc=as.data.frame(classify(pd)[3])
  
  #three matrices
  distq=mahal1(pdq)
  D2q=varg(varg(distq)+Dmr(temp))
  distb=jacc(pdb)
  D2b=varg(distb)
  distc=coin(pdc)
  D2c=varg(distc)
  D2=varg(D2q+D2c+D2b)
  one=rep(1, n)
  id=diag(n)
  H=id-1/n*one%*%t(one)
  B=-1/2*(H%*%D2%*%H)
  if (min(eigen(B)$values)<0){
    c=2*abs(min(eigen(B)$values))
    D2=D2+c
    D2=D2-c*diag(n)
    B=-1/2*H%*%D2%*%H}
  eigenv=round(eigen(B)$values, 5)#round to the fifth decimal
  A=diag(eigenv,n,n)
  U=eigen(B)$vectors
  X=U%*%(A^(0.5))
  coord=X[,c(1,2)]
  #par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  continents<-pd[,7]
  plot(coord[,1],coord[,2], xlab="", ylab="", pch=20, col=continents, ylim=c(-2,0.8), xlim=c(-0.6,2))
  text(coord[c(21,23,25,42,48),2] ~coord[c(21,23,25,42,48),1], col=c(2,1,1,2,3),offset = 0.6, labels=pd[c(21,23,25,42,48),1], cex= 0.7, pos=3)
  return(coord)}


#####data####
#setwd("C:/Users/lucia/Downloads/TFM/data")
temp<-read.xlsx("temp1.xlsx",1, header=TRUE)
notemp<-read.xlsx("notemp1.xlsx",1, header=TRUE)
temp<-temp[,-c(1,5)]
notemp<-notemp[-29,-1]
temp<-temp[which(temp[,1]!="Montenegro"),]
notemp$income<-as.factor(notemp$income)
notemp$law<-as.factor(notemp$law)
notemp$drinking_awareness=as.factor(notemp$drinking_awareness)
notemp$UHC_legislation=as.factor(notemp$UHC_legislation)
notemp$violence_treatment=as.factor(notemp$violence_treatment)
notemp$policy_plan<-as.factor(notemp$policy_plan)
notemp$continent<-as.factor(notemp$continent)
notemp$country<-as.factor(notemp$country)
temp$country<-as.factor(temp$country)

###analysis axes####

coord=robmds(notemp)
x=coord[,1]
y=coord[,2]

#time series
Xd=signfunct2(temp)
Xdax=cbind(x,cbind(y,Xd))
names<-names(temp[,c(3,5:15)])
names<-append(list("x", "y"), names)
Xdax=as.data.frame(Xdax)
colnames(Xdax)<-names
Xdaxc=cor(Xdax)
corrplot(Xdaxc,  method = 'color', order = 'alphabet', tl.cex = 0.7, type = "upper", tl.col = "darkblue", tl.srt=65)

#ordinal
notemp$income=as.numeric(notemp$income)
notemp$violence_treatment=as.numeric(notemp$violence_treatment)
data_cor <- cor(na.omit(cbind(y,cbind(x,notemp[ , c(2,14)]))),method = "spearman")
corydata <-sort(data_cor[1,-2], decreasing = T)
corxdata <-sort(data_cor[2,-1], decreasing = T)
barplot(corxdata, names.arg =names(corxdata), col="cornflowerblue",
        xaxt="n")
text(cex=0.9, x=linspace(0.7,3.4,3), y=-0.05, names(corxdata), xpd=TRUE, srt=45, pos=2)


barplot(corydata, names.arg =names(corydata), col="cornflowerblue",
        xaxt="n")
text(cex=0.9, x=linspace(0.7,3.4,3), y=-0.15, names(corydata)[1:3], xpd=TRUE, srt=45, pos=2)

#nominal and binary
xquant=quantile(x, prob=c(0,0.25,0.5,0.75,1))
yquant=quantile(y, prob=c(0,0.25,0.5,0.75,1))
categorizar<-function(x, xquant){
  xcat<-ifelse(x<xquant[2],0,ifelse(x<xquant[3],1,ifelse(x<xquant[4],2,3)))
  return(as.factor(xcat))}

xcat<-categorizar(x,xquant)
ycat<-categorizar(y, yquant)


lx=list()
lx[6]=cramerV(xcat, notemp[,4])
lx[2]=cramerV(xcat, notemp[,5])
lx[3]=cramerV(xcat, notemp[,7])
lx[4]=cramerV(xcat, notemp[,12])
lx[5]=cramerV(xcat, notemp[,13])
lx[1]=cramerV(xcat, xcat)
lx=unlist(lx)

ly=list()
ly[6]=cramerV(ycat, notemp[,4])
ly[2]=cramerV(ycat, notemp[,5])
ly[3]=cramerV(ycat, notemp[,7])
ly[4]=cramerV(ycat, notemp[,12])
ly[5]=cramerV(ycat, notemp[,13])
ly[1]=cramerV(ycat, ycat)
ly=unlist(ly)

barplot(lx, names.arg =append("x",names(notemp)[c(5,7,12,13,4)]), col="cornflowerblue", xaxt="n")
text(cex=0.9, x=linspace(0.9,7,6), y=-0.05, append("x",names(notemp)[c(5,7,12,13,4)]), xpd=TRUE, srt=45, pos=2)


barplot(ly, names.arg =append("y",names(notemp)[c(4,5,7,12,13)]), col="cornflowerblue",
        xaxt="n")
text(cex=0.9, x=linspace(0.9,7,6), y=-0.05, append("y",names(notemp)[c(5,7,12,13,4)]), xpd=TRUE, srt=45, pos=2)
