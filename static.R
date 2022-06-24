library(robust)
library(MASS)
library(robustbase)
library(modi)
library(xlsx)
#### functions implemented####
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
  return(list(pdq, pdb,pdc))}#splits dataframe by type of variable

varg=function(mat){
  n=ncol(mat)
  v=sum(mat)/(2*n^2)
  mat1=mat/v
  return(mat1)
}#puts geom. var. equal 1

gower=function(pd){
  n=nrow(pd)
  p=ncol(pd)-1
  pdq=as.data.frame(classify(pd)[1])
  pdb=as.data.frame(classify(pd)[2])
  pdc=as.data.frame(classify(pd)[3])
  s=matrix(ncol=n, nrow=n)
  for (i in 1:n){
    j=1
    while (j<=i){
      alf=sum(pdc[i,]==pdc[j,], na.rm=T)
      d=sum((pdb[i,]==0 & pdb[j,]==0), na.rm=T)
      a=sum((pdb[i,]==1 & pdb[j,]==1), na.rm=T)
      m=apply(apply(as.matrix(pdq),2,range, na.rm=T),2,diff, na.rm=T)
      s[i,j]=(sum(1-abs(as.numeric(pdq[i,])-as.numeric(pdq[j,]))/m, na.rm=T)+alf+a)/(p-d)
      s[j,i]=s[i,j]
      j=j+1
    }}
  D2=1-s
  D2=varg(D2)#geometric variability 1
  return(D2)}#computes gower

mahal=function(pd){
  cov1=covMcd(pd)$cov
  n=nrow(pd)
  D2=matrix(ncol=n, nrow = n)
  pd1=pd
  for (i in 1: ncol(pd)){
    pd1[,i][is.na(pd[, i])]<-mean(pd1[, i], na.rm = T)
  }
  for (i in 1:n){
    D2[i,] <- mahalanobis(pd1, center=as.numeric(pd1[i,]), cov=cov1, pairwise=TRUE)
    }
  return(D2)}#computes mahalannobis

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
}#jaccard similarities

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
}#coincidence similarity

robustdist=function(pd){
  pdq=as.data.frame(classify(pd)[1])
  pdb=as.data.frame(classify(pd)[2])
  pdc=as.data.frame(classify(pd)[3])
  
  #three matrices
  distq=mahal(pdq)
  D2q=varg(distq)
  distb=jacc(pdb)
  D2b=varg(distb)
  distc=coin(pdc)
  D2c=varg(distc)
  D2=D2q+D2c+D2b
  return(D2)}#robust distance combining the three functions above


mds<-function(pd, p){#If p=0 gower if p=1 robust
  n=nrow(pd)
  if (p==0){D2=gower(pd)}
  else {D2=robustdist(pd)}
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
  continents=pd$continent
  #we multiply by -1 to turn around the figure
  if (coord[25,1]<=0){coord[,1]=coord[,1]*(-1)}
  if (coord[25,2]<=0){coord[,2]=coord[,2]*(-1)}
  #par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  opar <- par(no.readonly = TRUE)
  
  # Cambiar los márgenes del gráfico (el cuarto es el margen derecho)
  par(mar =c(5.1, 4.1, 4.1, 7.1))
  plot(coord[,1],coord[,2], xlab="", ylab="", pch=20, col=continents)
  text(coord[c(1,23,25,26,43,49),2] ~coord[c(1,23,25,26,43,49),1], col=c(5,1,"black","red","red",3),offset = 0.6, labels=pd[c(1,23,25,26,43,49),1], cex= 0.7
       , pos=2)
}#mds algorithm for representation

####data######

#setwd("C:/Users/lucia/Downloads/TFM/data")
temp<-read.xlsx("temp1.xlsx",1, header=TRUE)
notemp<-read.xlsx("notemp1.xlsx",1, header=TRUE)
temp<-temp[,-1]
notemp<-notemp[,-1]
notemp$income<-as.factor(notemp$income)
notemp$law<-as.factor(notemp$law)
notemp$drinking_awareness=as.factor(notemp$drinking_awareness)
notemp$UHC_legislation=as.factor(notemp$UHC_legislation)
notemp$violence_treatment=as.factor(notemp$violence_treatment)
notemp$policy_plan<-as.factor(notemp$policy_plan)
notemp$continent<-as.factor(notemp$continent)

###MDS representation for static ####
mds(notemp,1)
mds(notemp,0)
