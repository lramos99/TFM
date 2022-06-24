library(shiny)
library(dplyr)
library(robust)
library(MASS)
library(robustbase)
library(modi)
library(xlsx)


byyear<-function(year){
  tempy<-temp[temp$year==year,-2]
  pd<-merge(tempy, notemp, by = "country")
  return(pd)
}

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
  return(D2)}

mahal<-function(pd){
  pd1=pd[,c(4:ncol(pd))]
  cov1=covMcd(pd1)$cov
  n1=nrow(pd1)
  D2=matrix(ncol=n1, nrow = n1)
  for (i in 1: ncol(pd1)){
    pd1[,i][is.na(pd1[, i])]<-median(pd1[, i], na.rm = T)
  }
  
  for (i in 1:n1){
    D2[i,] <- mahalanobis(pd1, center=as.numeric(pd1[i,]), cov=cov1, pairwise=TRUE)
  }
  D2<-varg(D2)
  return(D2)
}

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
  return(D2)}

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
  if (coord[21,1]>=coord[48,1]){coord[,1]=coord[,1]*(-1)}
  if (coord[42,2]<=0){coord[,2]=coord[,2]*(-1)}
  continents<-pd[,"continent"]
  plot(coord[,1],coord[,2], xlab="", ylab="", pch=20, col=continents, ylim=c(-2,2.3), xlim=c(-5.1,5.2))
  text(coord[c(21,23,25,42,48),2] ~coord[c(21,23,25,42,48),1], col=c(2,1,1,2,3),offset = 0.6, labels=pd[c(21,23,25,42,48),1], cex= 0.7, pos=3)
}


setwd("C:/Users/lucia/Downloads/TFM/data")
temp<-read.xlsx("temp1.xlsx",1, header=TRUE)
notemp<-read.xlsx("notemp1.xlsx",1, header=TRUE)
temp<-temp[,-1]
notemp<-notemp[,-1]
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


ui <- fluidPage(
    h2("Robust representation by year"),
    plotOutput("plot"),
    sliderInput("n", "Year", 2002, 2018, 2008)
)


server <- function(input, output, session) {
  output$plot <- renderPlot({
    mds(byyear(input$n),1)
  })
  }

shinyApp(ui = ui, server = server)
