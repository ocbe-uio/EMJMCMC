rm(list = ls(all = TRUE))

library(RCurl)
library(INLA)
library(MCMCpack)

#define your working directory, where the data files are stored
workdir<-""

#prepare data
simx <- read.table(text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/US%20Data/simcen-x1.txt"),sep = ",")
simy <- read.table(text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/US%20Data/simcen-y1.txt"))
data.example <- cbind(simy,simx)
names(data.example)[1]="Y"

formula1 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-c(1,3,5,7,9,11,13,15)],collapse = "+")))
formula2 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-c(1,3,7,9,13)],collapse = "+")))

M=100
res<-array(0,dim = c(100,5))
meanp<-0
for(i in 1:M)
{
  print(i)
  
  precp<-runif(1,0,1)
  res[i,5]<-precp
  #Model 1 MCMC
  m1C <- MCMCregress(formula1,burnin=1000,b0=meanp,data=data.example,B0=precp,c0 = 2, d0 = 2,marginal.likelihood="Chib95",mcmc=100000,verbose=0)
  mc2 <- MCMCregress(formula2,burnin=1000,b0=meanp,data=data.example,B0=precp,c0 = 2, d0 = 2,marginal.likelihood="Chib95",mcmc=100000,verbose=0)
  
  BF <- BayesFactor(m1C,mc2)
  #print(BF$BF.logmarglike)
  res[i,1]<-BF$BF.logmarglike[1]
  res[i,2]<-BF$BF.logmarglike[2]
  #Model 1 INLA
  
  M1 <- inla(formula = formula1,data=data.example,family="Gaussian",control.predictor=list(compute=T),control.family=list(hyper = list(prec = list(prior="loggamma",param=c(1, 1)))),
             control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
  
  M2 <- inla(formula = formula2,data=data.example,family="Gaussian",control.predictor=list(compute=T),control.family=list(hyper = list(prec = list(prior="loggamma",param=c(1, 1)))),
             control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
  
  #M1$mlik
  #M2$mlik
  
  res[i,3]<-M1$mlik[1]
  res[i,4]<-M2$mlik[1]
  
  print(res[i,])
}

write.csv(x = (res),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/comparison of mlik methods/prob3.csv")

plot(density(res[,1]), main="Compare Kernel Density of ML with perturbed priors M1")
lines(density(res[,3]))
polygon(density(res[,1]), col="white", border="blue") 
polygon(density(res[,3]), col="white", border="red") 

plot(density(res[,2]), main="Compare Kernel Density of ML with perturbed priors M2")
lines(density(res[,4]))
polygon(density(res[,2]), col="white", border="blue") 
polygon(density(res[,4]), col="white", border="red") 

plot(x = res[,1],y=res[,3],xlab="Chib's", ylab="INLA", main="Compare simulated ML with perturbed priors M1")
plot(x = res[,2],y=res[,4],xlab="Chib's", ylab="INLA", main="Compare simulated ML with perturbed priors M2")

plot(y=res[,1]-res[,3],x = res[,5],ylab = "difference",xlab="precision",main="Deifference between Chib's and INLA ML vs prior precision M1")
plot(y=res[,2]-res[,4],x = res[,5],ylab = "difference",xlab="precision",main="Deifference between Chib's and INLA ML vs prior precision M2")

plot(density(res[,1]-res[,3]), main="Compare Kernel Density of diff of ML with perturbed priors M1")
polygon(density(res[,1]-res[,3]), col="red", border="blue") 

plot(density(res[,2]-res[,4]), main="Compare Kernel Density of diff of ML with perturbed priors M2")
polygon(density(res[,2]-res[,4]), col="red", border="blue") 



M=100
res<-array(0,dim = c(100,5))
meanp<-0
for(i in 1:M)
{
  print(i)
  
  precp<-runif(1,0,10)
  res[i,5]<-precp
  #Model 1 MCMC
  m1C <- MCMCregress(formula1,burnin=1000,b0=meanp,data=data.example,B0=precp,c0 = 2, d0 = 2,marginal.likelihood="Chib95",mcmc=100000,verbose=0)
  mc2 <- MCMCregress(formula2,burnin=1000,b0=meanp,data=data.example,B0=precp,c0 = 2, d0 = 2,marginal.likelihood="Chib95",mcmc=100000,verbose=0)
  
  BF <- BayesFactor(m1C,mc2)
  #print(BF$BF.logmarglike)
  res[i,1]<-BF$BF.logmarglike[1]
  res[i,2]<-BF$BF.logmarglike[2]
  #Model 1 INLA
  
  M1 <- inla(formula = formula1,data=data.example,family="Gaussian",control.predictor=list(compute=T),control.family=list(hyper = list(prec = list(prior="loggamma",param=c(1, 1)))),
             control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
  
  M2 <- inla(formula = formula2,data=data.example,family="Gaussian",control.predictor=list(compute=T),control.family=list(hyper = list(prec = list(prior="loggamma",param=c(1, 1)))),
             control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
  
  #M1$mlik
  #M2$mlik
  
  res[i,3]<-M1$mlik[1]
  res[i,4]<-M2$mlik[1]
  
  print(res[i,])
}

write.csv(x = (res),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/comparison of mlik methods/prob4.csv")

plot(density(res[,1]), main="Compare Kernel Density of ML with perturbed priors M1")
lines(density(res[,3]))
polygon(density(res[,1]), col="white", border="blue") 
polygon(density(res[,3]), col="white", border="red") 

plot(density(res[,2]), main="Compare Kernel Density of ML with perturbed priors M2")
lines(density(res[,4]))
polygon(density(res[,2]), col="white", border="blue") 
polygon(density(res[,4]), col="white", border="red") 

plot(x = res[,1],y=res[,3],xlab="Chib's", ylab="INLA", main="Compare simulated ML with perturbed priors M1")
plot(x = res[,2],y=res[,4],xlab="Chib's", ylab="INLA", main="Compare simulated ML with perturbed priors M2")

plot(y=res[,1]-res[,3],x = res[,5],ylab = "difference",xlab="precision",main="Deifference between Chib's and INLA ML vs prior precision M1")
plot(y=res[,2]-res[,4],x = res[,5],ylab = "difference",xlab="precision",main="Deifference between Chib's and INLA ML vs prior precision M2")

plot(density(res[,1]-res[,3]), main="Compare Kernel Density of diff of ML with perturbed priors M1")
polygon(density(res[,1]-res[,3]), col="red", border="blue") 

plot(density(res[,2]-res[,4]), main="Compare Kernel Density of diff of ML with perturbed priors M2")
polygon(density(res[,2]-res[,4]), col="red", border="blue") 



M=100
res<-array(0,dim = c(100,5))
meanp<-0
for(i in 1:M)
{
  print(i)
  
  precp<-runif(1,0,100)
  res[i,5]<-precp
  #Model 1 MCMC
  m1C <- MCMCregress(formula1,burnin=1000,b0=meanp,data=data.example,B0=precp,c0 = 2, d0 = 2,marginal.likelihood="Chib95",mcmc=100000,verbose=0)
  mc2 <- MCMCregress(formula2,burnin=1000,b0=meanp,data=data.example,B0=precp,c0 = 2, d0 = 2,marginal.likelihood="Chib95",mcmc=100000,verbose=0)
  
  BF <- BayesFactor(m1C,mc2)
  #print(BF$BF.logmarglike)
  res[i,1]<-BF$BF.logmarglike[1]
  res[i,2]<-BF$BF.logmarglike[2]
  #Model 1 INLA
  
  M1 <- inla(formula = formula1,data=data.example,family="Gaussian",control.predictor=list(compute=T),control.family=list(hyper = list(prec = list(prior="loggamma",param=c(1, 1)))),
             control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
  
  M2 <- inla(formula = formula2,data=data.example,family="Gaussian",control.predictor=list(compute=T),control.family=list(hyper = list(prec = list(prior="loggamma",param=c(1, 1)))),
             control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
  
  #M1$mlik
  #M2$mlik
  
  res[i,3]<-M1$mlik[1]
  res[i,4]<-M2$mlik[1]
  
  print(res[i,])
}

write.csv(x = (res),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/comparison of mlik methods/prob5.csv")

plot(density(res[,1]), main="Compare Kernel Density of ML with perturbed priors M1")
lines(density(res[,3]))
polygon(density(res[,1]), col="white", border="blue") 
polygon(density(res[,3]), col="white", border="red") 

plot(density(res[,2]), main="Compare Kernel Density of ML with perturbed priors M2")
lines(density(res[,4]))
polygon(density(res[,2]), col="white", border="blue") 
polygon(density(res[,4]), col="white", border="red") 

plot(x = res[,1],y=res[,3],xlab="", ylab="", main="",yaxt="n",xaxt="n")
mtext("INLA", side=2, line=2.7, cex=1.7)
mtext("Chib's method", side=1, line=3.5, cex=1.7)
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)
plot(x = res[,2],y=res[,4],xlab="", ylab="", main="",yaxt="n",xaxt="n")
mtext("INLA", side=2, line=2.7, cex=1.7)
mtext("Chib's method", side=1, line=3.5, cex=1.7)
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)

plot(y=res[,1]-res[,3],x = res[,5],ylab = "difference",xlab="precision",main="Deifference between Chib's and INLA ML vs prior precision M1")
plot(y=res[,2]-res[,4],x = res[,5],ylab = "difference",xlab="precision",main="Deifference between Chib's and INLA ML vs prior precision M2")

plot(density(res[,1]-res[,3]), main="Compare Kernel Density of diff of ML with perturbed priors M1")
polygon(density(res[,1]-res[,3]), col="red", border="blue") 

plot(density(res[,2]-res[,4]), main="Compare Kernel Density of diff of ML with perturbed priors M2")
polygon(density(res[,2]-res[,4]), col="red", border="blue") 
