rm(list = ls(all = TRUE))

library(RCurl)
library(INLA)
library(MCMCpack)

#define your working directory, where the data files are stored
workdir<-""

#prepare data
simx <- read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/Simulated%20Logistic%20Data%20With%20Multiple%20Modes%20%28Example%203%29/sim3-X.txt"),sep = ",")
simy <-  read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/Simulated%20Logistic%20Data%20With%20Multiple%20Modes%20%28Example%203%29/sim3-Y.txt"),sep = ",")
data.example <- cbind(simy,simx)
names(data.example)[1]="Y1"

data.example$V2<-(data.example$V10+data.example$V14)*data.example$V9
data.example$V5<-(data.example$V11+data.example$V15)*data.example$V12
#fparam <- c("Const",colnames(data)[-1])
fparam.example <- colnames(data.example)[-1]
fobserved.example <- colnames(data.example)[1]


formula1 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-c(1,3,5,7,9,11,13,15,17,19)],collapse = "+")))
formula2 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-c(1,3,7,9,13,15,17,19)],collapse = "+")))

M=100
res<-array(0,dim = c(100,5))
meanp<-0     
for(i in 1:M)
{
  print(i)
  
  precp<-runif(1,0,1)
  res[i,5]<-precp
  #Model 1 MCMC
  m1C <- MCMCprobit(formula1,burnin=1000,b0=meanp,data=data.example,B0=precp,marginal.likelihood="Chib95",mcmc=100000,verbose=0)
  mc2 <- MCMCprobit(formula2,burnin=1000,b0=meanp,data=data.example,B0=precp,marginal.likelihood="Chib95",mcmc=100000,verbose=0)
  
  BF <- BayesFactor(m1C,mc2)
  #print(BF$BF.logmarglike)
  res[i,1]<-BF$BF.logmarglike[1]
  res[i,2]<-BF$BF.logmarglike[2]
  #Model 1 INLA
  
  n <- 2000
  M1 <- inla(formula = formula1,data=data.example,family="binomial",Ntrials=rep(1,n),
             control.family=list(link="probit"),control.predictor=list(compute=T),
             control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
  M2 <- inla(formula = formula2,data=data.example,family="binomial",Ntrials=rep(1,n),
             control.family=list(link="probit"),control.predictor=list(compute=T),
             control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
  #M1$mlik
  #M2$mlik
  
  res[i,3]<-M1$mlik[1]
  res[i,4]<-M2$mlik[1]
}

write.csv(x = (res),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/comparison of mlik methods/prob1.csv")

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
  m1C <- MCMCprobit(formula1,burnin=1000,b0=meanp,data=data.example,B0=precp,marginal.likelihood="Chib95",mcmc=100000,verbose=0)
  mc2 <- MCMCprobit(formula2,burnin=1000,b0=meanp,data=data.example,B0=precp,marginal.likelihood="Chib95",mcmc=100000,verbose=0)
  
  BF <- BayesFactor(m1C,mc2)
  #print(BF$BF.logmarglike)
  res[i,1]<-BF$BF.logmarglike[1]
  res[i,2]<-BF$BF.logmarglike[2]
  #Model 1 INLA
  
  n <- 2000
  M1 <- inla(formula = formula1,data=data.example,family="binomial",Ntrials=rep(1,n),
             control.family=list(link="probit"),control.predictor=list(compute=T),
             control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
  M2 <- inla(formula = formula2,data=data.example,family="binomial",Ntrials=rep(1,n),
             control.family=list(link="probit"),control.predictor=list(compute=T),
             control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
  #M1$mlik
  #M2$mlik
  
  res[i,3]<-M1$mlik[1]
  res[i,4]<-M2$mlik[1]
}

write.csv(x = (res),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/comparison of mlik methods/prob2.csv")

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

plot(x = res[,1],y=res[,3],xlab="Chib's", ylab="INLA", main="Compare simulated ML with perturbed priors M1")
plot(x = res[,2],y=res[,4],xlab="Chib's", ylab="INLA", main="Compare simulated ML with perturbed priors M2")

plot(y=res[,1]-res[,3],x = res[,5],ylab = "difference",xlab="precision",main="Deifference between Chib's and INLA ML vs prior precision M1")
plot(y=res[,2]-res[,4],x = res[,5],ylab = "difference",xlab="precision",main="Deifference between Chib's and INLA ML vs prior precision M2")

plot(density(res[,1]-res[,3]), main="Compare Kernel Density of diff of ML with perturbed priors M1")
polygon(density(res[,1]-res[,3]), col="red", border="blue") 

plot(density(res[,2]-res[,4]), main="Compare Kernel Density of diff of ML with perturbed priors M2")
polygon(density(res[,2]-res[,4]), col="red", border="blue") 

# now for a fixed prior precision let us discover the concergence
meanp<-1
precp<-0.2

res<-array(0,dim = c(10,10,2))

for(i in 1:10)
{

  for(j in 1:10)
  {  
    print(paste(i," and ",j))
    m1C <- MCMCregress(formula1,burnin=100,b0=meanp,data=data.example,B0=precp,c0 = 2, d0 = 2,marginal.likelihood="Chib95",mcmc=100*2^i,verbose=0)
    BF <- BayesFactor(m1C)
    res[i,j,1]<-BF$BF.logmarglike
    res[i,j,2]<-1000*2^i
  }
}
M1 <- inla(formula = formula1,data=data.example,family="Gaussian",control.predictor=list(compute=T),control.family=list(hyper = list(prec = list(prior="loggamma",param=c(1, 1)))),
           control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)),
           control.inla = list(strategy = "laplace", npoints = 100, diff.logdens= 2.5, int.strategy = "grid",dz=1.93,interpolator="gaussian")
           )
M2 <- inla(formula = formula1,data=data.example,family="Gaussian",control.predictor=list(compute=T),control.family=list(hyper = list(prec = list(prior="loggamma",param=c(1, 1)))),
           control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp))
           )

plot(y = res[,,1],x = res[,,2], type = "p",col = 2,ylim = c(min(min(res[,,1]),M1$mlik[1]),max(max(res[,,1]),M1$mlik[1])),
     main="Compare ML between INLA and Chib's for different number of MCMC steps",
     xlab = "number of MCMC iterations",
     ylab = "mlik")
abline(a = M1$mlik[1], b = 0,col = 5)
abline(a = M2$mlik[1], b = 0,col = 4)

# now for a fixed prior precision let us discover the concergence
meanp<-50
precp<-50

res<-array(0,dim = c(10,10,2))

for(i in 1:10)
{
  
  for(j in 1:10)
  {  
    print(paste(i," and ",j))
    m1C <- MCMCregress(formula1,burnin=100,b0=meanp,data=data.example,B0=precp,c0 = 2, d0 = 2,marginal.likelihood="Chib95",mcmc=100*2^i,verbose=0)
    BF <- BayesFactor(m1C)
    res[i,j,1]<-BF$BF.logmarglike
    res[i,j,2]<-1000*2^i
  }
}

M1 <- inla(formula = formula1,data=data.example,family="Gaussian",control.predictor=list(compute=T),control.family=list(hyper = list(prec = list(prior="loggamma",param=c(1, 1)))),
           control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)),
           control.inla = list(strategy = "laplace", npoints = 100, diff.logdens= 2.5, int.strategy = "grid",dz=3.85,interpolator="gaussian")
)
M2 <- inla(formula = formula1,data=data.example,family="Gaussian",control.predictor=list(compute=T),control.family=list(hyper = list(prec = list(prior="loggamma",param=c(1, 1)))),
           control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp))
)

plot(y = res[,,1],x = res[,,2], type = "p",col = 2,ylim = c(min(min(res[,,1]),min(M1$mlik[1],M2$mlik[1])),max(max(res[,,1]),max(M1$mlik[1],M2$mlik[1]))),
     main="Compare ML between INLA and Chib's for different number of MCMC steps",
     xlab = "number of MCMC iterations",
     ylab = "mlik")
abline(a = M1$mlik[1], b = 0,col = 5)
abline(a = M2$mlik[1], b = 0,col = 4)


# now a harmonic mean estimator is addressed with respect to INLA
# see https://radfordneal.wordpress.com/2008/08/17/the-harmonic-mean-of-the-likelihood-worst-monte-carlo-method-ever/

harmonic.mean.marg.lik <- function (x, s0, s1, n)
{ 
  post.prec <- 1/s0^2 + 1/s1^2
  t <- rnorm(n,(x/s1^2)/post.prec,sqrt(1/post.prec))
  lik <- dnorm(x,t,s1)
  return(log(1/mean(1/lik)))
}

true.marg.lik <- function (x,s0,s1)
{ 
  log(dnorm(x,0,sqrt(s0^2+s1^2)))
}

# case 1
spe<-2
s.h<-NULL
for(i in 1:10)
{
  print(i)
  s.h<-c(s.h,(harmonic.mean.marg.lik(spe,10,1,10^7)))
}

s.t<-(true.marg.lik(spe,10,1))

data<-as.data.frame(list(y = spe,x=rep(1,length(spe))))



out = inla(y ~ x-1, data = data,control.compute=list(mlik = T),
           control.family = list(
             hyper = list(
               prec = list(
                 initial = log(1), fixed = T
               )
             )
           ),
           control.fixed=list(
             mean = 0,
             prec = 0.01
           )
)

summary(out)

out$mlik
s.t
s.h

# case 2

spe<-2
s.h<-NULL
for(i in 1:10)
{
  print(i)
  s.h<-c(s.h,(harmonic.mean.marg.lik(spe,0.1,1,10^7)))
}
s.t<-(true.marg.lik(spe,0.1,1))

data<-as.data.frame(list(y = spe,x=rep(1,length(spe))))



out = inla(y ~ x-1, data = data,control.compute=list(mlik = T),
           control.family = list(
             hyper = list(
               prec = list(
                 initial = log(1), fixed = T
               )
             )
           ),
           control.fixed=list(
             mean = 0,
             prec = 100
           )
)

summary(out)

out$mlik
s.t
s.h

#m1C <- MCMCregress(formula = y~x-1,burnin=1000,b0=0,data=data,B0=1,c0 = 2, d0 = 2,marginal.likelihood="Chib95",mcmc=1000000,verbose=0)
#BayesFactor(m1C)$BF.logmarglike


#out = inla(y ~ x-1, data = data,control.compute=list(mlik = T),
#           control.family=list(hyper = list(prec = list(prior="loggamma",param=c(1, 1)))),
#           control.fixed=list(
#             mean = 0,
#             prec = 1
#           )
#)
#out$mlik
# case 3

spe<-2
s.h<-NULL
for(i in 1:10)
{
  print(i)
  s.h<-c(s.h,(harmonic.mean.marg.lik(spe,1000,1,10^7)))
}

s.t<-(true.marg.lik(spe,1000,1))

data<-as.data.frame(list(y = spe,x=rep(1,length(spe))))



out = inla(y ~ x-1, data = data,control.compute=list(mlik = T),
           control.family = list(
             hyper = list(
               prec = list(
                 initial = log(1), fixed = T
               )
             )
           ),
           control.fixed=list(
             mean = 0,
             prec = (1/1000)^2
           )
)

summary(out)

out$mlik
s.t
s.h


