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

# now for a fixed prior precision let us discover the concergence
res<-array(0,dim = c(3,6,2))
meanp<-0
precp<-(0.001)^2
i=10
for(j in 1:5)
{  
  seed = i+i*j + j*runif(n = 1,min = 0,max = 1000)
  print(paste(i," and ",j))
  m1C <- MCMCprobit(formula1,seed = seed,burnin=1000,b0=meanp,data=data.example,B0=precp,marginal.likelihood="Chib95",mcmc=100*2^i,verbose=0)
  BF1 <- BayesFactor(m1C)
  m2C <- MCMCprobit(formula2,seed = seed,burnin=1000,b0=meanp,data=data.example,B0=precp,marginal.likelihood="Chib95",mcmc=100*2^i,verbose=0)
  BF2 <- BayesFactor(m2C)
  res[1,j,1]<-BF1$BF.logmarglike
  res[1,j,2]<-BF2$BF.logmarglike
}

n <- 2000
M1 <- inla(formula = formula1,data=data.example,family="binomial",Ntrials=rep(1,n),
           control.family=list(link="probit"),control.predictor=list(compute=T),
           control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
res[1,6,1]<-M1$mlik[1]
M2 <- inla(formula = formula2,data=data.example,family="binomial",Ntrials=rep(1,n),
           control.family=list(link="probit"),control.predictor=list(compute=T),
           control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
res[1,6,2]<-M2$mlik[1]

meanp<-0
precp<-(0.1)^2
i=10
for(j in 1:5)
{  
  seed = i+i*j + j*runif(n = 1,min = 0,max = 1000)
  print(paste(i," and ",j))
  m1C <- MCMCprobit(formula1,seed = seed,burnin=1000,b0=meanp,data=data.example,B0=precp,marginal.likelihood="Chib95",mcmc=100*2^i,verbose=0)
  BF1 <- BayesFactor(m1C)
  m2C <- MCMCprobit(formula2,seed = seed,burnin=1000,b0=meanp,data=data.example,B0=precp,marginal.likelihood="Chib95",mcmc=100*2^i,verbose=0)
  BF2 <- BayesFactor(m2C)
  res[2,j,1]<-BF1$BF.logmarglike
  res[2,j,2]<-BF2$BF.logmarglike
}

n <- 2000
M1 <- inla(formula = formula1,data=data.example,family="binomial",Ntrials=rep(1,n),
           control.family=list(link="probit"),control.predictor=list(compute=T),
           control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
res[2,6,1]<-M1$mlik[1]
M2 <- inla(formula = formula2,data=data.example,family="binomial",Ntrials=rep(1,n),
           control.family=list(link="probit"),control.predictor=list(compute=T),
           control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
res[2,6,2]<-M2$mlik[1]

meanp<-0
precp<-(10)^2
i=10
for(j in 1:5)
{  
  seed = i+i*j + j*runif(n = 1,min = 0,max = 1000)
  print(paste(i," and ",j))
  m1C <- MCMCprobit(formula1,seed = seed,burnin=1000,b0=meanp,data=data.example,B0=precp,marginal.likelihood="Chib95",mcmc=100*2^i,verbose=0)
  BF1 <- BayesFactor(m1C)
  m2C <- MCMCprobit(formula2,seed = seed,burnin=1000,b0=meanp,data=data.example,B0=precp,marginal.likelihood="Chib95",mcmc=100*2^i,verbose=0)
  BF2 <- BayesFactor(m2C)
  res[3,j,1]<-BF1$BF.logmarglike
  res[3,j,2]<-BF2$BF.logmarglike
}

n <- 2000
M1 <- inla(formula = formula1,data=data.example,family="binomial",Ntrials=rep(1,n),
           control.family=list(link="probit"),control.predictor=list(compute=T),
           control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
res[3,6,1]<-M1$mlik[1]
M2 <- inla(formula = formula2,data=data.example,family="binomial",Ntrials=rep(1,n),
           control.family=list(link="probit"),control.predictor=list(compute=T),
           control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
res[3,6,2]<-M2$mlik[1]


# time measuring
system.time({
  i=10
m1C <- MCMCprobit(formula1,seed = seed,burnin=1000,b0=meanp,data=data.example,B0=precp,marginal.likelihood="Chib95",mcmc=100*2^i,verbose=0)
BF1 <- BayesFactor(m1C)
})

# now for a fixed prior precision let us discover the concergence
meanp<-1
precp<-0.2

res<-array(0,dim = c(10,10,2))

for(i in 1:10)
{
  
  for(j in 1:10)
  {  
    seed = i+i*j + j*runif(n = 1,min = 0,max = 1000)
    print(paste(i," and ",j))
    m1C <- MCMCprobit(formula1,seed = seed,burnin=1000,b0=meanp,data=data.example,B0=precp,marginal.likelihood="Chib95",mcmc=100*2^i,verbose=0)
    BF <- BayesFactor(m1C)
    res[i,j,1]<-BF$BF.logmarglike
    res[i,j,2]<-100*2^i
  }
}
n=2000

M1 <- inla(formula = formula1,data=data.example,family="binomial",Ntrials=rep(1,n),
           control.family=list(link="probit"),control.predictor=list(compute=T),
           control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))

plot(y = res[,,1],x = res[,,2], type = "p",col = 2,ylim = c(min(min(res[,,1]),M1$mlik[1]),max(max(res[,,1]),M1$mlik[1])),
     xlab = "",
     ylab = "",yaxt="n",xaxt="n")
mtext("MLIK", side=2, line=2.7, cex=1.7)
mtext("Number of MCMC iterations", side=1, line=3.5, cex=1.7)
axis(2,cex.axis=1.5)
axis(1,cex.axis=1.5)
abline(a = M1$mlik[1], b = 0,col = 4)

