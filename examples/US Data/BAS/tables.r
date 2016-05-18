
rm(list=ls())

install.packages("BAS")

library(BAS, lib.loc = "/mn/anatu/ansatte-u3/aliaksah/Desktop/competotrs/") # May need to install locally using lib.loc

source("allmodels.r")
source("bindec.r")

as.matrix.which <- function (x, which.models = NULL) 
           {
    namesx = x$namesx
    listobj = x$which
    if (!is.null(which.models)) 
        listobj = listobj[which.models]
    p = length(namesx)
    mat = t(sapply(listobj, function(x, dimp) {
        xx = rep(0, dimp)
        xx[x + 1] = 1
        xx
    }, p))
    colnames(mat) = namesx;
    mat
  }

  
setwd("/mn/anatu/ansatte-u3/aliaksah/Desktop/competotrs/")  
  
# 1. Begin: Calculate BIAS & MSE for inclusion probabilities (Tables 1 and 2)

simx <-  read.table("/mn/anatu/ansatte-u3/aliaksah/Desktop/competotrs/simcen-x.txt")
simy <- read.table("/mn/anatu/ansatte-u3/aliaksah/Desktop/competotrs/simcen-y.txt")
multidata <- cbind(simx,simy);
names(multidata)[ncol(simx)+1] <- "y";

multi.all <- bas.lm(y~.,data=multidata,n.models=2305,prior="g-prior",modelprior= uniform(),update=500,
                     alpha=nrow(multidata),initprobs="eplogp")
true.incprob <- multi.all$probne0[-1];

all <- as.matrix.which(multi.all);
all.dec <- apply(all[,-1],1,bindec);

truesummarg <- sum(exp(multi.all$logmarg))

sum(exp(multi.all$logmarg))/truesummarg

nsim <- 100; p <- 15; burnin <- 500;
#
incprob.rs.rm <- matrix(NA,nsim,p);
incprob.rst.rm <- matrix(NA,nsim,p);
incprob.mc3.rm <- matrix(NA,nsim,p);
#
incprob.rs.mc <- matrix(NA,nsim,p);
incprob.rst.mc <- matrix(NA,nsim,p);
incprob.mc3.mc <- matrix(NA,nsim,p);
#
incprob.bas.ep <- matrix(NA,nsim,p);
incprob.bas.unif <- matrix(NA,nsim,p);
incprob.srs <- matrix(NA,nsim,p);

for (i in 1:nsim)
  {
     print(i);
    # RS
    name <- paste("sim-rs",i,"dat",sep=".");
    multi.rs <- read.table(name);
    samp <- multi.rs[,-16];
    incprob.rs.mc[i,] <-  apply(samp[-c(1:burnin),],2,mean);
    dup <- duplicated(samp)
    samp.unique <- as.matrix(samp[dup==FALSE,]);
    marg.unique <- exp(multi.rs[dup==FALSE,16]);
    prob.unique <- marg.unique/sum(marg.unique);
    incprob <- prob.unique%*%samp.unique;
    incprob.rs.rm[i,] <- incprob;
    # RS-Thin
    name <- paste("sim-rs-thin",i,"dat",sep=".");
    multi.rst <- read.table(name);
    samp <- multi.rst[,-16];
    incprob.rst.mc[i,] <-  apply(samp[-c(1:burnin),],2,mean); 
    dup <- duplicated(samp)
    samp.unique <- as.matrix(samp[dup==FALSE,]);
    marg.unique <- exp(multi.rst[dup==FALSE,16]);
    prob.unique <- marg.unique/sum(marg.unique);
    incprob <- prob.unique%*%samp.unique;
    incprob.rst.rm[i,] <- incprob;
    # MC^3
    name <- paste("sim-mc3",i,"dat",sep=".");
    sim.mc3 <- read.table(name);
    samp <- sim.mc3[,-16];
    incprob.mc3.mc[i,] <-  apply(samp[-c(1:burnin),],2,mean); 
    dup <- duplicated(samp)
    samp.unique <- as.matrix(samp[dup==FALSE,]);
    marg.unique <- exp(sim.mc3[dup==FALSE,16]);
    prob.unique <- marg.unique/sum(marg.unique);
    incprob <- prob.unique%*%samp.unique;
    incprob.mc3.rm[i,] <- incprob;
    # BAS-eplogp
    name <- paste("sim.bas.ep",i,"Rdata",sep=".");
    load(name);
    incprob.bas.ep[i,] <- sim.bas.ep$probne0[-1];
    # BAS-uniform
    name <- paste("sim.bas.unif",i,"Rdata",sep=".");
    load(name);
    incprob.bas.unif[i,] <- sim.bas.unif$probne0[-1];
    # SRSWOR
    name <- paste("sim.srs",i,"Rdata",sep=".");
    load(name);
    incprob.srs[i,] <- sim.srs$probne0[-1]; 
  }

# MSE
mse.mc3.mc <- apply((incprob.mc3.mc-matrix(true.incprob,nrow=100,ncol=15,byrow=TRUE))^2,2,mean);
mse.rs.mc <- apply((incprob.rs.mc-matrix(true.incprob,nrow=100,ncol=15,byrow=TRUE))^2,2,mean);
mse.rst.mc <- apply((incprob.rst.mc-matrix(true.incprob,nrow=100,ncol=15,byrow=TRUE))^2,2,mean);
#
mse.mc3.rm <- apply((incprob.mc3.rm-matrix(true.incprob,nrow=100,ncol=15,byrow=TRUE))^2,2,mean);
mse.rs.rm <- apply((incprob.rs.rm-matrix(true.incprob,nrow=100,ncol=15,byrow=TRUE))^2,2,mean);
mse.rst.rm <- apply((incprob.rst.rm-matrix(true.incprob,nrow=100,ncol=15,byrow=TRUE))^2,2,mean);
#
mse.srs <- apply((incprob.srs-matrix(true.incprob,nrow=100,ncol=15,byrow=TRUE))^2,2,mean);
mse.bas.ep <- apply((incprob.bas.ep-matrix(true.incprob,nrow=100,ncol=15,byrow=TRUE))^2,2,mean);
mse.bas.unif <- apply((incprob.bas.unif-matrix(true.incprob,nrow=100,ncol=15,byrow=TRUE))^2,2,mean);

# BIAS
bias.mc3.mc <- apply((incprob.mc3.mc-matrix(true.incprob,nrow=100,ncol=15,byrow=TRUE)),2,mean);
bias.rs.mc <- apply((incprob.rs.mc-matrix(true.incprob,nrow=100,ncol=15,byrow=TRUE)),2,mean);
bias.rst.mc <- apply((incprob.rst.mc-matrix(true.incprob,nrow=100,ncol=15,byrow=TRUE)),2,mean);
#
bias.mc3.rm <- apply((incprob.mc3.rm-matrix(true.incprob,nrow=100,ncol=15,byrow=TRUE)),2,mean);
bias.rs.rm <- apply((incprob.rs.rm-matrix(true.incprob,nrow=100,ncol=15,byrow=TRUE)),2,mean);
bias.rst.rm <- apply((incprob.rst.rm-matrix(true.incprob,nrow=100,ncol=15,byrow=TRUE)),2,mean);
#
bias.srs <- apply((incprob.srs-matrix(true.incprob,nrow=100,ncol=15,byrow=TRUE)),2,mean);
bias.bas.ep <- apply((incprob.bas.ep-matrix(true.incprob,nrow=100,ncol=15,byrow=TRUE)),2,mean);
bias.bas.unif <- apply((incprob.bas.unif-matrix(true.incprob,nrow=100,ncol=15,byrow=TRUE)),2,mean);

ord.incprob <- order(true.incprob)

#MSE-incprob-matrix
mse.mat <- cbind(mse.bas.ep,mse.bas.unif,mse.mc3.mc,mse.mc3.rm,mse.rs.mc,mse.rs.rm,mse.rst.mc,mse.rst.rm,mse.srs)
rmse.mat.ord <- sqrt(mse.mat[ord.incprob,]);
View(round(cbind(true.incprob[ord.incprob],rmse.mat.ord*10^2),2))

#BIAS-incprob-matrix
bias.mat <- cbind(bias.bas.ep,bias.bas.unif,bias.mc3.mc,bias.mc3.rm,bias.rs.mc,bias.rs.rm,bias.rst.mc,bias.rst.rm,bias.srs)
bias.mat.ord <- bias.mat[ord.incprob,]
View(round(cbind(true.incprob[ord.incprob],bias.mat.ord*10^2),2))

# 1. End: Calculate BIAS & MSE for inclusion probabilities (Tables 1 and 2)


# 2. Begin: Calculate BIAS & MSE for model probabilities (Tables 1 and 2)

# Note mprobs.bias.bas actually saves avg bias^2

# BIAS 
mprobs.bias.mc3.mc <- rep(0,2^p);
mprobs.bias.rs.mc <- rep(0,2^p);
mprobs.bias.rst.mc <- rep(0,2^p);
#
mprobs.bias.mc3.rm <- rep(0,2^p);
mprobs.bias.rs.rm <- rep(0,2^p);
mprobs.bias.rst.rm <- rep(0,2^p);
#
mprobs.bias.srs <- rep(0,2^p);
mprobs.bias.bas.ep <- rep(0,2^p);
mprobs.bias.bas.unif <- rep(0,2^p);

# MSE
mprobs.mse.mc3.mc <- rep(0,2^p);
mprobs.mse.rs.mc <- rep(0,2^p);
mprobs.mse.rst.mc <- rep(0,2^p);
#
mprobs.mse.mc3.rm <- rep(0,2^p);
mprobs.mse.rs.rm <- rep(0,2^p);
mprobs.mse.rst.rm <- rep(0,2^p);
#
mprobs.mse.srs <- rep(0,2^p);
mprobs.mse.bas.ep <- rep(0,2^p);
mprobs.mse.bas.unif <- rep(0,2^p);

masses <- array(0,9)
unique.mods <- array(0,9)
total.mods <- array(0,9)

for (i in 1:nsim)
{
  print(i);
  # MC^3-MC
  sim.m <- read.table(file=paste("sim-mc3",i,"dat",sep="."))
  models.dec <- apply(sim.m[,-16],1,bindec);
  unique.mods[3]<-length(unique(models.dec));
  total.mods[3]<-length(models.dec);
  mprobs.mc <- table(models.dec)/length(models.dec);
  models.mc <- as.numeric(names(mprobs.mc));
  ind.match <- match(models.mc,all.dec);
  mprobs.all <- rep(0,2^p);
  mprobs.all[ind.match] <- mprobs.mc;
  mprobs.bias.mc3.mc <- mprobs.bias.mc3.mc+ (mprobs.all-multi.all$postprobs);
  mprobs.mse.mc3.mc <- mprobs.mse.mc3.mc+(multi.all$postprobs-mprobs.all)^2;
  dup <- duplicated(sim.m[,-16]);
  samp <- sim.m[dup==FALSE,-16];
  marg <- exp(sim.m[,16]);
  marg.unique <- marg[dup==FALSE];
  masses[3]<-masses[3] + sum(marg.unique)/truesummarg;
  rm(sim.m); rm(models.dec); rm(models.mc); rm(mprobs.mc); rm(ind.match); rm(mprobs.all);
  # RS-MC
  sim.rs <- read.table(file=paste("sim-rs",i,"dat",sep="."))
  models.dec <- apply(sim.rs[,-16],1,bindec);
  unique.mods[5]<-length(unique(models.dec));
  total.mods[5]<-length(models.dec);
  mprobs.mc <- table(models.dec)/length(models.dec);
  models.mc <- as.numeric(names(mprobs.mc));
  ind.match <- match(models.mc,all.dec);
  mprobs.all <- rep(0,2^p);
  mprobs.all[ind.match] <- mprobs.mc;
  mprobs.bias.rs.mc <-  mprobs.bias.rs.mc+(mprobs.all-multi.all$postprobs);
  mprobs.mse.rs.mc <-  mprobs.mse.rs.mc+(multi.all$postprobs-mprobs.all)^2;
  dup <- duplicated(sim.rs[,-16]);
  samp <- sim.rs[dup==FALSE,-16];
  marg <- exp(sim.rs[,16]);
  marg.unique <- marg[dup==FALSE];
  masses[5]<-masses[5] + sum(marg.unique)/truesummarg;
  rm(sim.rs); rm(models.dec); rm(models.mc); rm(mprobs.mc); rm(ind.match); rm(mprobs.all);
  # RS-Thin-MC
  sim.rst <- read.table(file=paste("sim-rs-thin",i,"dat",sep="."))
  models.dec <- apply(sim.rst[,-16],1,bindec);
  unique.mods[7]<-length(unique(models.dec));
  total.mods[7]<-length(models.dec);
  length(unique(models.dec));
  mprobs.mc <- table(models.dec)/length(models.dec);
  models.mc <- as.numeric(names(mprobs.mc));
  ind.match <- match(models.mc,all.dec);
  mprobs.all <- rep(0,2^p);
  mprobs.all[ind.match] <- mprobs.mc;
  mprobs.bias.rst.mc <-  mprobs.bias.rst.mc +(mprobs.all-multi.all$postprobs);
  mprobs.mse.rst.mc <-  mprobs.mse.rst.mc+(multi.all$postprobs-mprobs.all)^2;
  dup <- duplicated(sim.rst[,-16]);
  samp <- sim.rst[dup==FALSE,-16];
  marg <- exp(sim.rst[,16]);
  marg.unique <- marg[dup==FALSE];
  masses[7]<-masses[7] + sum(marg.unique)/truesummarg;
  rm(sim.rst); rm(models.dec); rm(models.mc); rm(mprobs.mc); rm(ind.match); rm(mprobs.all);
  # MC^3-RM
  sim.m <- read.table(file=paste("sim-mc3",i,"dat",sep="."))
  dup <- duplicated(sim.m[,-16]);
  samp <- sim.m[dup==FALSE,-16];
  marg <- exp(sim.m[,16]);
  marg.unique <- marg[dup==FALSE];
  models.dec <- apply(samp,1,bindec);
  unique.mods[4]<-length(unique(models.dec));
  total.mods[4]<-length(models.dec);
  ind.match <- match(models.dec,all.dec);
  mprobs.all <- rep(0,2^p);
  mprobs.all[ind.match] <- exp(sim.m[dup==FALSE,16])/sum(exp(sim.m[dup==FALSE,16]));
  mprobs.bias.mc3.rm <- mprobs.bias.mc3.rm+(mprobs.all-multi.all$postprobs);
  mprobs.mse.mc3.rm <- mprobs.mse.mc3.rm+(multi.all$postprobs-mprobs.all)^2;
  masses[4]<-masses[4]+sum(marg.unique)/truesummarg;
  rm(sim.m); rm(samp); rm(dup); rm(models.dec); rm(ind.match); rm(mprobs.all); 
  # RS-RM
  sim.rs <- read.table(file=paste("sim-rs",i,"dat",sep="."))
  dup <- duplicated(sim.rs[,-16]);
  samp <- sim.rs[dup==FALSE,-16];
  marg <- exp(sim.rs[,16]);
  marg.unique <- marg[dup==FALSE];
  models.dec <- apply(samp,1,bindec);
  unique.mods[6]<-length(unique(models.dec));
  total.mods[6]<-length(models.dec);
  ind.match <- match(models.dec,all.dec);
  mprobs.all <- rep(0,2^p);
  mprobs.all[ind.match] <- exp(sim.rs[dup==FALSE,16])/sum(exp(sim.rs[dup==FALSE,16]));
  mprobs.bias.rs.rm <- mprobs.bias.rs.rm+(mprobs.all-multi.all$postprobs);
  mprobs.mse.rs.rm <- mprobs.mse.rs.rm+(multi.all$postprobs-mprobs.all)^2;
  masses[6]<-masses[6]+sum(marg.unique)/truesummarg;
  rm(sim.rs); rm(samp); rm(dup); rm(models.dec); rm(ind.match); rm(mprobs.all); 
  # RS-Thin-RM
  sim.rst <- read.table(file=paste("sim-rs-thin",i,"dat",sep="."))
  dup <- duplicated(sim.rst[,-16]);
  samp <- sim.rst[dup==FALSE,-16];
  marg <- exp(sim.rst[,16]);
  marg.unique <- marg[dup==FALSE];
  models.dec <- apply(samp,1,bindec);
  unique.mods[8]<-length(unique(models.dec));
  total.mods[8]<-length(models.dec);
  ind.match <- match(models.dec,all.dec);
  mprobs.all <- rep(0,2^p);
  mprobs.all[ind.match] <- exp(sim.rst[dup==FALSE,16])/sum(exp(sim.rst[dup==FALSE,16]));
  mprobs.bias.rst.rm <- mprobs.bias.rst.rm+(mprobs.all-multi.all$postprobs);
  mprobs.mse.rst.rm <- mprobs.mse.rst.rm+(multi.all$postprobs-mprobs.all)^2;
  masses[8]<-masses[8]+sum(marg.unique)/truesummarg;
  rm(sim.rst); rm(samp); rm(dup); rm(models.dec); rm(ind.match); rm(mprobs.all); 
  # BAS-eplogp
  load(file=paste("sim.bas.ep",i,"Rdata",sep="."))
  models <- as.matrix.which(sim.bas.ep);
  models.dec <- apply(models[,-1],1,bindec);
  ind.match <- match(models.dec,all.dec);
  mprobs.all <- rep(0,2^p);
  mprobs.all[ind.match] <- sim.bas.ep$postprobs;
  mprobs.bias.bas.ep <- mprobs.bias.bas.ep+(mprobs.all-multi.all$postprobs);
  masses[1]<- masses[1]+sum(exp(sim.bas.ep$logmarg))/truesummarg;
  unique.mods[1]<-length(unique(models.dec));
  total.mods[1]<-length(models.dec);
  mprobs.mse.bas.ep <- mprobs.mse.bas.ep+(multi.all$postprobs-mprobs.all)^2;
  # BAS-uniform
  load(file=paste("sim.bas.unif",i,"Rdata",sep="."))
  models <- as.matrix.which(sim.bas.unif);
  models.dec <- apply(models[,-1],1,bindec);
  ind.match <- match(models.dec,all.dec);
  mprobs.all <- rep(0,2^p);
  mprobs.all[ind.match] <- sim.bas.unif$postprobs;
  mprobs.bias.bas.unif <- mprobs.bias.bas.unif+(mprobs.all-multi.all$postprobs);
  mprobs.mse.bas.unif <- mprobs.mse.bas.unif+(multi.all$postprobs-mprobs.all)^2;
  masses[2]<-masses[2]+sum(exp(sim.bas.unif$logmarg))/truesummarg;
  unique.mods[2]<-length(unique(models.dec));
  total.mods[2]<-length(models.dec);
}
# SRSWOR
load(file=paste("sim.srs",i,"Rdata",sep="."))
models <- as.matrix.which(sim.srs);
models.dec <- apply(models[,-1],1,bindec);
ind.match <- match(models.dec,all.dec);
mprobs.all <- rep(0,2^p);
mprobs.all[ind.match] <- sim.srs$postprobs;
mprobs.bias.srs <- mprobs.bias.srs+(mprobs.all-multi.all$postprobs);
mprobs.mse.srs <- mprobs.mse.srs+(multi.all$postprobs-mprobs.all)^2;
masses[9]<-100*sum(exp(sim.srs$logmarg))/truesummarg;
unique.mods[9]<-length(unique(models.dec));
total.mods[9]<-length(models.dec);
}

# BIAS
mprobs.bias.mc3.mc <- mprobs.bias.mc3.mc/nsim;
mprobs.bias.rs.mc <- mprobs.bias.rs.mc/nsim;
mprobs.bias.rst.mc <- mprobs.bias.rst.mc/nsim;
#
mprobs.bias.mc3.rm <- mprobs.bias.mc3.rm/nsim;
mprobs.bias.rs.rm <- mprobs.bias.rs.rm/nsim;
mprobs.bias.rst.rm <- mprobs.bias.rst.rm/nsim;
#
mprobs.bias.srs <- mprobs.bias.srs/nsim;
mprobs.bias.bas.ep <- mprobs.bias.bas.ep/nsim;
mprobs.bias.bas.unif <- mprobs.bias.bas.unif/nsim;

#MSE
mprobs.mse.mc3.mc <- mprobs.mse.mc3.mc/nsim;
mprobs.mse.rs.mc <- mprobs.mse.rs.mc/nsim;
mprobs.mse.rst.mc <- mprobs.mse.rst.mc/nsim;
#
mprobs.mse.mc3.rm <- mprobs.mse.mc3.rm/nsim;
mprobs.mse.rs.rm <- mprobs.mse.rs.rm/nsim;
mprobs.mse.rst.rm <- mprobs.mse.rst.rm/nsim;
#
mprobs.mse.srs <- mprobs.mse.srs/nsim;
mprobs.mse.bas.ep <- mprobs.mse.bas.ep/nsim;
mprobs.mse.bas.unif <- mprobs.mse.bas.unif/nsim;


mprobs.bias.vec <- sqrt(c(mean(mprobs.bias.bas.ep^2),mean(mprobs.bias.bas.unif^2),
                          mean(mprobs.bias.mc3.mc^2), mean(mprobs.bias.mc3.rm^2),mean(mprobs.bias.rs.mc^2),
                          mean(mprobs.bias.rs.rm^2),mean(mprobs.bias.rst.mc^2),mean(mprobs.bias.rst.rm^2),mean(mprobs.bias.srs^2)))

mprobs.rmse.vec <- sqrt(c(mean(mprobs.mse.bas.ep),mean(mprobs.mse.bas.unif),mean(mprobs.mse.mc3.mc),mean(mprobs.mse.mc3.rm),
                          mean(mprobs.mse.rs.mc),mean(mprobs.mse.rs.rm),mean(mprobs.mse.rst.mc),mean(mprobs.mse.rst.rm),
                          mean(mprobs.mse.srs)))

View(t(round(cbind(mod.bias = mprobs.bias.vec*10^(5),mod.rmse = mprobs.rmse.vec*10^(5),mod.masses = masses/nsim,unique.mods,total.mods),2)))

# 2. End: Calculate BIAS & MSE for model probabilities (Tables 1 and 2)


# 3. Begin: Calculate BIAS & MSE for BMA estimate of Y (Tables 1 and 2)
simy <- as.matrix(simy);
simy <- as.vector(simy);
simx <- as.matrix(simx);

g <- 100;
xbeta.all <- matrix(NA,nrow=2^p,ncol=100)

for (i in 1:2^p)
  {
    if(i%%1000==0) {print(i)}
    if (sum(all[i,-1])==0)
      {
        ols <- coef(lm(simy~1));
        alphapost <- ols[1];
        xbeta.all[i,] <- rep(alphapost,100);
      } else 
       {
         ols <- coef(lm(simy~simx[,all[i,-1]==1]))
         alphapost <- ols[1];
         betagpost <- g*(ols[-1])/(g+1);
         xbeta.all[i,] <- rep(alphapost,100)+simx[,all[i,-1]==1]%*%as.matrix(betagpost);
       }
  }

## Compute True E(Y) under BMA
yhatbma.true <- fitted.bma(multi.all,type="BMA",no.top=2^p)
## Compute True E(Y) under BMA

### BIAS and MSE for BAS/SRS estimates ####
nsim <- 100; n <- 100;
yhat.bias.srs <- rep(0,n);
yhat.bias.bas.ep <- rep(0,n);
yhat.bias.bas.unif <- rep(0,n);

yhat.mse.srs <- rep(0,n);
yhat.mse.bas.ep <- rep(0,n);
yhat.mse.bas.unif <- rep(0,n);

  
for (i in 1:nsim)
  {
     print(i);
     # SRS
     load(file=paste("sim.srs",i,"Rdata",sep="."))
     fit <- fitted.bma(sim.srs,type="BMA",no.top=length(sim.srs$logmarg))
     yhat.bias.srs <- yhat.bias.srs + (fit-yhatbma.true);
     yhat.mse.srs <- yhat.mse.srs+ (fit-yhatbma.true)^2;
     # BAS-eplogp
     load(file=paste("sim.bas.ep",i,"Rdata",sep="."))
     fit <- fitted.bma(sim.bas.ep,type="BMA",no.top=length(sim.bas.ep$logmarg))
     yhat.bias.bas.ep <- yhat.bias.bas.ep + (fit-yhatbma.true);
     yhat.mse.bas.ep <- yhat.mse.bas.ep + (fit-yhatbma.true)^2;
     # BAS-uniform
     load(file=paste("sim.bas.unif",i,"Rdata",sep="."))
     fit <- fitted.bma(sim.bas.unif,type="BMA",no.top=length(sim.bas.unif$logmarg))
     yhat.bias.bas.unif <- yhat.bias.bas.unif + (fit-yhatbma.true);
     yhat.mse.bas.unif <- yhat.mse.bas.unif+ (fit-yhatbma.true)^2;
  }

yhat.bias.srs <- yhat.bias.srs/nsim;
yhat.bias.bas.ep <- yhat.bias.bas.ep/nsim;
yhat.bias.bas.unif <- yhat.bias.bas.unif/nsim;

yhat.mse.srs <- yhat.mse.srs/nsim;
yhat.mse.bas.ep <- yhat.mse.bas.ep/nsim;
yhat.mse.bas.unif <- yhat.mse.bas.unif/nsim;


### BIAS and MSE for BAS/SRS estimates ####
    
### BIAS and MSE for MCMC(MC^3,RS,RS-Thin) estimates ####
yhat.mse.mc3.rm <- rep(0,n);
yhat.mse.rs.rm <- rep(0,n);
yhat.mse.rst.rm <- rep(0,n);

yhat.bias.mc3.rm <- rep(0,n);
yhat.bias.rs.rm <- rep(0,n);
yhat.bias.rst.rm <- rep(0,n);

yhat.mse.mc3.mc <- rep(0,n);
yhat.mse.rs.mc <- rep(0,n);
yhat.mse.rst.mc <- rep(0,n);

yhat.bias.mc3.mc <- rep(0,n);
yhat.bias.rs.mc <- rep(0,n);
yhat.bias.rst.mc <- rep(0,n);


for (i in 1:nsim)
  {
    print(i);
    # MC^3-MC
    sim.m <- read.table(file=paste("sim-mc3",i,"dat",sep="."))
    sampmod <- as.matrix(sim.m[,-16]);
    sampmod.dec <- apply(sampmod,1,bindec);
    ind.match <- match(sampmod.dec,all.dec);
    yhat.mc <- apply(xbeta.all[ind.match,],2,mean);
    yhat.bias.mc3.mc <- yhat.bias.mc3.mc+(yhat.mc-yhatbma.true);
    yhat.mse.mc3.mc <- yhat.mse.mc3.mc+(yhat.mc-yhatbma.true)^2;
    # MC^3-RM
    dup <- duplicated(sampmod.dec);
    sampmod.dec.u <- sampmod.dec[dup==FALSE];
    postprobs.u <- exp(sim.m[dup==FALSE,16])/sum(exp(sim.m[dup==FALSE,16]));
    ind.match <- match(sampmod.dec.u,all.dec);
    yhat.renorm <- t(t(postprobs.u)%*%xbeta.all[ind.match,]);
    yhat.bias.mc3.rm <- yhat.bias.mc3.rm+(yhat.renorm-yhatbma.true);
    yhat.mse.mc3.rm <- yhat.mse.mc3.rm+(yhat.renorm-yhatbma.true)^2;
    rm(sim.m); rm(sampmod); rm(sampmod.dec); rm(ind.match);rm(yhat.mc);
    rm(dup); rm(sampmod.dec.u); rm(postprobs.u); rm(yhat.renorm);
    # RS-MC
    sim.rs <- read.table(file=paste("sim-rs",i,"dat",sep="."))
    sampmod <- as.matrix(sim.rs[,-16]);
    sampmod.dec <- apply(sampmod,1,bindec);
    ind.match <- match(sampmod.dec,all.dec);
    yhat.mc <- apply(xbeta.all[ind.match,],2,mean);
    yhat.bias.rs.mc <- yhat.bias.rs.mc+(yhat.mc-yhatbma.true);
    yhat.mse.rs.mc <- yhat.mse.rs.mc+(yhat.mc-yhatbma.true)^2;
    # RS-RM
    dup <- duplicated(sampmod.dec);
    sampmod.dec.u <- sampmod.dec[dup==FALSE];
    postprobs.u <- exp(sim.rs[dup==FALSE,16])/sum(exp(sim.rs[dup==FALSE,16]));
    ind.match <- match(sampmod.dec.u,all.dec);
    yhat.renorm <- t(t(postprobs.u)%*%xbeta.all[ind.match,]);
    yhat.bias.rs.rm <- yhat.bias.rs.rm+(yhat.renorm-yhatbma.true);
    yhat.mse.rs.rm <- yhat.mse.rs.rm+(yhat.renorm-yhatbma.true)^2;
    rm(sim.rs); rm(sampmod); rm(sampmod.dec); rm(ind.match);rm(yhat.mc);
    rm(dup); rm(sampmod.dec.u); rm(postprobs.u); rm(yhat.renorm);
    #RS-Thin-MC
    sim.rst <- read.table(file=paste("sim-rs-thin",i,"dat",sep="."))
    sampmod <- as.matrix(sim.rst[,-16]);
    sampmod.dec <- apply(sampmod,1,bindec);
    ind.match <- match(sampmod.dec,all.dec);
    yhat.mc <- apply(xbeta.all[ind.match,],2,mean);
    yhat.bias.rst.mc <- yhat.bias.rst.mc+(yhat.mc-yhatbma.true);
    yhat.mse.rst.mc <- yhat.mse.rst.mc+(yhat.mc-yhatbma.true)^2;
    # RS-Thin-RM
    dup <- duplicated(sampmod.dec);
    sampmod.dec.u <- sampmod.dec[dup==FALSE];
    postprobs.u <- exp(sim.rst[dup==FALSE,16])/sum(exp(sim.rst[dup==FALSE,16]));
    ind.match <- match(sampmod.dec.u,all.dec);
    yhat.renorm <- t(t(postprobs.u)%*%xbeta.all[ind.match,]);
    yhat.bias.rst.rm <- yhat.bias.rst.rm+(yhat.renorm-yhatbma.true);
    yhat.mse.rst.rm <- yhat.mse.rst.rm+(yhat.renorm-yhatbma.true)^2;
    rm(sim.rst); rm(sampmod); rm(sampmod.dec); rm(ind.match);rm(yhat.mc);
    rm(dup); rm(sampmod.dec.u); rm(postprobs.u); rm(yhat.renorm);   
  }


yhat.mse.mc3.rm <- yhat.mse.mc3.rm/nsim
yhat.mse.rs.rm <- yhat.mse.rs.rm/nsim
yhat.mse.rst.rm <- yhat.mse.rst.rm/nsim

yhat.bias.mc3.rm <- yhat.bias.mc3.rm/nsim
yhat.bias.rs.rm <- yhat.bias.rs.rm/nsim
yhat.bias.rst.rm <- yhat.bias.rst.rm/nsim

yhat.mse.mc3.mc <- yhat.mse.mc3.mc/nsim
yhat.mse.rs.mc <- yhat.mse.rs.mc/nsim
yhat.mse.rst.mc <- yhat.mse.rst.mc/nsim

yhat.bias.mc3.mc <- yhat.bias.mc3.mc/nsim
yhat.bias.rs.mc <- yhat.bias.rs.mc/nsim
yhat.bias.rst.mc <- yhat.bias.rst.mc/nsim
  
### BIAS and MSE for MCMC(MC^3,RS,RS-Thin) estimates ####

yhat.bias.vec <- sqrt(c(mean(yhat.bias.bas.ep^2),mean(yhat.bias.bas.unif^2),mean(yhat.bias.mc3.mc^2), mean(yhat.bias.mc3.rm^2),
               mean(yhat.bias.rs.mc^2), mean(yhat.bias.rs.rm^2), mean(yhat.bias.rst.mc^2),mean(yhat.bias.rst.rm^2),
                        mean(yhat.bias.srs^2)))

yhat.rmse.vec <- sqrt(c(mean(yhat.mse.bas.ep),mean(yhat.mse.bas.unif),mean(yhat.mse.mc3.mc),mean(yhat.mse.mc3.rm),mean(yhat.mse.rs.mc),
                       mean(yhat.mse.rs.rm),mean(yhat.mse.rst.mc),mean(yhat.mse.rst.rm),mean(yhat.mse.srs)))

round(yhat.bias.vec*10^(3),2)

round(yhat.rmse.vec*10^(3),2)


yhat.bias.vec.mc <- sqrt(c(mean(yhat.bias.mc3.mc^2),mean(yhat.bias.rs.mc^2),mean(yhat.bias.rst.mc^2)))

yhat.rmse.vec.mc <- sqrt(c(mean(yhat.mse.mc3.mc),mean(yhat.mse.rs.mc),mean(yhat.mse.rst.mc)))

round(yhat.bias.vec.mc*10^(3),2)

round(yhat.rmse.vec.mc*10^(3),2)
### BIAS and MSE for MCMC(MC^3,RS,RS-Thin) estimates ####

# 3. End: Calculate BIAS & MSE for BMA estimate of Y (Tables 1 and 2)

# 4. Begin: Calculate CVRMSE for protein data (Table 3) #######

# Begin: Calculate RS-Thin CV
protein.full.uncen <- read.table("proteinuncen.withnames.txt",header=TRUE) # uncentered protein data for n=96 cases
nrow.protein <- nrow(protein.full.uncen);
n.top <- 10000;         # n.top: no. of top models to be used for BMA 
g.pred <- 95;  # value of g in g-prior used for prediction                                                      
burnin <- 10000;
nsim <- 96;
cv.rst.mc <- rep(NA,nsim);
cv.rst.rm <- rep(NA,nsim);

for(i in 1:nsim)
  {
   print(i);
  # Read RST samples for ith case deleted
   name <- paste("rsthin-pred",i,"dat",sep=".")
   protein.rst <- read.table(name);
   protein.rst <- protein.rst[-c(1:burnin),]
   # Determine n.top models out of 2^{20} post-burnin
   dup.rst <- duplicated(protein.rst[,-89]);
   lmargunique.rst <- protein.rst[dup.rst==FALSE,89];
   postprobs.rst <- exp(lmargunique.rst)/sum(exp(lmargunique.rst));
   sampmod.rst.u <- protein.rst[dup.rst==FALSE,-89]; 
   best.ord <-  order(-postprobs.rst); # positions of ordered models (best to worst)
   best.ord.top <- best.ord[1:n.top];
   models.top <- sampmod.rst.u[best.ord.top,];
# Read ith case deleted (centered afterwards) data & find Bayes estimates for top models
 proteindata.name <- paste("proteindata",i,"cen","txt",sep=".");
 proteindata <- read.table(file=proteindata.name);
 design.colmean <- mean(protein.full.uncen[-i,-89]);
 olsmat <- apply(models.top,1,function (x) coef(lm(proteindata[,89]~as.matrix(proteindata[,-89][,x==1]))))

# Predict ith case based on BMA-RM for top models
   postprobs.top <- postprobs.rst[best.ord.top];
   postprobs.top <- postprobs.top/sum(postprobs.top);
 rm(dup.rst);rm(lmargunique.rst);rm(postprobs.rst); rm(sampmod.rst.u);
   
 Ypredvec <-  rep(NA,n.top);
     
  for(j in 1:n.top)
    {
      alphapost <- olsmat[[j]][1]; betagpost <- g.pred*(olsmat[[j]][-1])/(g.pred+1);
      Ypredvec[j] <- alphapost+sum((protein.full.uncen[i,-89][,models.top[j,]==1]-design.colmean[models.top[j,]==1])*betagpost);
    }

  cv.rst.rm[i] <- sum(Ypredvec*postprobs.top);

# Predict ith case based on BMA-MC for top models

    gamma.char <- apply(protein.rst[,-89],1,function(x) paste(x,collapse="")); 
    tab.char <- table(gamma.char);
    prob.mc <- tab.char/sum(tab.char);
    models.mc <- names(tab.char);
    dup <- duplicated(gamma.char)
    gamma.char.u <- gamma.char[dup==FALSE];
    samp.posi <- match(gamma.char.u,models.mc);
    gamma.u <- protein.rst[dup==FALSE,-89]; #same as sampmod.rst.u
    prob.mc.u <- prob.mc[samp.posi];

   
   postprobs.top <- prob.mc.u[best.ord.top];
   postprobs.top <- postprobs.top/sum(postprobs.top);
 rm(gamma.char); rm(gamma.char.u); rm(tab.char); rm(prob.mc);
 rm(models.mc); rm(samp.posi); rm(gamma.u);  
   
 Ypredvec <-  rep(NA,n.top);
     
  for(j in 1:n.top)
    {
      alphapost <- olsmat[[j]][1]; betagpost <- g.pred*(olsmat[[j]][-1])/(g.pred+1);
      Ypredvec[j] <- alphapost+sum((protein.full.uncen[i,-89][,models.top[j,]==1]-design.colmean[models.top[j,]==1])*betagpost);
    }
      cv.rst.mc[i] <- sum(Ypredvec*postprobs.top);
   }

# End: Calculate RS-Thin CV

# Load BAS/SRS cv results 
 load(file="cv.bas.ep.Rdata")
 load(file="cv.bas.unif.Rdata")
 load(file="cv.srs.Rdata")
 load(file="cv.bas.rst.mc.Rdata")
 load(file="cv.bas.rst.rm.Rdata")

# Calculate cvrmse as in Table 3
cvrmse <- data.frame(matrix(NA,nrow=1,ncol=7))
names(cvrmse) <- c("SRSWOR","BAS-uniform","BAS-eplogp","BAS-RST-MC","BAS-RST-RM","RST-MC","RST-RM")


cvrmse$SRSWOR <- sqrt(mean((cv.srs-protein.full.uncen$prot.act4)^2))
cvrmse$"BAS-uniform" <- sqrt(mean((cv.bas.unif-protein.full.uncen$prot.act4)^2))
cvrmse$"BAS-eplogp" <- sqrt(mean((cv.bas.ep-protein.full.uncen$prot.act4)^2))
cvrmse$"BAS-RST-MC" <- sqrt(mean((cv.bas.rst.mc-protein.full.uncen$prot.act4)^2))
cvrmse$"BAS-RST-RM" <- sqrt(mean((cv.bas.rst.rm-protein.full.uncen$prot.act4)^2))
cvrmse$"RST-MC" <- sqrt(mean((cv.rst.mc-protein.full.uncen$prot.act4)^2))
cvrmse$"RST-RM" <- sqrt(mean((cv.rst.rm-protein.full.uncen$prot.act4)^2))

round(cvrmse,2)

# 4. End: Calculate CVRMSE for protein data (Table 3) #######
