
rm(list=ls())
library(BAS) # May need to install locally using lib.loc

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

# 1. Begin: Calculate BIAS & MSE for inclusion probabilities (Tables 1 and 2)

simx <-  read.table("simcen-x.txt")
simy <- read.table("simcen-y.txt")
multidata <- cbind(simx,simy);
names(multidata)[ncol(simx)+1] <- "y";
system.time(
multi.all <-bas.glm(y~.,data=multidata,n.models=NULL,betaprior=aic.prior(),method = "BAS",family = binomial(),update=2^20,modelprior= uniform(),initprobs="Uniform")
)
save(multi.all,file=paste("sim.all.models",1,"Rdata",sep="."))
multidata$V2<-(multidata$V10+multidata$V14)*multidata$V9
multidata$V5<-(multidata$V11+multidata$V15)*multidata$V12

View(cor(multidata))
pairs(multidata)

formula = as.formula(y~1+ V1+V4+V7+V9+V11+V13+V14+V16+V17+V18+V19)
X <- model.matrix(object = formula,data = multidata)
out <- bayesglm.fit(x = X, y =  multidata[,21], family=binomial(),coefprior=aic.prior())
out$logmarglik
true.incprob <- multi.all$probne0[-1];
truesummarg <- sum(exp(multi.all$logmarg))
fake500 <- sum(exp(sort(multi.all$logmarg,decreasing = TRUE)[1:2^20]-1))/truesummarg;
all <- as.matrix.which(multi.all);
all.dec <- apply(all[,-1],1,bindec);

length(multi.all$logmarg);

which(multi.all$logmarg == max(multi.all$logmarg))

system.time(
  sim.bas.unif <- bas.glm(y~.,data=multidata,n.models=10000, betaprior=aic.prior(),method = "BAS",family = binomial(),update=1000,modelprior= uniform(),initprobs="Uniform")
)
sum(exp(sim.bas.unif$logmarg))/truesummarg
sort(sim.bas.unif$logmarg,decreasing = TRUE)[1:10]
# BAS-MCMC
system.time(
  sim.bas.mcbas <- bas.glm(y~.,data=multidata,n.models=10000,betaprior=aic.prior(),method =  "MCMC+BAS" ,family = binomial(),update=1000,modelprior= uniform(),initprobs="Uniform")
)
sum(exp(sim.bas.mcbas$logmarg))/truesummarg

# MCMC
system.time(
  sim.bas.mcmc <- bas.glm(y~.,data=multidata,n.models=10000,betaprior=aic.prior(),method =  "MCMC" ,family = binomial(),update=1000,modelprior= uniform(),initprobs="Uniform")
)
sum(exp(sim.bas.mcmc$logmarg))/truesummarg

nsim <- 100; p <- 20; burnin <- 500;
#
incprob.rs.rm <- matrix(NA,nsim,p);
incprob.rs.mc <- matrix(NA,nsim,p);
incprob.bas.mc <- matrix(NA,nsim,p);
incprob.bas.unif <- matrix(NA,nsim,p);

for (i in 1:nsim)
{
  print(i);
 
  # BAS-mc
  name <- paste("sim.bas.mcbas",i,"Rdata",sep=".");
  load(name);
  incprob.bas.mc[i,] <- sim.bas.mcbas$probne0[-1];
  # BAS-uniform
  name <- paste("sim.bas.unif",i,"Rdata",sep=".");
  load(name);
  incprob.bas.unif[i,] <- sim.bas.unif$probne0[-1];
  # mcmc
  name <- paste("sim.bas.mcmc",i,"Rdata",sep=".");
  load(name);
  incprob.rs.rm [i,] <- sim.bas.mcmc$probne0[-1]; 
  incprob.rs.mc[i,]<-sim.bas.mcmc$probs.MCMC[-1]
}

mse.bas.mc <- apply((incprob.bas.mc-matrix(true.incprob,nrow=nsim,ncol=p,byrow=TRUE))^2,2,mean);
mse.rs.rm <- apply(( incprob.rs.rm-matrix(true.incprob,nrow=nsim,ncol=p,byrow=TRUE))^2,2,mean);
mse.rs.mc <- apply(( incprob.rs.mc-matrix(true.incprob,nrow=nsim,ncol=p,byrow=TRUE))^2,2,mean);
mse.bas.unif <- apply((incprob.bas.unif-matrix(true.incprob,nrow=nsim,ncol=p,byrow=TRUE))^2,2,mean);


bias.bas.mc <- apply((incprob.bas.mc-matrix(true.incprob,nrow=nsim,ncol=p,byrow=TRUE)),2,mean);
bias.rs.rm <- apply((incprob.rs.rm-matrix(true.incprob,nrow=nsim,ncol=p,byrow=TRUE)),2,mean);
bias.rs.mc <- apply((incprob.rs.mc-matrix(true.incprob,nrow=nsim,ncol=p,byrow=TRUE)),2,mean);
bias.bas.unif <- apply((incprob.bas.unif-matrix(true.incprob,nrow=nsim,ncol=p,byrow=TRUE)),2,mean);

ord.incprob <- order(true.incprob)

#MSE-incprob-matrix
mse.mat <- cbind(mse.bas.unif,mse.bas.mc, mse.rs.rm,mse.rs.mc)
rmse.mat.ord <- sqrt(mse.mat[ord.incprob,]);
View(round(cbind(true.incprob[ord.incprob],rmse.mat.ord*10^2),2))

#BIAS-incprob-matrix
bias.mat <- cbind(bias.bas.unif,bias.bas.mc, bias.rs.rm, bias.rs.mc)
bias.mat.ord <- bias.mat[ord.incprob,]
View(round(cbind(true.incprob[ord.incprob],bias.mat.ord*10^2),2))

# 1. End: Calculate BIAS & MSE for inclusion probabilities (Tables 1 and 2)


# 2. Begin: Calculate BIAS & MSE for model probabilities (Tables 1 and 2)

# Note mprobs.bias.bas actually saves avg bias^2

# BIAS 
mprobs.bias.bas.mc <- rep(0,2^p);
mprobs.bias.bas.unif <- rep(0,2^p);
mprobs.bias.rs.rm <- rep(0,2^p);
mprobs.bias.rs.mc <- rep(0,2^p);

#
mprobs.mse.bas.mc <- rep(0,2^p);
mprobs.mse.bas.unif <- rep(0,2^p);
mprobs.mse.rs.rm <- rep(0,2^p);
mprobs.mse.rs.mc <- rep(0,2^p);

masses <- array(0,4)
unique.mods <- array(0,4)
total.mods <- array(0,4)

for (i in 1:nsim)
{
  print(i);
  # BAS-mc
  load(file=paste("sim.bas.mcbas",i,"Rdata",sep="."))
  models <- as.matrix.which(sim.bas.mcbas);
  models.dec <- apply(models[,-1],1,bindec);
  ind.match <- match(models.dec,all.dec);
  mprobs.all <- rep(0,2^p);
  mprobs.all[ind.match] <- sim.bas.mcbas$postprobs;
  mprobs.bias.bas.mc <- mprobs.bias.bas.mc+(mprobs.all-multi.all$postprobs);
  masses[2]<- masses[2]+sum(exp(sim.bas.mcbas$logmarg))/truesummarg;
  unique.mods[2]<-length(unique(models.dec));
  total.mods[2]<-length(models.dec);
  mprobs.mse.bas.mc <- mprobs.mse.bas.mc+(multi.all$postprobs-mprobs.all)^2;
  # BAS-uniform
  load(file=paste("sim.bas.unif",i,"Rdata",sep="."))
  models <- as.matrix.which(sim.bas.unif);
  models.dec <- apply(models[,-1],1,bindec);
  ind.match <- match(models.dec,all.dec);
  mprobs.all <- rep(0,2^p);
  mprobs.all[ind.match] <- sim.bas.unif$postprobs;
  mprobs.bias.bas.unif <- mprobs.bias.bas.unif+(mprobs.all-multi.all$postprobs);
  mprobs.mse.bas.unif <- mprobs.mse.bas.unif+(multi.all$postprobs-mprobs.all)^2;
  masses[1]<-masses[1]+sum(exp(sim.bas.unif$logmarg))/truesummarg;
  unique.mods[1]<-length(unique(models.dec));
  total.mods[1]<-length(models.dec);
  # rs-rm
  load(file=paste("sim.bas.mcmc",i,"Rdata",sep="."))
  models <- as.matrix.which(sim.bas.mcmc);
  models.dec <- apply(models[,-1],1,bindec);
  ind.match <- match(models.dec,all.dec);
  mprobs.all <- rep(0,2^p);
  mprobs.all[ind.match] <- sim.bas.mcmc$postprobs;
  mprobs.bias.rs.rm <- mprobs.bias.rs.rm+(mprobs.all-multi.all$postprobs);
  mprobs.mse.rs.rm <- mprobs.mse.rs.rm+(multi.all$postprobs-mprobs.all)^2;
  masses[3]<-masses[3]+sum(exp(sim.bas.mcmc$logmarg))/truesummarg;
  unique.mods[3]<-length(unique(models.dec));
  total.mods[3]<-(sum(sim.bas.mcmc$freq));
  unique.mods[4]<-length(unique(models.dec));
  total.mods[4]<-(sum(sim.bas.mcmc$freq));
  masses[4]<-masses[3]
  mprobs.all <- rep(0,2^p);
  mprobs.all[ind.match] <- sim.bas.mcmc$freq/total.mods[4]
  mprobs.bias.rs.mc <- mprobs.bias.rs.mc+(mprobs.all-multi.all$postprobs);
  mprobs.mse.rs.mc <- mprobs.mse.rs.mc+(multi.all$postprobs-mprobs.all)^2;
}
 

# BIAS
mprobs.bias.bas.mc <- mprobs.bias.bas.mc/nsim;
mprobs.bias.rs.rm <- mprobs.bias.rs.rm/nsim;
mprobs.bias.rs.mc <- mprobs.bias.rs.mc/nsim;
mprobs.bias.bas.unif <- mprobs.bias.bas.unif/nsim;

#
mprobs.mse.bas.mc <- mprobs.mse.bas.mc/nsim;
mprobs.mse.rs.rm <- mprobs.mse.rs.rm/nsim;
mprobs.mse.rs.mc <- mprobs.mse.rs.mc/nsim;
mprobs.mse.bas.unif <- mprobs.mse.bas.unif/nsim;


mprobs.bias.vec <- sqrt(c(mean(mprobs.bias.bas.unif^2),mean(mprobs.bias.bas.mc^2),
                          mean(mprobs.bias.rs.rm^2),mean(mprobs.bias.rs.mc^2)))

mprobs.rmse.vec <- sqrt(c(mean(mprobs.mse.bas.unif),mean(mprobs.mse.bas.mc),mean(mprobs.mse.rs.rm),mean(mprobs.mse.rs.mc)))

View(t(round(cbind(mod.bias = mprobs.bias.vec*10^(5),mod.rmse = mprobs.rmse.vec*10^(5),mod.masses = masses/nsim,unique.mods,total.mods),2)))

  median(masses[2])

View(t(ord.incprob))
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
