rm(list=ls())
install.packages("BAS")

library(BAS)

source("repeatsamp.r")


setwd("/mn/anatu/ansatte-u3/aliaksah/Desktop/competotrs/originalexperiments/")  
########### 1. Run MCMC, BAS and SRSWOR for simulated data  #################
nsim <- 100;

# RS
system("gcc -g -O2 -o sim-rs sim-rs.c -llapack -lblas -lm")

system(paste(getwd(),"/sim-rs",sep = "",collapse = ""))


# RS-Thin
system("gcc -g -O2 -o sim-rs-thin sim-rs-thin.c -llapack -lblas -lm")
system(paste(getwd(),"/sim-rs-thin",sep = "",collapse = ""))


# MC^3
system("gcc -g -O2 -o sim-mc3 sim-mc3.c -llapack -lblas -lm")
system(paste(getwd(),"/sim-mc3",sep = "",collapse = ""))


# BAS (eplogp and uniform)
simx <- read.table("simcen-x.txt")
simy <- read.table("simcen-y.txt")
multidata <- cbind(simx,simy);
names(multidata)[ncol(simx)+1] <- "y";

for (isim in 1:nsim)
{
  print(isim);
  # BAS-eplogp
  sim.bas.ep <- bas.lm(y~.,data=multidata,n.models=3276,prior="g-prior",update=500,alpha=nrow(multidata),modelprior= uniform(),initprobs="eplogp")
  save(sim.bas.ep,file=paste("sim.bas.ep",isim,"Rdata",sep="."))
  # BAS-uniform
  sim.bas.unif <- bas.lm(y~.,data=multidata,n.models=3276,prior="g-prior",update=500,alpha=nrow(multidata),modelprior= uniform(),initprobs="Uniform")
  save(sim.bas.unif,file=paste("sim.bas.unif",isim,"Rdata",sep="."))
  # SRSWOR
  sim.srs <- bas.lm(y~.,data=multidata,n.models=3276,prior="g-prior",update=2^15,alpha=nrow(multidata),initprobs="Uniform")
  save(sim.srs,file=paste("sim.srs",isim,"Rdata",sep="."))
}


# SRSWOR and BAS over a grid 1,2... 10% of 2^15
bas.100 <- repeatsamp.bas(100)
save(bas.100,file="bas.100.Rdata")

srs.100 <- repeatsamp.srs(100)
save(srs.100,file="srs.100.Rdata")


########### 2. Run MCMC and BAS for crime data (2^15=32768 iterations) #################

# RS-Thin
system("gcc -g -O2 -o crime-rs-thin.exe crime-rs-thin.c -llapack -lblas -lm")
system("./crime-rs-thin.exe")


# BAS-eplogp
library(MASS);
data(UScrime);
UScrime[,-2] = log(UScrime[,-2])
 
  
crime.bas.ep = bas.lm(y~.,data=UScrime,n.models=2^15,prior="g-prior",update=500,alpha=nrow(UScrime),initprobs="eplogp")
save(crime.bas.ep,file="crime.bas.ep.Rdata")

########### 2. Run MCMC and BAS for crime data (2^15=32768 iterations) #################

########### 3. Run MCMC and BAS for inference for protein data (2^20 iterations) #################
nsim <- 10;

# RS-Thin
system("gcc -g -O2 -o protein-rs-thin.exe protein-rs-thin.c -llapack -lblas -lm")
system("./protein-rs-thin.exe")

# BAS-ep & BAS-unif & SRSWOR
data(protein)

protein.lm <- lm(prot.act4 ~ 
	buf + buf*pH + buf*NaCl + buf*con + buf*ra +
	 buf*det + buf*MgCl2 + buf*temp + 
	pH + pH*NaCl + pH*con + pH*ra + pH*det +
	 pH*MgCl2 +pH*temp +
   	NaCl + NaCl*con + NaCl*ra + NaCl*det +
	 NaCl*MgCl2 + NaCl*temp +
	con + con*ra + con*det + 
	 con*MgCl2 +con*temp +
	ra + ra*det + ra*MgCl2 + ra*temp +
	det + det*MgCl2 + det*temp +
	MgCl2 + MgCl2*temp + I(NaCl^2) + I(pH^2) + I(con^2) + I(temp^2), 
	data=protein,x=T)

protein.designmat <- protein.lm$x
proteindata <-  data.frame(cbind(protein.designmat[,-1],protein$prot.act4));
names(proteindata)[89] <- "prot.act4"
write.table(proteindata,file="proteinuncen.withnames.txt")


for (isim in 1:nsim)
  {
   print(isim);
  protein.ep.bas = bas.lm(prot.act4 ~.,data=proteindata, n.models=2^(20),prior="g-prior",
  alpha=nrow(proteindata),initprobs="eplogp",update=10000);
name <- paste("protein.ep.bas",isim,"Rdata",sep=".");
save(protein.ep.bas,file=name)
rm(protein.ep.bas)
#
  protein.unif.bas = bas.lm(prot.act4 ~.,data=proteindata, n.models=2^(20),prior="g-prior",
  alpha=nrow(proteindata),initprobs="Uniform",update=10000);
name <- paste("protein.unif.bas",isim,"Rdata",sep=".");
save(protein.unif.bas,file=name)
rm(protein.unif.bas)

}

# BAS-RST-MC and BAS-RST-RM
incprob.rstmc.half.mat <- matrix(NA,nsim,ncol(proteindata)-1);
incprob.rstrm.half.mat <- matrix(NA,nsim,ncol(proteindata)-1);
burnin <- 10000;

 for (isim in 1:nsim)
  {
    print(isim)
    name <- paste("protein-rs-thin.",isim,".dat",sep="")
    protein.rst <- read.table(name,header=F);
    protein.rst <- protein.rst[c(1:2^(19)),]; 
    dup.rst <- duplicated(protein.rst[,-89])
    lmarg.rst.unique <- protein.rst[dup.rst==FALSE,89]
    postprobs.rst <- exp(lmarg.rst.unique)/sum(exp(lmarg.rst.unique));
    sampmod.rst <- protein.rst[dup.rst==FALSE,-89]
    incprob.rst <- postprobs.rst%*%as.matrix(sampmod.rst);
    incprob.rstrm.half.mat[isim,] <- incprob.rst;
    #
    gam.rst <- protein.rst[,-89];
    incprob.rstmc.half.mat[isim,] <- apply(gam.rst[-c(1:burnin),],2,mean);
    #
    rm(dup.rst);rm(name);rm(protein.rst);rm(sampmod.rst);
    rm(lmarg.rst.unique); rm(postprobs.rst);   
  }


incprob.rstmc.half.mat[incprob.rstmc.half.mat>0.975] <- 0.975;
incprob.rstmc.half.mat[incprob.rstmc.half.mat<0.025] <- 0.025;
#
incprob.rstrm.half.mat[incprob.rstrm.half.mat>0.975] <- 0.975;
incprob.rstrm.half.mat[incprob.rstrm.half.mat<0.025] <- 0.025;



for (isim in 1:nsim)
  {
    print(isim)
protein.rstmchalf.bas = bas.lm(prot.act4 ~.,data=proteindata, n.models=2^19, prior="g-prior",alpha=nrow(proteindata),
  initprobs=c(1,incprob.rstmc.half.mat[isim,]),update=10000);
name <- paste("protein.rstmchalf.bas",isim,"Rdata",sep=".");
save(protein.rstmchalf.bas,file=name)
rm(protein.rstmchalf.bas)
    #
protein.rstrmhalf.bas = bas.lm(prot.act4 ~.,data=proteindata, n.models=2^19, prior="g-prior",alpha=nrow(proteindata),
  initprobs=c(1,incprob.rstrm.half.mat[isim,]),update=10000);
name <- paste("protein.rstrmhalf.bas",isim,"Rdata",sep=".");
save(protein.rstrmhalf.bas,file=name)
rm(protein.rstrmhalf.bas)    
      }
########### 3. Run MCMC and BAS for inference for protein data (2^20 iterations) #################


########### 4. Run MCMC and BAS for out of sample predictions for protein data (2^20 iterations) #################
  
  nsim <- 96; p <- 88;
  cv.srs <- rep(NA,nsim);
  cv.bas.ep <- rep(NA,nsim);
  cv.bas.unif <- rep(NA,nsim);
  
  # 4a. Begin: Run BAS/SRS for 2^20 iterations for ith case-deleted dataset, i=1 to nsim
  
  for (i in 1:nsim)  
    {
   print(i);
   data.out = as.matrix(proteindata[i,-89], nrows=1);
   # BAS-eplogp
   bas.ep <- bas.lm(prot.act4 ~.,data=proteindata[-i,], n.models=2^20, prior="g-prior",alpha=96,update=10000,initprobs="eplogp");
   cv.bas.ep[i] = as.numeric(predict(bas.ep,data.out,top=10000)$Ybma);
   rm(bas.ep);
   # BAS-uniform
   bas.unif = bas.lm(prot.act4 ~.,data=proteindata[-i,], n.models=2^20, prior="g-prior",alpha=96,update=10000,initprobs="Uniform");
   cv.bas.unif[i] = as.numeric(predict(bas.unif,data.out,top=10000)$Ybma);
   rm(bas.unif);
   # SRSWOR
   #srs = bas.lm(prot.act4 ~.,data=proteindata[-i,], n.models=2^20, prior="g-prior",alpha=96,update=NULL,initprobs="Uniform");
   #cv.srs[i] = as.numeric(predict(srs,data.out,top=10000)$Ybma);
   #rm(srs); 
   }
  
   save(cv.bas.ep,file="cv.bas.ep.Rdata")
   save(cv.bas.unif,file="cv.bas.unif.Rdata")
   #save(cv.srs,file="cv.srs.Rdata")
  
  sqrt(mean((cv.bas.ep-proteindata$prot.act4)^2))# 
  sqrt(mean((cv.bas.unif-proteindata$prot.act4)^2))#
  sqrt(mean((cv.srs-proteindata$prot.act4)^2))#
  
  # 4a. End: Run BAS/SRS for 2^20 iterations for ith case-deleted dataset, i=1 to nsim
  
  # 4b. Begin: RS-Thin for 2^20 iterations for ith case-deleted dataset, i=1 to nsim
  
  # Create case number for input to C function
  for (i in 1:nsim)
    {
      name <- paste("case",i,"txt",sep=".");
      caseno <- i;
      write.table(caseno,file=name,col.names=F,row.names=F);
      rm(caseno);
    }
  
  # Create case deleted & centered (afterwards) protein data for input to C functions
  protein.design.woint <- protein.designmat[,-1];
  
  for (i in 1:nsim)
    {
      protein.i.cen <- scale(protein.design.woint[-i,],center=T,scale=F);
      proteindata.i.cen <- cbind(protein.i.cen,protein$prot.act4[-i]);
      name <- paste("proteindata",i,"cen","txt",sep=".");
      write.table(proteindata.i.cen,file=name,col.names=F,row.names=F);
      rm(proteindata.i.cen);
    }
  
  # Run RS-Thin for case deleted protein data for every case
  for(i in 1:nsim)
    {
      print(i);
  casename <- paste("case",i,"txt",sep=".")
  system("gcc -g -O2 -o protein-rsthin-pred.exe protein-rsthin-pred.c -llapack -lblas -lm")  # compile C code
  system(paste("./protein-rsthin-pred.exe",casename,sep=" ")) # Run for ith case deleted
    }
  # 4b. End: RS-Thin for 2^20 iterations for ith case-deleted dataset, i=1 to nsim
  
  # 4c. Begin: BAS-RS-Thin for 2^19 iterations for ith case-deleted dataset, i=1 to nsim
  
  # Calculate initial inclusion probs (MC and RM) for BAS using 2^19 samples from RS-Thin
  incprob.cv.rst.mc <- matrix(NA,nrow=nsim,ncol=p);
  incprob.cv.rst.rm <- matrix(NA,nrow=nsim,ncol=p);
  
  burnin <- 10000;
  
  for (i in 1:nsim)
    {
      print(i);
      name <- paste("rsthin-pred",i,"dat",sep=".")
      protein.rst <- read.table(name,header=F);
      protein.rst <- protein.rst[c(1:2^19),];
      gam.rst <- protein.rst[,-89];
      incprob.cv.rst.mc[i,] <- apply(gam.rst[-c(1:burnin),],2,mean);
      #
      dup.rst <- duplicated(protein.rst[,-89]);
      lmarg.rst <- protein.rst[dup.rst==FALSE,89];
      probs.rst <- exp(lmarg.rst)/sum(exp(lmarg.rst));
      samp.rst <- as.matrix(protein.rst[dup.rst==FALSE,-89]);
      incprob.cv.rst.rm[i,] <- t(probs.rst)%*%samp.rst;
      #
      rm(protein.rst); rm(gam.rst); rm(dup.rst);
      rm(lmarg.rst); rm(probs.rst); rm(samp.rst);
    }
  
  # Run BAS-RST
  cv.bas.rst.rm <- rep(NA,nsim);
  samp.prob.rm <- incprob.cv.rst.rm;
  samp.prob.rm[samp.prob.rm<0.025] <- 0.025;  samp.prob.rm[samp.prob.rm>0.975] <- 0.975;
  
  cv.bas.rst.mc <- rep(NA,nsim);
  samp.prob.mc <- incprob.cv.rst.mc;
  samp.prob.mc[samp.prob.mc<0.025] <- 0.025;  samp.prob.mc[samp.prob.mc>0.975] <- 0.975;
  
  
  for (i in 1:nsim)
    {
      print(i);
    data.out = as.matrix(proteindata[i,-89], nrows=1);   
   # BAS-RST-MC  
   bas.mc = bas.lm(prot.act4 ~.,data=proteindata[-i,], n.models=2^19, prior="g-prior",alpha=95,initprobs=c(1,samp.prob.mc[i,]),update=10000);
   cv.bas.rst.mc[i] = as.numeric(predict(bas.mc,data.out,top=10000)$Ybma);
   rm(bas.mc);    
   # BAS-RST-RM   
   bas.rm = bas.lm(prot.act4 ~.,data=proteindata[-i,], n.models=2^19, prior="g-prior",alpha=95,initprobs=c(1,samp.prob.rm[i,]),update=10000);
   cv.bas.rst.rm[i] = as.numeric(predict(bas.rm,data.out,top=10000)$Ybma);
   rm(bas.rm);
    }
  
  save(cv.bas.rst.mc,file="cv.bas.rst.mc.Rdata")
  save(cv.bas.rst.rm,file="cv.bas.rst.rm.Rdata")
  
  sqrt(mean((cv.bas.rst.mc-proteindata$prot.act4)^2)) # 
  sqrt(mean((cv.bas.rst.rm-proteindata$prot.act4)^2)) # 


########### 4. Run MCMC and BAS for out of sample predictions for protein data (2^20 iterations) #################
