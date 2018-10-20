rm(list=ls())


library(BAS) # may need to specify path via "lib.loc" if installed locally

# Begin: Figure 2a
load("bas.100.Rdata");load("srs.100.Rdata");
#postscript("boxbas.multi.ps",height=6,width=4)
boxplot(as.data.frame(bas.100$trueprob),xlab="Percent Sampled",ylab="True Probability of Unsampled Models",
        ylim=c(0,1),main="BAS",names=as.character(c(1,2,3,4,5,6,7,8,9,10)),cex.lab=1.5,cex.axes=1.5)
#dev.off()
# End: Figure 2a

# Begin: Figure 2b
#postscript("boxsrs.multi.ps",height=6,width=4)
boxplot(as.data.frame(srs.100$trueprob[,1:10]),xlab="Percent Sampled",ylab="True Probability of Unsampled Models",main="SRSWOR",ylim=c(0,1),names=as.character(c(1,2,3,4,5,6,7,8,9,10)),cex.lab=1.5,cex.axes=1.5)
#dev.off()
# End: Figure 2b


# Begin: Figure 3
nsim <- 100;
simx <- read.table("simcen-x.txt")
simy <- read.table("simcen-y.txt")
multidata <- cbind(simx,simy);
names(multidata)[ncol(simx)+1] <- "y";

 multi.all <- bas.lm(y~.,data=multidata,n.models=2^15,prior="g-prior",update=500,
                     alpha=100,initprobs="eplogp")
 truesummarg <- sum(exp(multi.all$logmarg))
 trueprob.bas.ep <- rep(NA,nsim);
 trueprob.bas.unif <- rep(NA,nsim);
 trueprob.rs <- rep(NA,nsim);
 trueprob.rst <- rep(NA,nsim);
 trueprob.mc3 <-rep(NA,nsim);
  
for (isim in 1:nsim)
  {
    print(isim);
    # BAS-ep
   load(file=paste("sim.bas.ep",isim,"Rdata",sep="."))
   marg.bas.ep <- exp(sim.bas.ep$logmarg);
   trueprob.bas.ep[isim] <- 1-(sum(marg.bas.ep)/truesummarg);
   rm(multi.bas.ep); rm(marg.bas.ep);
    # BAS-unif
   load(file=paste("sim.bas.unif",isim,"Rdata",sep="."))
   marg.bas.unif <- exp(sim.bas.unif$logmarg);
   trueprob.bas.unif[isim] <- 1-(sum(marg.bas.unif)/truesummarg);
   rm(multi.bas.unif); rm(marg.bas.unif);
     # RS
  rs <- read.table(paste("sim-rs",isim,"dat",sep="."));
  samp <- rs[,-16];
  marg <- exp(rs[,16]);
  dup <- duplicated(samp);
  marg.unique <- marg[dup==FALSE];
  trueprob.rs[isim] <- 1-(sum(marg.unique)/truesummarg);
  rm(rs); rm(samp); rm(marg); rm(dup); rm(marg.unique);
     # RS-Thin
  rst <- read.table(paste("sim-rs-thin",isim,"dat",sep="."));
  samp <- rst[,-16];
  marg <- exp(rst[,16]);
  dup <- duplicated(samp);
  marg.unique <- marg[dup==FALSE];
  trueprob.rst[isim] <- 1-(sum(marg.unique)/truesummarg);
  rm(rst); rm(samp); rm(marg); rm(dup); rm(marg.unique);
  # MC^3
  mc3.itn <- read.table(paste("sim-mc3",isim,"dat",sep="."));
  samp <- mc3.itn[,-16];
  marg <- exp(mc3.itn[,16]);
  dup <- duplicated(samp);
  marg.unique <- marg[dup==FALSE];
  trueprob.mc3[isim] <- 1-(sum(marg.unique)/truesummarg); 
  rm(mc3.itn); rm(samp); rm(marg); rm(dup); rm(marg.unique);     
}

View(t(1-c(mean(trueprob.bas.ep),mean(trueprob.bas.unif),mean(trueprob.mc3),mean(trueprob.mc3),mean(trueprob.rs),mean(trueprob.rs),mean(trueprob.rst),mean(trueprob.rst))))
#postscript("box100p15bimodal.itns.ps",width=6,height=6)
boxplot(data.frame(cbind(trueprob.bas.ep,trueprob.bas.unif,trueprob.mc3,trueprob.rs,trueprob.rst)),
        names=rep(NA,5),cex.lab=1.3,cex.axes=1.3,ylab="True Probability of Unsampled Models")
mtext(side=1,line=c(2),at=c(1:2,5),text=c("BAS\neplogp","BAS\nuniform","RS\nThin"));
mtext(side=1,lin=1,at=3:4,text=c(expression(MC^3),"RS"));
title("3276 Iterations",cex.main=1.5)
#dev.off()
# End: Figure 3

# Begin: Figure 4
 load("crime.bas.ep.Rdata")
 crime.rst <- read.table("crime-rs-thin.1.dat")
 sample.index <- c(1:(2^15));
 crime.rst.lmarg <- crime.rst[,16];
 dup.rst <- duplicated(crime.rst[,-c(16)]);
 crime.rst.lmarg.unique <- crime.rst.lmarg[dup.rst==FALSE]
 samp.ind.rst.unique <-  sample.index[dup.rst==FALSE]
 crime.rst.lmarg.rep <- crime.rst.lmarg[dup.rst==TRUE]
 samp.ind.rst.rep <-  sample.index[dup.rst==TRUE]

#postscript("uscrime.bas.ep.enumerate.ps",horizontal=TRUE,width=13,height=3)
 plot(x=1:2^15,y=crime.bas.ep$logmarg,pch=".",xlab="Iterations",ylab="log(Marg.Lik)",main="BAS-eplogp",
      cex.lab=1.7,cex.axis=1.4,cex.main=1.3,omd=c(0,0.7,0,0.7))
#dev.off()

# postscript("../output/crime-out/crime.rst.traceplot.ps",width=13,height=3,horizontal=TRUE)
 plot(x=sample.index,y=crime.rst.lmarg,type="n",xlab="Iterations",ylab="log(Marg.Lik)",main="RS-Thin",
      cex.lab=1.7,cex.axis=1.4,cex.main=1.3,omd=c(0,0.7,0,0.7))
 points(x=samp.ind.rst.rep,y=crime.rst.lmarg.rep,col="grey",pch=".")
 points(x=samp.ind.rst.unique,y=crime.rst.lmarg.unique,pch=".")                                        
# dev.off()
1-sum(crime.rst.marg)/norm.const    # 0.03682799: unsampled mass
# End: Figure 4

# Begin: Figure 5

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

sortunique.rsthin.full <- vector("list",10)
sortunique.bas.ep <- vector("list",10)
sortunique.bas.unif <- vector("list",10)
sortunique.bas.rst.mc <- vector("list",10) # no comb just BAS samples
sortunique.bas.rst.rm <- vector("list",10) # no comb just BAS samples
masses<- array(data = 0,dim = c(10,5))

 for (i in 1:10)
  {
    print(i)
     # RS-Thin
    name1 <- paste("protein-rs-thin",i,"dat",sep=".");
    protein.rst <- read.table(name1,header=F);
    dup <- duplicated(protein.rst[,-89]);
    lmargrst.unique <- protein.rst[dup==FALSE,89];
    lmargrst.ord <-  sort(lmargrst.unique,decreasing=T);
    sortunique.rsthin.full[[i]] <- lmargrst.ord;
    masses[i,3]<-sum(exp(sortunique.rsthin.full[[i]]))
    rm(dup);
    rm(name1); rm(protein.rst);
    rm(lmargrst.unique);rm(lmargrst.ord);
    # BAS-RST-MC  
    name2 <- paste("protein.rstmchalf.bas",i,"Rdata",sep=".")
    load(name2);
    sortunique.bas.rst.mc[[i]] <- sort(protein.rstmchalf.bas$logmarg,decreasing=TRUE)
    masses[i,5]<-sum(exp(sortunique.bas.rst.mc[[i]]))
    rm(name2);rm(protein.rstmchalf.bas);
    # BAS-RST-RM  
    name2 <- paste("protein.rstrmhalf.bas",i,"Rdata",sep=".")
    load(name2);
    sortunique.bas.rst.rm[[i]] <- sort(protein.rstrmhalf.bas$logmarg,decreasing=TRUE)
    masses[i,4]<-sum(exp(sortunique.bas.rst.rm[[i]]))
    rm(name2);rm(protein.rstrmhalf.bas);
    # BAS-eplogp
    name <- paste("protein.ep.bas",i,"Rdata",sep=".")
    load(name);
    sortunique.bas.ep[[i]] <- sort(protein.ep.bas$logmarg,decreasing=TRUE); 
    masses[i,1]<-sum(exp(sortunique.bas.ep[[i]]))
    rm(name);rm(protein.ep.bas);
      # BAS-unif
    name <- paste("protein.unif.bas",i,"Rdata",sep=".")
    load(name);
    sortunique.bas.unif[[i]] <- sort(protein.unif.bas$logmarg,decreasing=TRUE); 
    masses[i,2]<-sum(exp(sortunique.bas.unif[[i]]))
    rm(name);rm(protein.unif.bas);
  }

masses1 <- cbind(masses[,4],masses[,5],masses[,1],masses[,3],masses[,2])

topmcmc <- function(n.top,unique.sorted)
  {
    topunique <- matrix(NA,n.top,10);
    for(i in 1:10)
      {print(i)
        topunique[,i] <- unique.sorted[[i]][1:n.top]
      }
    topunique <- data.frame(topunique)
    return(topunique)
        }

top100000.rsthin.full <- topmcmc(100000,sortunique.rsthin.full)
top100000.bas.rst.rm <- topmcmc(100000,sortunique.bas.rst.rm)
top100000.bas.rst.mc <- topmcmc(100000,sortunique.bas.rst.mc)
top100000.bas.ep <- topmcmc(100000,sortunique.bas.ep)
top100000.bas.unif <- topmcmc(100000,sortunique.bas.unif)

setwd("/mn/anatu/ansatte-u3/aliaksah/Desktop/Untitled Folder/")

masses.mine<-array(0,10)
mliks<-array(0,c(2^20-1,10))
for (i in 1:10)
{
  print(i)
  mliks[,i]<- sort(read.csv(paste(" mliks",i,".csv",sep = "",collapse = ""))[[1]][1:(2^20-1)],decreasing=TRUE)
  masses.mine[i] <- sum(exp(mliks[,i]))
}

top100000.mjmcmc<-mliks[1:100000,]

top100000.mjmcmc[1:10,]

max(top100000.rsthin.full)
max(top100000.bas.rst.rm,na.rm = T)
max(top100000.bas.rst.mc,na.rm = T)
max(top100000.bas.ep,na.rm = T)
max(top100000.bas.unif,na.rm = T)


top100000.rsthin.full[1:10,]
top100000.bas.rst.rm[1:10,]
top100000.bas.rst.mc[1:10,]
top100000.bas.ep[1:10,]
top100000.bas.unif[1:10,]

mat <- cbind(top100000.mjmcmc,NA,top100000.bas.ep,NA,top100000.bas.unif,NA,top100000.rsthin.full,NA,top100000.bas.rst.rm,NA,top100000.bas.rst.mc)


write(mliks[,6], file = paste(" mliks6.csv"),
      ncolumns = 1,
      append = FALSE, sep = " ")
write(pp6.rs, file = paste(" pp6.rs.csv"),
      ncolumns = 1,
      append = FALSE, sep = " ")
write(pp6.mc, file = paste(" pp6.mc.csv"),
      ncolumns = 1,
      append = FALSE, sep = " ")

dim(mat)
#postscript("proteintop100000.hybrid.bash.basep.ps",width=13,height=5,horizontal=TRUE)
boxplot(mat[,1:65],xaxt="n",ylab="log(Marginal Likelihood)",xlab="Replicates",ylim=c(27,42),horizontal=FALSE,pch=".",cex.lab=1.7,cex.axis=1.5,
        omd=c(0,0.7,0,0.7))
mtext(side=3,at=seq(6,6*11,by=11),line=0.5,las=1,text=c("MJMCMC ","BAS\neplogp","BAS\nuniform","RST","BAS\nRST-RM","BAS\nRST-MC"),cex=1.3)
axis(side=1,at=(1:65)[-11*(1:4)],labels=FALSE,tick=TRUE,cex=0.2)
at.vec <- c(seq(1,10,by=2),seq(12,21,by=2),seq(23,32,by=2),seq(34,43,by=2),seq(45,54,by=2),seq(56,65,by=2))
mtext(side=1,line=0.65,at=at.vec,text=rep(c(1,3,5,7,9),5),las=1,cex=1.3)
clip(x1 = 0, x2 = 100, y1 = 26, y2 = 43) 
abline(v=11*(1:5))

masses<-cbind(masses.mine,masses)

#postscript("proteintop100000.hybrid.bash.basep.ps",width=13,height=5,horizontal=TRUE)
boxplot(masses,xaxt="n",ylab="sum(exp(Marginal Likelihoods))",horizontal=FALSE,pch=".",cex.lab=1.7,cex.axis=1.5,
        omd=c(0,0.7,0,0.7))
mtext(side=1,at=seq(1,6,by=1),line=2,las=0,text=c("MJMCMC ","BAS\neplogp","BAS\nuniform","RST","BAS\nRST-RM","BAS\nRST-MC"),cex=1.3)

axis(side=1,at=(1:65)[-11*(1:4)],labels=FALSE,tick=TRUE,cex=0.2)
at.vec <- c(seq(1,10,by=2),seq(12,21,by=2),seq(23,32,by=2),seq(34,43,by=2),seq(45,54,by=2),seq(56,65,by=2))
mtext(side=1,line=0.65,at=at.vec,text=rep(c(1,3,5,7,9),5),las=1,cex=1.3)
clip(x1 = 0, x2 = 100, y1 = 27, y2 = 42) 
abline(v=11*(1:5))

setwd("/mn/anatu/ansatte-u3/aliaksah/Desktop/Untitled Folder/")
pp6.rs<-t(as.matrix(read.csv(paste(" pp6.rs",".csv",sep = "",collapse = ""))))
pp6.mc<-t(as.matrix(read.csv(paste(" pp6.mc",".csv",sep = "",collapse = ""))))
barplot(pp6.mc,density = 46,border="black",main = "Marginal Inclusion (MC)",xlab = "Covariates",ylab="Probability")
barplot(pp6.rs,density = 46,border="black",main = "Marginal Inclusion (RM)",xlab = "Covariates",ylab="Probability")
barplot(x2,) 
barplot()
#dev.off()
# End: Figure 5
barplot(pp6.mc,density = 71,border="grey")

barplot(pp6.rs,border="black",density=99)

colMeans(masses)

# Begin: Figure 6
burnin <- 10000;

nsim <- 10; p <- 88;
incprob.srs <- matrix(NA,nsim,p)
incprob.basep <- matrix(NA,nsim,p);
incprob.basunif <- matrix(NA,nsim,p);
incprob.bas.rstmc <-  matrix(NA,10,88);
incprob.bas.rstrm <-  matrix(NA,10,88);
incprob.rstmc <- matrix(NA,10,88);
incprob.rstrm <- matrix(NA,10,88);


for (i in 1:10)
  {
    print(i);
    #RST-RM
    name <- paste("protein-rs-thin",i,"dat",sep=".")
    protein.rst <- read.table(name,header=F);
    dup.rst <- duplicated(protein.rst[,-89])
    lmarg.rst.unique <- protein.rst[dup.rst==FALSE,89]
    postprobs.rst <- exp(lmarg.rst.unique)/sum(exp(lmarg.rst.unique));
    sampmod.rst <- protein.rst[dup.rst==FALSE,-89]
    incprob.rst <- postprobs.rst%*%as.matrix(sampmod.rst);
    incprob.rstrm[i,] <- incprob.rst;
    rm(dup.rst);rm(name);rm(sampmod.rst);
    rm(lmarg.rst.unique); rm(postprobs.rst);   
    # RST-MC
    incprob.rstmc[i,] <- apply(protein.rst[-c(1:burnin),-89],2,mean);
    rm(protein.rst);
    # BAS-RST-MC
    name2 <- paste("protein.rstmchalf.bas",i,"Rdata",sep=".")
    load(name2);
    incprob.bas.rstmc[i,] <- protein.rstmchalf.bas$probne0[-1];
    rm(protein.rstmchalf.bas);
    # BAS-RST-RM
    name2 <- paste("protein.rstrmhalf.bas",i,"Rdata",sep=".")
    load(name2);
    incprob.bas.rstrm[i,] <- protein.rstrmhalf.bas$probne0[-1];
    rm(protein.rstrmhalf.bas);
    # BAS-eplogp
    name <- paste("protein.ep.bas",i,"Rdata",sep=".")
    load(name);
    incprob.basep[i,] <- protein.ep.bas$probne0[-1];
    rm(protein.ep.bas);
    # BAS-unif
    name <- paste("protein.unif.bas",i,"Rdata",sep=".")
    load(name);
    incprob.basunif[i,] <- protein.unif.bas$probne0[-1];
    rm(protein.unif.bas);
    # SRSWOR
     name <- paste("protein.srs",i,"Rdata",sep=".")
    load(name);
    incprob.srs[i,] <- protein.srs$probne0[-1];
    rm(protein.srs);
  }


SRSWOR <- apply(incprob.srs,2,mean)
RST.MC <- apply(incprob.rstmc,2,mean)
RST.RM <- apply(incprob.rstrm,2,mean)
BAS.RST.MC <- apply(incprob.bas.rstmc,2,mean)
BAS.RST.RM <- apply(incprob.bas.rstrm,2,mean)
BAS.unif <- apply(incprob.basunif,2,mean)

incprobs <- data.frame(SRSWOR,BAS.unif,RST.MC,BAS.RST.MC,RST.RM,BAS.RST.RM)
names(incprobs)[2:6] <- c("BAS-unif","RST-MC","BAS-RST-MC","RST-RM","BAS-RST-RM")
names(incprobs)

panel.hist <- function(x, ...)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(usr[1:2], 0, 1.5) )
         h <- hist(x, plot = FALSE)
         breaks <- h$breaks; nB <- length(breaks)
         y <- h$counts; y <- y/max(y)
         rect(breaks[-nB], 0, breaks[-1], y, col="grey", ...)
     }

panel.quad <- function(x,y, ...){
        points(x, y, xlim=c(0,1), ylim=c(0,1), ...);
        abline(0,1);
        abline(h=.5, lty=2); abline(v=.5,lty=2)
      }

textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) mtext(txt, cex = .75, font = font, side=3)


postscript("protein-incprobs.ps", height=11, width=11, horizontal=T)
pairs(incprobs, xlim=c(0,1), ylim=c(0,1),
    panel=panel.quad, text.panel=textPanel, cex=.75,
    # diag.panel=panel.hist
    diag.panel = function(x, ...) {
      lines((1:length(x))/length(x), x, type="h", ...)
    })
dev.off()
# End: Figure 6
