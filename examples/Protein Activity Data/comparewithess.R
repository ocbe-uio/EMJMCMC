

formula.prot.1 = as.formula(paste0(fobserved.example,"~",paste(fparam.example[c(11,38,46,68)],collapse = "+")))

formula.prot.2 = as.formula(paste0(fobserved.example,"~",paste(fparam.example[c(5,6,10,14,16,20,27,38,46,54,69,72)],collapse = "+")))

formula.prot.3 = as.formula(paste0(fobserved.example,"~",paste(fparam.example[c(3,6,7,8,10,11,16,26,40,58,64,69,76,86)],collapse = "+")))

formula.prot.4 = as.formula(paste0(fobserved.example,"~",paste(fparam.example[c(6,10,11,38,46,69,86)],collapse = "+")))


estimate.bas.lm(formula = formula.prot.1,data = data.example,prior = 3, n = 96, g = 96)$mlik + 149.3995468795

estimate.bas.lm(formula = formula.prot.2,data = data.example,prior = 3, n = 96, g = 96)$mlik + 150.5537667327

estimate.bas.lm(formula = formula.prot.3,data = data.example,prior = 3, n = 96, g = 96)$mlik + 151.3277283846

estimate.bas.lm(formula = formula.prot.4,data = data.example,prior = 3, n = 96, g = 96)$mlik + 136.9203682823



marg.inclusions.ess <- read.table("Z:/ESS/Prot1/Output/prot_Example_3000000_sweeps_output_marg_prob_incl.txt",header = T)[,2]

marg.inclusions.mjm <-  read.table("Z:/ESS/Prot/Output/pp3.rs.csv",header = F)[,1]

plot(marg.inclusions.ess,marg.inclusions.mjm)


models.ess<- read.table("Z:/ESS/Prot1/Output/prot_Example_3000000_sweeps_output_models_history.txt",header = T,sep = ",")[,1]


models.ess<-stri_split_fixed(str = models.ess,pattern = "\t")  

mliks.ess<-array(0,length(models.ess))
mliks.mjm<-array(0,length(models.ess))

for(i in 1:length(models.ess))
{
  mliks.ess[i]<-as.numeric(models.ess[i][[1]][3])
  covs<-stri_split_fixed(str = models.ess[i][[1]][length(models.ess[i][[1]])],pattern = " ")[[1]]
  covs<-covs[-length(covs)]
  formula.prot = as.formula(paste0(fobserved.example,"~",paste(fparam.example[as.numeric(covs)],collapse = "+")))
  
  mliks.mjm[i]<-estimate.bas.lm(formula = formula.prot,data = data.example,prior = 3, n = 96, g = 96)$mlik
}
  

plot(mliks.ess,mliks.mjm)

mean(mliks.ess-mliks.mjm)
var(mliks.ess-mliks.mjm)


max(mliks.ess-mliks.mjm)
min(mliks.ess-mliks.mjm)

cor(mliks.ess,mliks.mjm)



setwd("Z:/EMJMCMC2016/examples/Protein Activity Data/MJMCMC Results")

masses.mine<-array(0,10)
mliks<-array(0,c(2^20-1,10))
for (i in 1:10)
{
  print(i)
  res<-sort(read.csv(paste(" mliks",i,".csv",sep = "",collapse = ""))[[1]][1:(2^20-1)],decreasing=TRUE)
  mliks[1:length(res),i]<-res 
  masses.mine[i] <- sum(exp(mliks[,i]))
}

top100000.mjmcmc<-mliks[1:100000,]

top100000.mjmcmc[1:10,]


top100000.ess<-array(0,dim = c(100000,10))
masses.ess<-array(0,10)

for(i in 1:10)
{
  models.ess<- read.csv(paste0("Z:/ESS/Prot",i,"/Output/prot_Example_30000000_sweeps_output_models_history.txt"),header = T,sep = "\t",na.strings = " ")[,2]
  print(length(models.ess))
  models.ess<- unique(models.ess)
  print(length(models.ess))
  cut<-min(length(models.ess),2^20)
  models.ess<-models.ess[1:cut]
  models.ess<-as.numeric(as.character(models.ess))+174.6478 
  models.ess<-sort(models.ess,decreasing = T)
  top100000.ess[,i]<-models.ess[1:100000]
  masses.ess[i]<-sum(exp(models.ess))
}


mat <- cbind(top100000.mjmcmc,NA,top100000.mjmcmc,NA,top100000.mjmcmc,NA,top100000.mjmcmc,NA,top100000.mjmcmc,NA,top100000.ess)

models.ess<- read.csv(paste0("Z:/ESS/Prot",10,"/Output/prot_Example_3000000_sweeps_output_models_history.txt"),header = T,sep = "\t",na.strings = " ")[,2]+174.6478 
ms<-NULL
mus<-NULL
for(i in 1:1000)
{
  ms<-c(ms,i*3000)
  mus<-c(mus,length(unique(models.ess[1:ms[i]])))
}  

  
plot(x = ms,y=ms)
points(x=ms,y=mus)

ms<-NULL
mas<-NULL
for(i in 1:1000)
{
  ms<-c(ms,i*3000)
  mas<-c(mas,sum(exp(unique(models.ess[1:ms[i]]))))
}  

plot(x = ms,y=mas)
points(x=ms,y=mas)


#postscript("proteintop100000.hybrid.bash.basep.ps",width=13,height=5,horizontal=TRUE)
boxplot(mat[,1:65],xaxt="n",ylab="log(posterior probability)",xlab="Replicates",ylim=c(27,42),horizontal=FALSE,pch=".",cex.lab=1.7,cex.axis=1.5,
        omd=c(0,0.7,0,0.7))
mtext(side=3,at=seq(6,6*11,by=11),line=0.5,las=1,text=c("MJMCMC ","BAS\neplogp","BAS\nuniform","RST","BAS\nRST-RM","ESS"),cex=1.3)
axis(side=1,at=(1:65)[-11*(1:4)],labels=FALSE,tick=TRUE,cex=0.2)
at.vec <- c(seq(1,10,by=2),seq(12,21,by=2),seq(23,32,by=2),seq(34,43,by=2),seq(45,54,by=2),seq(56,65,by=2))
mtext(side=1,line=0.65,at=at.vec,text=rep(c(1,3,5,7,9),5),las=1,cex=1.3)
clip(x1 = 0, x2 = 100, y1 = 26, y2 = 43) 
abline(v=11*(1:5))

masses<-cbind(masses.mine,masses.mine,masses.mine,masses.mine,masses.mine,masses.ess)

#postscript("proteintop100000.hybrid.bash.basep.ps",width=13,height=5,horizontal=TRUE)
boxplot(masses,xaxt="n",ylab="total probability mass",horizontal=FALSE,pch=".",cex.lab=1.7,cex.axis=1.5,
        omd=c(0,0.7,0,0.7))
mtext(side=1,at=seq(1,6,by=1),line=2,las=0,text=c("MJMCMC ","BAS\neplogp","BAS\nuniform","RST","BAS\nRST-RM","ESS"),cex=1.3)

axis(side=1,at=(1:65)[-11*(1:5)],labels=FALSE,tick=TRUE,cex=0.2)
at.vec <- c(seq(1,10,by=2),seq(12,21,by=2),seq(23,32,by=2),seq(34,43,by=2),seq(45,54,by=2),seq(56,65,by=2))
mtext(side=1,line=0.65,at=at.vec,text=rep(c(1,3,5,7,9),6),las=1,cex=1.3)
clip(x1 = 0, x2 = 100, y1 = 27, y2 = 42) 
abline(v=11*(1:6))





