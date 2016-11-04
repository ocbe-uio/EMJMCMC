rm(list = ls(all = TRUE))

library(INLA)
library(MASS)

df = rbind(Pima.tr,Pima.te)
df$type = as.integer(df$type=="Yes")

formula1 = type ~ 1 + npreg + glu + bmi + ped
formula2 = type ~ 1 + npreg + glu + bmi + ped + age

# standartize the variables
for(i in 1:7)
{
  df[,i]<-(df[,i]-mean(df[,i]))/sqrt(var(df[,i]))
}

meanp<-0
precp<-0.01
n <- 532
M11 <- inla(formula = formula1,data=df,family="binomial",Ntrials=rep(1,n),
           control.family=list(link="logit"),control.predictor=list(compute=T),
           control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
summary(M11)
M11$mlik

M12 <- inla(formula = formula2,data=df,family="binomial",Ntrials=rep(1,n),
           control.family=list(link="logit"),control.predictor=list(compute=T),
           control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
summary(M12)
M12$mlik



meanp<-0
precp<-1
n <- 532
M11 <- inla(formula = formula1,data=df,family="binomial",Ntrials=rep(1,n),
            control.family=list(link="logit"),control.predictor=list(compute=T),
            control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
summary(M11)
M11$mlik

M12 <- inla(formula = formula2,data=df,family="binomial",Ntrials=rep(1,n),
            control.family=list(link="logit"),control.predictor=list(compute=T),
            control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))
summary(M12)
M12$mlik
