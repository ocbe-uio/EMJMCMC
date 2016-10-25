rm(list = ls(all = TRUE))

library(RCurl)
library(INLA)
library(MCMCpack)

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


