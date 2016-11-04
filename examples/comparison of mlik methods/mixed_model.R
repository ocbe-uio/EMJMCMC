library(geepack)
library(INLA)

data(seizure)

seizure$id<-unique(1:59)
# following Chib and Jeliazkov we remove observation 49
seizure<-seizure[-49,]
# then we fit Poisson Mixed Longitudinal Regression for this data set
meanp = 0
precp = 0.01

y = c(seizure$base,seizure$y1,seizure$y2,seizure$y3,seizure$y4)
trt = rep(seizure$trt,5)
id  = rep(seizure$id,5)
id2  = rep(seizure$id,5)
const.rand = rep(1,length(seizure$trt)*5)
bsl = c(rep(0,length(seizure$trt)),rep(1,length(seizure$trt)*4))
bsl2 = c(rep(0,length(seizure$trt)),rep(1,length(seizure$trt)*4))
N = length(bsl)*2
id3=rep(seizure$id+max(seizure$id),5)
i = 1:(N/2)
ofs = c(rep(8,length(seizure$trt)),rep(2,length(seizure$trt)*4))
X = data.frame(y, trt, bsl, ofs, const.rand,id,id2,bsl2,id3)


# for M1 as defined in http://www.sciencedirect.com/science/article/pii/S0304407697001085
formula1 = y~  offset(ofs) + I(trt) + I(trt*bsl) +I(bsl) +f(id2, model="iid2d", n=N) +f(id3, bsl, copy = "id2")

# default inla
M11 <- inla(formula = formula1,data=X,family="poisson",
            control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)),
           control.family=list(link="log"),control.predictor=list(compute=T))
  
summary(M11)
M11$mlik

#log marginal-likelihood (integration) -914.6904
#log marginal-likelihood (Gaussian)    -915.6116
# Chib's -915.404
# Chib's-Jeliazkov -915.23
gc()
# tuned inla
M12 <- inla(formula = formula1,data=X,family="poisson",
            control.family=list(link="log"),control.predictor=list(compute=T),
            control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)),
            control.inla = list(strategy = "laplace", npoints = 100, diff.logdens= 0.01,fast =F, int.strategy = "grid",dz=15,interpolator="gaussian")
)
summary(M12)
M12$mlik

#log marginal-likelihood (integration) -915.6603
#log marginal-likelihood (Gaussian)    -915.6116
# Chib's -915.404
# Chib's-Jeliazkov -915.023

# for M2 as defined in http://www.sciencedirect.com/science/article/pii/S0304407697001085
formula2 = y~  offset(ofs) + I(trt) + I(trt*bsl) +I(bsl) +f(id2, model="iid2d", n=N/2)

# default inla
M21 <- inla(formula = formula2,data=X,family="poisson",
            control.family=list(link="log"),control.predictor=list(compute=T),
            control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)))

summary(M21)
M21$mlik
#log marginal-likelihood (integration) -967.9906
#log marginal-likelihood (Gaussian)    -969.6394
# Chib's -969.824

# tuned inla
M22 <- inla(formula = formula2,data=X,family="poisson",
            control.family=list(link="log"),control.predictor=list(compute=T),
            control.fixed=list(mean=c(default=meanp),mean.intercept=meanp,prec.intercept=precp,prec=c(default=precp)),
            control.inla = list(strategy = "laplace", npoints = 100, diff.logdens= 0.01,fast =F, int.strategy = "grid",dz=11,interpolator="gaussian")
)
summary(M22)
M22$mlik
#log marginal-likelihood (integration) -969.9984
#log marginal-likelihood (Gaussian)    -969.6394
# Chib's -969.824


