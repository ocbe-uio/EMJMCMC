set.seed(040590)
X1= as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))
Y1=rnorm(n = 1000,mean = 1+0.7*(X1$V1*X1$V4) + 0.8896846*(X1$V8*X1$V11)+1.434573*(X1$V5*X1$V9),sd = 1)
X1$Y1=Y1

#specify the initial formula
formula1 = as.formula(paste(colnames(X1)[51],"~ 1 +",paste0(colnames(X1)[-c(51)],collapse = "+")))
data.example = as.data.frame(X1)


#run the inference with robust g prior
n_cores <- parallel::detectCores() - 1
res4G = EMJMCMC:::LogicRegr(formula = formula1,data = data.example,family = "Gaussian",prior = "G",report.level = 0.5,d = 15,cmax = 2,kmax = 15,p.and = 0.9,p.not = 0.01,p.surv = 0.2,ncores = n_cores)
print(res4G$feat.stat)
#run the inference with Jeffrey's prior
res4J = EMJMCMC:::LogicRegr(formula = formula1,data = data.example,family = "Gaussian",prior = "J",report.level = 0.5,d = 15,cmax = 2,kmax = 15,p.and = 0.9,p.not = 0.01,p.surv = 0.2,ncores = n_cores)
print(res4J$feat.stat)
