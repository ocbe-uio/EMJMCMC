X1<- as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))
X1$V9 <- X1$V9*X1$V8
X1$V4 <- (X1$V5)*X1$V1
X1$V2 <- (X1$V8)
View(cor(X1))
Y1<-rnorm(n = 1000,mean = 1+10*(X1$V1) + 0.8896846*(X1$V8)+1.434573*(X1$V5),sd = 1)
X1$Y1<-Y1

set.seed(16112016)
system.time({

  formula1 = as.formula(paste(colnames(X1)[11],"~ 1 +",paste0(colnames(X1)[-c(11)],collapse = "+")))

  res = runemjmcmc(formula = formula1,data = X1,estimator = estimate.bas.lm,estimator.args =  list(data = data.example,n = 1000, prior = 3, g = 1000),recalc_margin = 250, save.beta = T,interact = F,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=1,mutation_rate = 300,last.mutation = 5000, max.tree.size = 4, Nvars.max =15,p.allow.replace=0.9,p.allow.tree=0.2,p.nor=0,p.and = 0.9),n.models = 1024,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,print.freq = 100,advanced.param = list(
    max.N.glob=as.integer(5),
    min.N.glob=as.integer(3),
    max.N=as.integer(1),
    min.N=as.integer(1),
    printable = F))
})
View(t(values(hashStat)))
write.csv(x = (X1),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/runemjmcm/simdata.csv")
max(values(hashStat)[,1],na.rm = T)

