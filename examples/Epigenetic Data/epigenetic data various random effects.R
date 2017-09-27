library("RCurl")
#define the working directory

workdir<-""

# get the data
M<-5
size<-1

data.example <- read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/Epigenetic%20Data/epigen.txt"),sep = ",",header = T)[,2:30]
data.example<-data.example[sample.int(dim(data.example)[1],200),]


data.example$pos1 = data.example$pos
data.example$pos2 = data.example$pos
data.example$pos3 = data.example$pos

fparams <-c(colnames(data.example )[c(8:10,12:17,21:24,29)],"f(data.example$pos,model=\"ar1\")","f(data.example$pos1,model=\"rw1\")","f(data.example$pos2,model=\"iid\")","f(data.example$pos3,model=\"ou\")")

fobservs <- colnames(data.example)[5]
#create MySearch object with default parameters. N/B default estimator is INLA!
args<-list(family = "poisson",control.compute = list(dic = TRUE, waic = TRUE, mlik = TRUE))


system.time({

  formula1 = as.formula(paste(fobservs,"~ 1 +",paste0(fparams,collapse = "+")))

  res = runemjmcmc(formula = formula1,data = data.example,latnames = c("f(data.example$pos,model=\"ar1\")","f(data.example$pos1,model=\"rw1\")","f(data.example$pos2,model=\"iid\")","f(data.example$pos3,model=\"ou\")"),recalc_margin = 199,estimator =estimate.inla.poisson,estimator.args =  list(data=data.example),save.beta = F,interact = T,relations = c("sin","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=3,mutation_rate = 200, max.tree.size = 200000, Nvars.max = 15,p.allow.replace=0.7,p.allow.tree=0.1,p.nor=0.3,p.and = 0.7),n.models = 10000,unique = T,max.cpu = 4,max.cpu.glob = 2,create.table = F,create.hash = T,pseudo.paral = F,burn.in = 100,print.freq = 10,advanced.param = list(
    max.N.glob=as.integer(10),
    min.N.glob=as.integer(5),
    max.N=as.integer(1),                                                                                                                                                                                                                                                                                                      min.N=as.integer(1),
    printable = F))
  print(res$p.post)
})

res=mySearch$post_proceed_results_hash(hashStat)

res$p.post
mySearch$fparam

# specify some INLA realted parameters
mySearch$estimator = inla

data.example$pos1<-data.example$pos
data.example$pos2<-data.example$pos
data.example$pos3<-data.example$pos

lambda = c(inla.pc.ar.lambda(p = 2, b = 0.5), rep(1, 10))
initial = c(inla.models()$latent$ar$hyper$theta2$to.theta(pacf), rep(0, 10))
#example of the underlying model within INLA
formula2 <- as.formula("methylated_bases ~ 1+ (data.example$pos*f(data.example$pos,model=\"ar1\"))+f(data.example$pos2,model=\"iid\")")#  +f(data.example$pos1,model=\"rw1\")+f(data.example$pos2,model=\"iid\")")


# +f(data.example$pos2,model=\"rw2\")+f(data.example$pos3,model=\"crw2\")")
estimate.inla.poisson(formula = formula2,data = data.example)
summary(fm4)

