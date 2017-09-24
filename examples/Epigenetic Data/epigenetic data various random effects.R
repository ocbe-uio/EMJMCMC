rm(list = ls(all = TRUE))
# install the required packges if needed
#install.packages("INLA", repos="http://www.math.ntnu.no/inla/R/testing")
#install.packages("bigmemory")
#install.packages("snow")
#install.packages("Rmpi")
#install.packages("ade4")
#install.packages("sp")
#install.packages("BAS")
#install.packages("https://github.com/aliaksah/EMJMCMC2016/files/270429/EMJMCMC_1.2.tar.gz", repos = NULL, type="source")
#install.packages("RCurl")
#install.packages("hash")

library(hash)
library(RCurl)
library(EMJMCMC)
library(sp)
library(INLA)
library(parallel)
library(bigmemory)
library(snow)
library(MASS)
library(ade4)
library(copula)
library(compiler)
library(BAS)
require(stats)

#define the working directory

workdir<-""

# get the data
M<-5
size<-1

data.example <- read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/Epigenetic%20Data/epigen.txt"),sep = ",",header = T)[,2:30]
workdir = ""

estimate.inla.refchoice<-function(formula,data,refs)
{
  fmla.proc<-as.character(formula)[2:3]
  fobserved <- fmla.proc[1]
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam.tmp <-stri_split_fixed(str = fmla.proc[2],pattern = "+",omit_empty = F)[[1]]
  refs.id <- which(fparam.tmp %in% refs)

}
fparam.example <-c(colnames(data.example )[c(8:10,12:17,21:24,29)],"f(data.example$pos,model=\"ar1\")","f(data.example$pos1,model=\"rw1\")","f(data.example$pos2,model=\"iid\")","f(data.example$pos3,model=\"ou\")")

fobserved.example <- colnames(data.example)[5]
#create MySearch object with default parameters. N/B default estimator is INLA!
args<-list(family = "poisson",data = data.example)
args$control.compute = list(dic = TRUE, waic = TRUE, mlik = TRUE)

system.time({

  formula1 = as.formula(paste(fobserved.example,"~ 1 +",paste0(fparam.example,collapse = "+")))

  res = runemjmcmc(formula = formula1,data = data.example,recalc_margin = 200,estimator =inla,estimator.args =  args,save.beta = F,interact = F,relations = c("","sin","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=2,mutation_rate = 100, max.tree.size = 200000, Nvars.max = 95,p.allow.replace=0.9,p.allow.tree=0.5,p.nor=0.3,p.and = 0.7),n.models = 10000,unique = T,max.cpu = 4,max.cpu.glob = 2,create.table = F,create.hash = T,pseudo.paral = F,burn.in = 100,print.freq = 100,advanced.param = list(
    max.N.glob=as.integer(10),
    min.N.glob=as.integer(5),
    max.N=as.integer(1),                                                                                                                                                                                                                                                                                                      min.N=as.integer(1),
    printable = F))
  print(res$p.post)
})

# specify some INLA realted parameters
mySearch$estimator = inla

data.example$pos1<-data.example$pos
data.example$pos2<-data.example$pos
data.example$pos3<-data.example$pos

lambda = c(inla.pc.ar.lambda(p = 2, b = 0.5), rep(1, 10))
initial = c(inla.models()$latent$ar$hyper$theta2$to.theta(pacf), rep(0, 10))
#example of the underlying model within INLA
formula2 <- as.formula("methylated_bases ~  f(data.example$pos,model=\"ar1\")+f(data.example$pos2,model=\"iid\")")#  +f(data.example$pos1,model=\"rw1\")+f(data.example$pos2,model=\"iid\")")


# +f(data.example$pos2,model=\"rw2\")+f(data.example$pos3,model=\"crw2\")")
fm4<-do.call(inla, c(args,formula = formula2))
summary(fm4)

