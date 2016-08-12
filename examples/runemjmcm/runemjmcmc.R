library(sp)
library(INLA)
library(parallel)
library(bigmemory)
library(snow)
library(MASS)
library(ade4)
library(hash)
library(RCurl)
library(compiler)
library(BAS)
require(stats)


#define your working directory, where the data files are stored
workdir<-"/results"

#prepare data
simx <- read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/Protein%20Activity%20Data/proteincen.txt"),sep = " ")
data.example <- as.data.frame(simx)
names(data.example)[89]="Y"

system.time({

formula1 = as.formula(paste(colnames(data.example)[89],"~ 1 +",paste0(colnames(data.example)[-89],collapse = "+")))

res = runemjmcmc(formula = formula1,data = data.example,estimator =estimate.bas.lm,estimator.args =  list(data = data.example,prior = 3, g = 96 ,n=96),save.beta = F,ineract = T,relations = c("","sin","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=2,mutation_rate = 100, max.tree.size = 10000, Nvars.max = 110,p.allow.replace=0.8,p.allow.tree=0.5,p.nor=0.3,p.and = 0.7),n.models = 200000,unique = T,max.cpu = 10,max.cpu.glob = 10,create.table = F,create.hash = T,pseudo.paral = F,burn.in = 100,print.freq = 100,advanced.param = list(
                                                                                                                                                                                                                                                                                                               max.N.glob=as.integer(20),
                                                                                                                                                                                                                                                                                                               min.N.glob=as.integer(5),
                                                                                                                                                                                                                                                                                                               max.N=as.integer(3),
                                                                                                                                                                                                                                                                                                               min.N=as.integer(1),
                                                                                                                                                                                                                                                                                                               printable = F))
})

View(statistics1[1:1000,])
