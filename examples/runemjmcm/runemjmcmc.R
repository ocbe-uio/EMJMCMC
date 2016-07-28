library(sp)
library(INLA)
library(parallel)
library(bigmemory)
library(snow)
library(MASS)
library(ade4)
#library(copula)
library(compiler)
library(BAS)
require(stats)


#define your working directory, where the data files are stored
workdir<-"/results"

#prepare data
simx <- read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/Protein%20Activity%20Data/proteincen.txt"),sep = " ")
data.example <- as.data.frame(simx)
names(data.example)[89]="Y"


#fparam <- c("Const",colnames(data)[-1])
fparam.example <- colnames(data.example)[-89]
fobserved.example <- colnames(data.example)[89]

for(i in 1:length(fparam.example))
{
  fparam.example[i]=paste("I(V",i,")",sep = "")
}

system.time({

formula1 = as.formula(paste(colnames(data.example)[89],"~ 1 +",paste0(colnames(data.example)[-89],collapse = "+")))

res = runemjmcmc(formula = formula1,data = data.example,estimator =estimate.bas.lm,estimator.args =  list(data = data.example,prior = 3, g = 96 ,n=96),save.beta = F,n.models = 20000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = T,create.hash = F,pseudo.paral = F,burn.in = 100,print.freq = 100,advanced.param = list(
                                                                                                                                                                                                                                                                                                               max.N.glob=as.integer(20),
                                                                                                                                                                                                                                                                                                               min.N.glob=as.integer(5),
                                                                                                                                                                                                                                                                                                               max.N=as.integer(3),
                                                                                                                                                                                                                                                                                                               min.N=as.integer(1),
                                                                                                                                                                                                                                                                                                               printable = F))
})

View(statistics1[1:1000,])
