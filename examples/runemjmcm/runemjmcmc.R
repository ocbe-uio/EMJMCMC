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

#protein data
simx <- read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/Protein%20Activity%20Data/proteincen.txt"),sep = " ")
data.example <- as.data.frame(simx)
names(data.example)[89]="Y"

system.time({

formula1 = as.formula(paste(colnames(data.example)[89],"~ 1 +",paste0(colnames(data.example)[-89],collapse = "+")))

res = runemjmcmc(formula = formula1,data = data.example,recalc_margin = 2^10,estimator =estimate.bas.lm,estimator.args =  list(data = data.example,prior = 3, g = 96 ,n=96),save.beta = T,interact = T,relations = c("","sin","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=2,mutation_rate = 100, max.tree.size = 200000, Nvars.max = 95,p.allow.replace=0.9,p.allow.tree=0.5,p.nor=0.3,p.and = 0.7),n.models = 50000,unique = T,max.cpu = 10,max.cpu.glob = 10,create.table = F,create.hash = T,pseudo.paral = F,burn.in = 100,print.freq = 100,advanced.param = list(
                                                                                                                                                                                                                                                                                                               max.N.glob=as.integer(20),
                                                                                                                                                                                                                                                                                                               min.N.glob=as.integer(5),
                                                                                                                                                                                                                                                                                                               max.N=as.integer(3),
                                                                                                                                                                                                                                                                                                               min.N=as.integer(1),
                                                                                                                                                                                                                                                                                                               printable = F))
})


#NEO asteroid data

#define your working directory, where the data files are stored
workdir<-""

#prepare the test set data
simx <- read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/asteroid%20data/Recognize/NEAs.txt"),sep = ",",header = T,fill=TRUE)
simy <-  read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/asteroid%20data/Recognize/NotNeas8%2B.txt"),sep = ",",header = T,fill=TRUE)
simx$neo<-1
simy$neo<-0
#data.example <- as.data.frame(t(cbind(t(simy[sample.int(size = 400,n = 6621,replace = T),]),t(simx[sample.int(size = 600,n = 14099,replace = T),]))),stringsAsFactors = F)
#data.example$epoch<-factor(data.example$epoch,labels = c(0,1))


data.example <- as.data.frame(t(cbind(t(simy),t(simx))),stringsAsFactors = F)

transform<-colnames(data.example)[-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25)]

for(i in 1:length(transform))
{
  print(i)
  data.example[[transform[i]]]<-as.numeric(as.character(data.example[[transform[i]]]))
}

data.example$esuar<-data.example$eccentricity^2
data.example$asuar<-data.example$absolute_magnitude^2
data.example$rsuar<-data.example$semi_major_axis^2
data.example$rcube<-data.example$semi_major_axis^3
data.example$anoml<-data.example$mean_anomaly*data.example$semi_major_axis
data.example$anoms<-data.example$mean_anomaly*data.example$semi_major_axis^2
data.example$anomv<-data.example$mean_anomaly*data.example$semi_major_axis^3
data.example$perihell<-data.example$argument_of_perihelion*data.example$semi_major_axis
data.example$perihels<-data.example$argument_of_perihelion*data.example$semi_major_axis^2
data.example$perihelv<-data.example$argument_of_perihelion*data.example$semi_major_axis^3
data.example$longitudel<-data.example$longitude_of_the_ascending.node*data.example$semi_major_axis
data.example$longitudes<-data.example$longitude_of_the_ascending.node*data.example$semi_major_axis^2
data.example$longitudev<-data.example$longitude_of_the_ascending.node*data.example$semi_major_axis^3


data.example1<-data.example

#prepare the training set data
simx <- read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/asteroid%20data/Teach/NeoPHA.txt"),sep = ",",header = T,fill=TRUE)
simy <-  read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/asteroid%20data/Teach/NotNeo-Type7.txt"),sep = ",",header = T,fill=TRUE)
simx$neo<-1
simy$neo<-0

data.example <- as.data.frame(t(cbind(t(simy),t(simx))),stringsAsFactors = F)

transform<-colnames(data.example)[-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25)]

for(i in 1:length(transform))
{
  print(i)
  data.example[[transform[i]]]<-as.numeric(as.character(data.example[[transform[i]]]))
}

data.example$esuar<-data.example$eccentricity^2
data.example$asuar<-data.example$absolute_magnitude^2
data.example$rsuar<-data.example$semi_major_axis^2
data.example$rcube<-data.example$semi_major_axis^3
data.example$anoml<-data.example$mean_anomaly*data.example$semi_major_axis
data.example$anoms<-data.example$mean_anomaly*data.example$semi_major_axis^2
data.example$anomv<-data.example$mean_anomaly*data.example$semi_major_axis^3
data.example$perihell<-data.example$argument_of_perihelion*data.example$semi_major_axis
data.example$perihels<-data.example$argument_of_perihelion*data.example$semi_major_axis^2
data.example$perihelv<-data.example$argument_of_perihelion*data.example$semi_major_axis^3
data.example$longitudel<-data.example$longitude_of_the_ascending.node*data.example$semi_major_axis
data.example$longitudes<-data.example$longitude_of_the_ascending.node*data.example$semi_major_axis^2
data.example$longitudev<-data.example$longitude_of_the_ascending.node*data.example$semi_major_axis^3

#define the covariates and theobservations
fparam.example <- colnames(data.example)[-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)]
fobserved.example <- colnames(data.example)[1]

system.time({

  formula1 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)],collapse = "+")))

  res = runemjmcmc(formula = formula1,data = data.example,estimator =estimate.bas.glm,estimator.args =  list(data = data.example,prior = aic.prior(),family = binomial(), logn = log(64)),save.beta = T,interact = T,relations = c("","sin","cos","sigmoid","tanh","atan","erf"),relations.prob =c(2.4,0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=2,mutation_rate =500, max.tree.size = 200000, Nvars.max = 25,p.allow.replace=0.95,p.allow.tree=0.1,p.nor=0.3,p.and = 0.7),n.models = 200000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 100,advanced.param = list(
    max.N.glob=as.integer(10),
    min.N.glob=as.integer(5),
    max.N=as.integer(3),
    min.N=as.integer(1),
    printable = F))
})

data.example1[1:20720,-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)][(is.na(data.example1[1:20720,-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)]))]<-0
g<-function(x)
{
  return((x = 1/(1+exp(-x))))
}


system.time({

res<-mySearch$forecast.matrix(link.g = g,covariates = (data.example1[1:20720,-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)]))$forecast

})

res<-as.integer(res>=0.5)
length(which(res>=0.5))
length(which(res<0.5))
length(res)
length(which(data.example1$neo==1))

(1-sum(abs(res-data.example1$neo),na.rm = T)/20720)

