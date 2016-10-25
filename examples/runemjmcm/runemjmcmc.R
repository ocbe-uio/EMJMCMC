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
                                                                                                                                                                                                                                                                                                               max.N=as.integer(3),                                                                                                                                                                                                                                                                                                      min.N=as.integer(1),
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

  res = runemjmcmc(formula = formula1,data = data.example,estimator =estimate.bas.glm,estimator.args =  list(data = data.example,prior = aic.prior(),family = binomial(), logn = log(64)),recalc_margin = 50, save.beta = T,interact = F,relations = c("","sin","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=2,mutation_rate = 200, max.tree.size = 200000, Nvars.max = 30,p.allow.replace=0.10,p.allow.tree=0.1,p.nor=0.3,p.and = 0.7),n.models = 10000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 100,advanced.param = list(
    max.N.glob=as.integer(10),
    min.N.glob=as.integer(5),
    max.N=as.integer(3),
    min.N=as.integer(1),
    printable = F))
})


ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
ppp$p.post

mySearch$g.results[,]
mySearch$fparam

g<-function(x)
{
  return((x = 1/(1+exp(-x))))
}


Nvars<-mySearch$Nvars.max
linx <-mySearch$Nvars.max+4
lHash<-length(hashStat)
mliks <- values(hashStat)[which((1:(lHash * linx)) %% linx == 1)]
betas <- values(hashStat)[which((1:(lHash * linx)) %% linx == 4)]
for(i in 1:(Nvars-1))
{
  betas<-cbind(betas,values(hashStat)[which((1:(lHash * linx)) %% linx == (4+i))])
}
betas<-cbind(betas,values(hashStat)[which((1:(lHash * linx)) %% linx == (0))])


# create 5 of na for all data columns in the test data

for(c in names(data.example1)[2:38])
{
  idna<-sample.int(n = 20720, size = 0.02*20720,replace = F)
  data.example1[[c]][idna]<-NA
}

is.na<-matrix(0,nrow = dim(data.example1)[1],ncol =  dim(data.example1)[2])
is.na[which(is.na(data.example1))]<-1
length(which(rowSums(is.na)>0))

system.time({

  res<-mySearch$forecast.matrix.na(link.g = g,covariates = (data.example1[1:20720,-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)]),betas = betas,mliks.in = mliks)$forecast

})

summary(res)

length(res)
res<-as.integer(res>=0.5)
length(which(res>=0.5))
length(which(res<0.5))
length(res)
length(which(data.example1$neo==1))

(1-sum(abs(res-data.example1$neo),na.rm = T)/20720)*100


#FNR
ps<-which(data.example1$neo==1)
sum(abs(res[ps]-data.example1$neo[ps]))/(sum(abs(res[ps]-data.example1$neo[ps]))+length(ps))*100

#FPR
ns<-which(data.example1$neo==0)
sum(abs(res[ns]-data.example1$neo[ns]))/(sum(abs(res[ns]-data.example1$neo[ns]))+length(ns))*100


classify <- read.table("/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/classify/unlabled.txt",fill=TRUE)

names(classify)<-names(data.example)[2:25]
View(head(classify))

transform<-colnames(data.example)[-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25)][1:9]
for(i in 1:length(transform))
{
  print(i)
  classify[[transform[i]]]<-as.numeric(as.character(classify[[transform[i]]]))
}

classify$esuar<-classify$eccentricity^2
classify$asuar<-classify$absolute_magnitude^2
classify$rsuar<-classify$semi_major_axis^2
classify$rcube<-classify$semi_major_axis^3
classify$anoml<-classify$mean_anomaly*classify$semi_major_axis
classify$anoms<-classify$mean_anomaly*classify$semi_major_axis^2
classify$anomv<-classify$mean_anomaly*classify$semi_major_axis^3
classify$perihell<-classify$argument_of_perihelion*classify$semi_major_axis
classify$perihels<-classify$argument_of_perihelion*classify$semi_major_axis^2
classify$perihelv<-classify$argument_of_perihelion*classify$semi_major_axis^3
classify$longitudel<-classify$longitude_of_the_ascending.node*classify$semi_major_axis
classify$longitudes<-classify$longitude_of_the_ascending.node*classify$semi_major_axis^2
classify$longitudev<-classify$longitude_of_the_ascending.node*classify$semi_major_axis^3

View(head(classify))



system.time({
  res<-NULL
  for(i in 1:24)
  {
    print(i)
    res<-c(res,as.array(mySearch$forecast.matrix(link.g = g,covariates = classify[((i-1)*50000+1):(i*50000),])$forecast))
  }
  
})
res<-c(res,as.array(mySearch$forecast.matrix(link.g = g,covariates = classify[((i)*50000+1):(1201528),])$forecast))


length(res)
length(which(res>0.5))

clasi<-NULL
clasi<-c(clasi,length(which(res>0.01))/length(res)*100)
         clasi<-c(clasi,length(which(res>0.10))/length(res)*100)
                  clasi<-c(clasi,length(which(res>0.20))/length(res)*100)
                           clasi<-c(clasi,length(which(res>0.30))/length(res)*100)
                                    clasi<-c(clasi,length(which(res>0.40))/length(res)*100)
                                             clasi<-c(clasi,length(which(res>0.50))/length(res)*100)
                                                      clasi<-c(clasi,length(which(res>0.60))/length(res)*100)
                                                               clasi<-c(clasi,length(which(res>0.70))/length(res)*100)
                                                                        clasi<-c(clasi,length(which(res>0.80))/length(res)*100)
                                                                                 clasi<-c(clasi,length(which(res>0.90))/length(res)*100)
                                                                                          clasi<-c(clasi,length(which(res>0.99))/length(res)*100)
print(clasi)
                                                                                          
                                                                                          
length(which(data.example1$neo>0.5))/length(data.example1$neo)
