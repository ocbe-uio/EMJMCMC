rm(list = ls(all = TRUE))

#install the required packges if needed
install.packages("INLA", repos="http://www.math.ntnu.no/inla/R/testing")
install.packages("bigmemory")
install.packages("snow")
install.packages("Rmpi")
install.packages("ade4")
install.packages("sp")
install.packages("BAS")
#install.packages("https://github.com/aliaksah/EMJMCMC2016/files/270429/EMJMCMC_1.3.tar.gz", repos = NULL, type="source")
install.packages("RCurl")
install.packages("hash")
install.packages("copula")

library(hash)
library(RCurl)
#library(EMJMCMC) bulid the sourse file https://github.com/aliaksah/EMJMCMC2016/blob/master/R/the_mode_jumping_package2.r instead, since the latest published build EMJMCMC_1.2.tar.gz doen not have prediction functionality (coming in the package soon with other cookies ;) ) available
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
# view the correlation structure

View(cor(data.example[,-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25)],use = "complete.obs"))


#create MySearch object with default parameters
mySearch = EMJMCMC2016()
# load functions for MLIK estimation
mySearch$estimator = estimate.bas.glm
mySearch$estimator.args = list(data = data.example,prior = aic.prior(),family = binomial(), logn = log(64))
mySearch$save.beta=T



# set parameters of the search and stuff
mySearch$max.cpu = as.integer(4)
mySearch$locstop.nd=FALSE
mySearch$max.cpu.glob = as.integer(4)
mySearch$switch.type=as.integer(1)
mySearch$switch.type.glob=as.integer(1)
#mySearch$printable.opt = TRUE
mySearch$max.N.glob=as.integer(10)
mySearch$min.N.glob=as.integer(5)
mySearch$max.N=as.integer(2)
mySearch$min.N=as.integer(1)
mySearch$recalc.margin = as.integer(2^20)
distrib_of_proposals = c(76.91870,71.25264,87.68184,60.55921,15812.39852)
#distrib_of_proposals = с(0,0,0,0,10)
distrib_of_neighbourhoods=t(array(data = c(7.6651604,16.773326,14.541629,12.839445,2.964227,13.048343,7.165434,
                                           0.9936905,15.942490,11.040131,3.200394,15.349051,5.466632,14.676458,
                                           1.5184551,9.285762,6.125034,3.627547,13.343413,2.923767,15.318774,
                                           14.5295380,1.521960,11.804457,5.070282,6.934380,10.578945,12.455602,
                                           6.0826035,2.453729,14.340435,14.863495,1.028312,12.685017,13.806295),dim = c(7,5)))
#carry the search (training out)
statistics1 <- big.matrix(nrow = 2 ^(length(fparam.example))+1, ncol = 16+length(fparam.example),init = NA, type = "double")
statistics <- describe(statistics1)
mySearch$g.results[4,1]<-0
mySearch$g.results[4,2]<-0
mySearch$p.add = array(data = 0.5,dim = length(fparam.example))
initsol=rbinom(n = length(fparam.example),size = 1,prob = 0.5)
resm<-mySearch$modejumping_mcmc(list(varcur=initsol,statid=5, distrib_of_proposals =distrib_of_proposals,distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = 0.000000001, trit = 999000, trest = 20000, burnin = 3, max.time = 30, maxit = 100000, print.freq =50))

# save the results (if needed)
write.big.matrix(x = statistics1,filename = "cosmoneo.csv")


#post proceed data and visualize it (if needed)

ppp<-mySearch$post_proceed_results(statistics1 = statistics1)
truth = ppp$p.post # make sure it is equal to Truth column from the article
truth.m = ppp$m.post
truth.prob = ppp$s.mass
ordering = sort(ppp$p.post,index.return=T)
print("pi truth")
sprintf("%.10f",truth[ordering$ix])
sprintf(fparam.example[ordering$ix])

par(mar = c(10,4,4,2) + 4.1)
barplot(resm$bayes.results$p.post,density = 46,border="black",main = "Marginal Inclusion (RM)",ylab="Probability",names.arg =fparam.example,las=2)
barplot(resm$p.post,density = 46,border="black",main = "Marginal Inclusion (MC)",ylab="Probability",names.arg =fparam.example,las=2)

template = "test"
# n/b visualization iisues on Windows! To be fixed!
#mySearch$visualize_results(statistics1, "test", 2^10-1, crit=list(mlik = T, waic = T, dic = T),draw_dist =F)


# check how many models were visited
idn<-which(!is.na(statistics1[,1]))
length(idn)

# filter some models out (if needed!!!!), improves prediction speed but can also induce additional errors
statistics1[which(statistics1[,15]<0.0001),]<-NA

# check how many models were visited
idn<-which(!is.na(statistics1[,1]))
length(idn)

# check the filtering effect (if required)
ppp<-mySearch$post_proceed_results(statistics1 = statistics1)
truth = ppp$p.post # make sure it is equal to Truth column from the article
truth.m = ppp$m.post
truth.prob = ppp$s.mass
ordering = sort(ppp$p.post,index.return=T)
print("pi truth")
sprintf("%.10f",truth[ordering$ix])
sprintf(fparam.example[ordering$ix])

par(mar = c(10,4,4,2) + 4.1)
barplot(resm$bayes.results$p.post,density = 46,border="black",main = "Marginal Inclusion (RM)",ylab="Probability",names.arg =fparam.example,las=2)
barplot(resm$p.post,density = 46,border="black",main = "Marginal Inclusion (MC)",ylab="Probability",names.arg =fparam.example,las=2)

# check how many models were visited
idn<-which(is.na(statistics1[,1]))
length(idn)



# carry classification out and estimate the errors

g<-function(x)
{
  return((x = 1/(1+exp(-x))))
}

system.time({
tot.er<-0
fp<-0
fn<-0
test.size<-20720
ids<-sample.int(size = test.size,n = 20720,replace = F)
ns<-0
ps<-0
for(i in ids)
{
  X<-c(1,as.numeric(data.example1[i,-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)]))
  y.hat<-round(mySearch$foreast(covariates = X,nvars = 36,link.g = g)$forecast)
  error<-as.numeric(data.example1[i,1])-y.hat
  if(error == 1)
    fn<-fn+1
  else
    if(error == -1)
      fp<-fp+1
  if(as.numeric(data.example1[i,1])==1)
    ps<-ps+1
  else
    ns<-ns+1
  print(paste(i," ",y.hat," ",(as.numeric(data.example1[i,1])-y.hat)))
}
tot.er<-(fp+fn)
# error
print(100*tot.er/test.size)
# precision
print((1-(tot.er/test.size))*100)
# false positive rate
print(fp*100/(ns+fp))
# false negative rate
print(fn*100/(ps+fn))
})

data.example1[1:20720,-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)][(is.na(data.example1[1:20720,-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)]))]<-0
#vectorized classification (much faster)
res<-mySearch$foreast.matrix(nvars = 20,ncases = 20720,link.g = g,covariates = data.example1[1:20720,-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)])$forecast
res<-as.integer(res>=0.5)
length(which(res>=0.5))
length(which(res<0.5))
length(res)
length(which(data.example1$neo==1))

(1-sum(abs(res-data.example1$neo),na.rm = T)/20720)

