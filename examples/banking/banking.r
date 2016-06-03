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
simx <- read.table("/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/banking/test (3).csv",sep = ",",header = T,fill=TRUE)

data.example <- simx

transform<-colnames(data.example)[-1]

for(i in 1:length(transform))
{
  print(i)
  data.example[[transform[i]]]<-as.numeric(as.character(data.example[[transform[i]]]))
}

data.example$logincome<-log(data.example$MonthlyIncome+1)
data.example$logincome2<-data.example$logincome^2
data.example$RevolvingUtilizationOfUnsecuredLines2<-data.example$RevolvingUtilizationOfUnsecuredLines^2
data.example$DebtRatio2<-data.example$DebtRatio^2
data.example$RevDebt<-data.example$DebtRatio*data.example$RevolvingUtilizationOfUnsecuredLines


data.example$missincome<-0
data.example$missdepend<-0

data.example$missincome[which(is.na(data.example$MonthlyIncome))]<-1
data.example$missdepend[which(is.na(data.example$NumberOfDependents))]<-1

data.example$MonthlyIncome[which(is.na(data.example$MonthlyIncome))]<-0
data.example$logincome[which(is.na(data.example$logincome))]<-0
data.example$logincome2[which(is.na(data.example$logincome2))]<-0
data.example$NumberOfDependents[which(is.na(data.example$NumberOfDependents))]<-0


data.example$incomrat<-data.example$RevolvingUtilizationOfUnsecuredLines*data.example$MonthlyIncome
data.example$incomrat2<-data.example$DebtRatio*data.example$MonthlyIncome
data.example$incomrat3<-data.example$DebtRatio*data.example$logincome2


data.example1<-data.example


simy<-read.table("/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/banking/train (4).csv",sep = ",",header = T,fill=TRUE)
data.example <- simy

transform<-colnames(data.example)[-1]

for(i in 1:length(transform))
{
  print(i)
  data.example[[transform[i]]]<-as.numeric(as.character(data.example[[transform[i]]]))
}

data.example$logincome<-log(data.example$MonthlyIncome+1)
data.example$logincome2<-data.example$logincome^2
data.example$RevolvingUtilizationOfUnsecuredLines2<-data.example$RevolvingUtilizationOfUnsecuredLines^2
data.example$DebtRatio2<-data.example$DebtRatio^2
data.example$RevDebt<-data.example$DebtRatio*data.example$RevolvingUtilizationOfUnsecuredLines


data.example$missincome<-0
data.example$missdepend<-0

data.example$missincome[which(is.na(data.example$MonthlyIncome))]<-1
data.example$missdepend[which(is.na(data.example$NumberOfDependents))]<-1

data.example$MonthlyIncome[which(is.na(data.example$MonthlyIncome))]<-0
data.example$logincome[which(is.na(data.example$logincome))]<-0
data.example$logincome2[which(is.na(data.example$logincome2))]<-0
data.example$NumberOfDependents[which(is.na(data.example$NumberOfDependents))]<-0


data.example$incomrat<-data.example$RevolvingUtilizationOfUnsecuredLines*data.example$MonthlyIncome
data.example$incomrat2<-data.example$DebtRatio*data.example$MonthlyIncome
data.example$incomrat3<-data.example$DebtRatio*data.example$logincome2

data.example[,1]<-data.example[,2]

names(data.example)[1]<-"Y"



#define the covariates and theobservations
fparam.example <- colnames(data.example)[-c(1,2)]
fobserved.example <- colnames(data.example)[1]
# view the correlation structure

estimate.bas.glm(Y~1,data = data.example,family = binomial(),prior = g.prior(g = 1000),logn=200)

#create MySearch object with default parameters
mySearch = EMJMCMC2016()
# load functions for MLIK estimation
mySearch$estimator = estimate.bas.glm
mySearch$estimator.args = list(data = data.example,prior = g.prior(g = 105000),family = binomial(), logn = log(105000))
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
resm<-mySearch$modejumping_mcmc(list(varcur=initsol,statid=5, distrib_of_proposals =distrib_of_proposals,distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = 0.000000001, trit = 999000, trest = 2500, burnin = 3, max.time = 30, maxit = 100000, print.freq =10))

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
workdir=""
# n/b visualization iisues on Windows! To be fixed!
mySearch$visualize_results(statistics1, "test", 11, crit=list(mlik = T, waic = T, dic = T),draw_dist =T)


# check how many models were visited
idn<-which(!is.na(statistics1[,1]))
length(idn)

# filter some models out (if needed!!!!), improves prediction speed but can also induce additional errors
statistics1[which(statistics1[,15]<0.001),]<-NA

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
which(statistics1[,1]==max(statistics1[,1],na.rm = T))
statistics1[,15]<-NA
statistics1[425915,15]<-1
length(idn)


# carry classification out and estimate the errors

g<-function(x)
{
  return((x = 1/(1+exp(-x))))
}
# plug in missing values for income



system.time({res<-mySearch$foreast.matrix(nvars = 20,ncases = 45000,link.g = g,covariates = data.example1[1:45000,-c(1,22)])$forecast})


length(which(res>=0.5))
length(which(res<0.5))
length(res)

xxx<-cbind(data.example1[,1],t(res))
write.csv(xxx, file = "MyData.csv",row.names = F,col.names = F)

