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

#define your working directory, where the data files are stored
workdir<-""

#prepare data
simx <- read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/Simulated%20Logistic%20Data%20With%20Multiple%20Modes%20%28Example%203%29/sim3-X.txt"),sep = ",")
simy <-  read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/Simulated%20Logistic%20Data%20With%20Multiple%20Modes%20%28Example%203%29/sim3-Y.txt"),sep = ",")
data.example <- cbind(simy,simx)
names(data.example)[1]="Y1"

data.example$V2<-(data.example$V10+data.example$V14)*data.example$V9
data.example$V5<-(data.example$V11+data.example$V15)*data.example$V12
#fparam <- c("Const",colnames(data)[-1])
fparam.example <- colnames(data.example)[-1]
fobserved.example <- colnames(data.example)[1]

#dataframe for results; n/b +1 is required for the summary statistics
statistics1 <- big.matrix(nrow = 2 ^(length(fparam.example))+1, ncol = 15,init = NA, type = "double")
statistics <- describe(statistics1)

#create MySearch object with default parameters
mySearch = EMJMCMC2016()
# load functions for MLIK estimation
mySearch$estimator = estimate.bas.glm
mySearch$estimator.args = list(data = data.example,prior = aic.prior(),family = binomial(), logn = log(2000))

#play around with various methods in order to get used to them and see how they work

# carry out full enumeration (it is still feasible)
system.time(
  FFF<-mySearch$full_selection(list(statid=6, totalit =2^20+1, ub = 36*20,mlikcur=-Inf,waiccur =100000))
)
# completed in   7889  for 1048576 models whilst BAS took 6954.101 seonds and thus now advantage of using C versus R is clearly seen as neglectible  (14688.209 user seconds)
# BAS completed the same job in

# check that all models are enumerated during the full search procedure
idn<-which(is.na(statistics1[,1]))
length(idn)

mySearch$visualize_results(statistics1, "test3", 1024, crit=list(mlik = T, waic = T, dic = T),draw_dist = FALSE)
# once full search is completed, get the truth for the experiment
ppp<-mySearch$post_proceed_results(statistics1 = statistics1)
truth = ppp$p.post # make sure it is equal to Truth column from the article
truth.m = ppp$m.post
truth.prob = ppp$s.mass
ordering = sort(ppp$p.post,index.return=T)
fake500 <- sum(exp(x = (sort(statistics1[,1],decreasing = T)[1:2^20] + 1)),na.rm = TRUE)/truth.prob
print("pi truth")
sprintf("%.10f",truth[ordering$ix])

#estimate best performance ever
min(statistics1[,1],na.rm = T)
idn<-which(is.na(statistics1[,1]))
2^20-length(idn)
statistics1[idn,1]<- -100000
iddx <- sort(statistics1[,1],decreasing = T,index.return=T,na.last = NA)$ix
# check that all models are enumerated during the full search procedure




# see the obtained maximum and minimum

min(statistics1[,1],na.rm = TRUE)
max(statistics1[,1],na.rm = TRUE)

# look at the best possible performance
statistics1[as.numeric(iddx[10001:2^20]),1:15]<-NA
ppp.best<-mySearch$post_proceed_results(statistics1 = statistics1)
best = ppp.best$p.post # make sure it is equal to Truth column from the article
bset.m = ppp.best$m.post
best.prob = ppp.best$s.mass/truth.prob
print("pi best")
sprintf("%.10f",best[ordering$ix])
# notice some interesting details on the posterior mass and number of models visited
# 50000 best models contain 100.0000% of mass 100.0000%
# 48300 best models contain 99.99995% of mass 100.0000%
# 48086 best models contain 99.99995% of mass 100.0000%
# 10000 best models contain 99.99990% of mass 99.99991%
# 5000  best models contain 93.83923% of mass 94.72895%
# 3500  best models contain 85.77979% of mass 87.90333%
# 1500  best models contain 63.33376% of mass 67.71380%
# 1000  best models contain 53.47534% of mass 57.91971%
# 500   best models contain 37.72771% of mass 42.62869%
# 100   best models contain 14.76030% of mass 17.71082%
# 50    best models contain 14.76030% of mass 11.36970%
# 10    best models contain 14.76030% of mass 3.911063%
# 5     best models contain 14.76030% of mass 2.351454%
# 1     best models contain 14.76030% of mass 0.595301%

best.bias.m<-sqrt(mean((bset.m - truth.m)^2,na.rm = TRUE))*100000
best.rmse.m<-sqrt(mean((bset.m - truth.m)^2,na.rm = TRUE))*100000

best.bias<- best - truth
best.rmse<- abs(best - truth)
# view results for the best possible performance model
View((cbind(best.bias[ordering$ix],best.rmse[ordering$ix])*100))

# proceed with the experiment

# set parameters of the search

mySearch$switch.type=as.integer(1)
mySearch$switch.type.glob=as.integer(1)
#mySearch$printable.opt = TRUE
mySearch$max.N.glob=as.integer(5)
mySearch$min.N.glob=as.integer(3)
mySearch$max.N=as.integer(1)
mySearch$min.N=as.integer(1)
mySearch$recalc.margin = as.integer(2^20)
distrib_of_proposals = c(76.91870,71.25264,87.68184,60.55921,15812.39852)
#distrib_of_proposals = с(0,0,0,0,10)
distrib_of_neighbourhoods=t(array(data = c(7.6651604,16.773326,14.541629,12.839445,2.964227,13.048343,7.165434,
                                           0.9936905,15.942490,11.040131,3.200394,15.349051,5.466632,14.676458,
                                           1.5184551,9.285762,6.125034,3.627547,13.343413,2.923767,15.318774,
                                           14.5295380,1.521960,11.804457,5.070282,6.934380,10.578945,12.455602,
                                           6.0826035,2.453729,14.340435,14.863495,1.028312,12.685017,13.806295),dim = c(7,5)))
distrib_of_neighbourhoods[7]=distrib_of_neighbourhoods[7]/100
distrib_of_neighbourhoods = array(data = 0, dim = c(5,7))
distrib_of_neighbourhoods[,3]=10
# proceed with the search
Niter <- 100
thining<-1
system.time({

  vect <-array(data = 0,dim = c(length(fparam.example),Niter))
  vect.mc <-array(data = 0,dim = c(length(fparam.example),Niter))
  inits <-array(data = 0,dim = Niter)
  freqs <-array(data = 100,dim = c(5,Niter))
  freqs.p <-array(data = 100,dim = c(5,7,Niter))
  masses <- array(data = 0,dim = Niter)
  biases.m <- array(data = 0,dim = 2 ^(length(fparam.example))+1)
  biases.m.mc <- array(data = 0,dim = 2 ^(length(fparam.example))+1)
  rmse.m <- array(data = 0,dim = Niter)
  rmse.m.mc <- array(data = 0,dim = Niter)
  iterats <- array(data = 0,dim = c(2,Niter))

  for(i in 1:Niter)
  {
    statistics1 <- big.matrix(nrow = 2 ^(length(fparam.example))+1, ncol = 15,init = NA, type = "double")
    statistics <- describe(statistics1)
    mySearch$g.results[4,1]<-0
    mySearch$g.results[4,2]<-0
    mySearch$p.add = array(data = 0.5,dim = length(fparam.example))
    print("BEGIN ITERATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print(i)
    set.seed(10*i)
    initsol=rbinom(n = length(fparam.example),size = 1,prob = 0.5)
    inits[i] <- mySearch$bittodec(initsol)
    freqs[,i]<- distrib_of_proposals
    resm<-mySearch$modejumping_mcmc(list(varcur=initsol,statid=5, distrib_of_proposals =distrib_of_proposals,distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = 0.000000001, trit = 999000, trest = 10000, burnin = 3, max.time = 30, maxit = 100000, print.freq =1000))
    vect[,i]<-resm$bayes.results$p.post
    vect.mc[,i]<-resm$p.post
    masses[i]<-resm$bayes.results$s.mass/truth.prob
    print(masses[i])
    freqs.p[,,i] <- distrib_of_neighbourhoods
    cur.p.post <- resm$bayes.results$m.post
    cur.p.post[(which(is.na(cur.p.post)))]<-0
    rmse.m[i]<-mean((cur.p.post - truth.m)^2,na.rm = TRUE)
    biases.m<-biases.m + (cur.p.post - truth.m)
    cur.p.post.mc <- resm$m.post
    cur.p.post.mc[(which(is.na(cur.p.post.mc)))]<-0
    rmse.m.mc[i]<-mean((cur.p.post.mc - truth.m)^2,na.rm = TRUE)
    biases.m.mc<-biases.m.mc + (cur.p.post.mc - truth.m)
    iterats[1,i]<-mySearch$g.results[4,1]
    iterats[2,i]<-mySearch$g.results[4,2]
    print("COMPLETE ITERATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! with")
    print(iterats[2,i])
    remove(statistics1)
    remove(statistics)
  }
}
)



Nlim <- 1
order.deviat <- sort(masses,decreasing = TRUE,index.return=T)


print("model bias rm")
sqrt(mean((biases.m/Niter)^2,na.rm = TRUE))*100000
print("model rmse rm")
sqrt(mean(rmse.m))*100000

print("model bias mc")
sqrt(mean((biases.m.mc/Niter)^2,na.rm = TRUE))*100000
print("model rmse mc")
sqrt(mean(rmse.m.mc))*100000


print("model coverages")
mean(masses)
median(masses)
print("mean # of iterations")# even smaller on average than in BAS
mean(iterats[1,])
print("mean # of estimations")# even smaller on average than in BAS
mean(iterats[2,])

hist(masses)


# correlation between the MSE and the masses, obviously almost minus 1
cor(rmse.m,masses)
cor(rmse.m.mc,masses)
cor(iterats[2,],masses)

truth.buf <- array(data = 0,dim = c(length(fparam.example),Niter))
truth.buf[,1:Niter]<-truth
bias <- vect - truth.buf
bias.mc <- vect.mc - truth.buf
rmse <- (vect^2 +truth.buf^2 - 2*vect*truth.buf)
rmse.mc <- (vect.mc^2 +truth.buf^2 - 2*vect.mc*truth.buf)
bias.avg.rm<-rowMeans(bias)
rmse.avg.rm <-sqrt(rowMeans(rmse))
bias.avg.mc<-rowMeans(bias.mc)
rmse.avg.mc <-sqrt(rowMeans(rmse.mc))


print("pi biases rm")
sprintf("%.10f",bias.avg.rm[ordering$ix]*100)
print("pi rmse rm")
sprintf("%.10f",rmse.avg.rm[ordering$ix]*100)

print("pi biases mc")
sprintf("%.10f",bias.avg.mc[ordering$ix]*100)
print("pi rmse mc")
sprintf("%.10f",rmse.avg.mc[ordering$ix]*100)

# view the final results
View((cbind(bias.avg.rm[ordering$ix],rmse.avg.rm[ordering$ix],bias.avg.mc[ordering$ix],rmse.avg.mc[ordering$ix])*100))
