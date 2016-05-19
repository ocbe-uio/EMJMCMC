rm(list = ls(all = TRUE))
# install the required packges if needed
#install.packages("INLA", repos="http://www.math.ntnu.no/inla/R/testing")
#install.packages("bigmemory")
#install.packages("snow")
#install.packages("Rmpi")
#install.packages("ade4")
#install.packages("sp")
#install.packages("BAS")
#install.packages("/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC_1.2.tar.gz", repos = NULL, type="source")
#install.packages("RCurl")

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





# set up the parameters of the simulation or optimization

#define your working directory, where the data files are stored
workdir<-"/results"

#prepare data
simx <- read.table(paste(workdir,"proteincen.txt",sep = "",collapse = ""),sep = " ")
data.example <- as.data.frame(simx)
names(data.example)[89]="Y"


#fparam <- c("Const",colnames(data)[-1])
fparam.example <- colnames(data.example)[-89]
fobserved.example <- colnames(data.example)[89]

for(i in 1:length(fparam.example))
{
  fparam.example[i]=paste("I(V",i,")",sep = "")
}

# create either a standard hash table (default for now)
hashStat <- hash()

# or the one based on bigmemory package

#dataframe for results; n/b +1 is required for the summary statistics

#statistics1 <- big.matrix(nrow = 2 ^(23)+1, ncol = 16,init = NA, type = "double")
#statistics <- describe(statistics1)

#dataframe for results; n/b +1 is required for the summary statistics
#hash.keys1 <- big.matrix(nrow = 2 ^(23)+1, ncol = 88,init = 0, type = "char")
#hash.keys <- describe(hash.keys1)


remove(statistics1)
remove(statistics)
remove(hash.keys1)
remove(hash.keys)

#create MySearch object with default parameters
mySearch = EMJMCMC2016()

#statistics1 <- array(data = NA,dim = c(2 ^(20)+1,16))
#hash.keys1 <- array(data = NA,dim = c(2 ^(20)+1,88))


# load functions as in BAS article by Clyde, Ghosh and Littman to reproduce their first example
mySearch$estimator = estimate.bas.lm
mySearch$estimator.args = list(data = data.example,prior = 3, g = 96 ,n=96)
mySearch$parallelize = lapply

# this one MUST be completed before moving to the experiments
system.time(
  FFF<-mySearch$full_selection(list(statid=6, totalit =32769, ub = 13600,mlikcur=-Inf,waiccur =100000))
)
# check that all models are enumerated during the full search procedure
idn<-which(!is.na(statistics1[,1]))
length(idn)
hashStat

length(hashStat)

mySearch$Nvars

# see the best current results and total number of iterations
mySearch$g.results[,]
#compare to
which(statistics1[,1]==max(statistics1[,1],na.rm = TRUE))

View(statistics1[2113,])

which(statistics1[868,1]<= -10)

# get graphical output (only makes sence after EMJMCMC procedure carried out)
mySearch$visualize_results(statistics1 = statistics1, template = "test/full1BASclasstest", mds_size = 1024, crit = list(waic=TRUE,mlik=TRUE,dic=TRUE), draw_dist = TRUE)


# once full search is completed, get the truth for the experiment
ppp<-mySearch$post_proceed_results(statistics1)
truth = ppp$p.post # make sure it is equal to Truth column from the article
truth.m = ppp$m.post
truth.prob = ppp$s.mass
ordering = sort(ppp$p.post,index.return=T)
fake500 <- sum(exp(x = sort(statistics1[,1],decreasing = T)[1:5000]),na.rm = TRUE)/truth.prob

print("pi truth")
sprintf("%.10f",truth[ordering$ix])

plot(truth,type = "h")

mySearch$printable.opt=F
mySearch$max.cpu = as.integer(10)
mySearch$locstop.nd=FALSE
mySearch$max.cpu.glob = as.integer(10)
mySearch$max.N.glob=as.integer(12)
mySearch$min.N.glob=as.integer(6)
mySearch$max.N=as.integer(4)
mySearch$min.N=as.integer(1)
mySearch$recalc.margin = (50000)
distrib_of_proposals = c(76.91870,71.25264,1.68184,90.55921,17812.39852)
distrib_of_neighbourhoods=t(array(data = c(7.6651604,16.773326,14.541629,12.839445,12.964227,13.048343,7.165434,
                                           0.9936905,15.942490,11.040131,3.200394,15.349051,15.466632,4.676458,
                                           1.5184551,9.285762,6.125034,3.627547,13.343413,12.923767,5.318774,
                                           14.5295380,1.521960,11.804457,5.070282,6.934380,10.578945,2.455602,
                                           26.0826035,2.453729,14.340435,14.863495,10.028312,12.685017,53.806295),dim = c(7,5)))
mySearch$hash.length<-as.integer(23)
mySearch$double.hashing<-F

#distrib_of_neighbourhoods[,7]<-1

Niter <- 1
thining<-1
system.time({

  hashlist <- list()
  vect <-array(data = 0,dim = c(length(fparam.example),Niter))
  vect.mc <-array(data = 0,dim = c(length(fparam.example),Niter))
  inits <-array(data = 0,dim = Niter)
  freqs <-array(data = 100,dim = c(5,Niter))
  freqs.p <-array(data = 100,dim = c(5,7,Niter))
  masses <- array(data = 0,dim = Niter)
  iterats <- array(data = 0,dim = c(2,Niter))

  for(i in 1:Niter)
  {
    #statistics1 <- big.matrix(nrow = 2 ^(length(fparam.example))+1, ncol = 15,init = NA, type = "double")
    #statistics <- describe(statistics1)


    #hashStat <- hash()

    mySearch$g.results[4,1]<-0
    mySearch$g.results[4,2]<-0
    mySearch$p.add = array(data = 0.5,dim = length(fparam.example))
    #distrib_of_neighbourhoods=array(data = runif(n = 5*7,min = 0, max = 20),dim = c(5,7))
    #distrib_of_proposals = runif(n = 5,min = 0, max = 100)
    #distrib_of_proposals[5]=sum(distrib_of_proposals[1:4])*runif(n = 1,min = 50, max = 150)
    print("BEGIN ITERATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print(i)
    set.seed(10*i)
    initsol=rbinom(n =  length(fparam.example),size = 1,prob = 0.5)
    inits[i] <- mySearch$bittodec(initsol)
    freqs[,i]<- distrib_of_proposals
    resm<-mySearch$modejumping_mcmc(list(varcur=NULL,statid=-1, distrib_of_proposals =distrib_of_proposals,distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = 0.000000000001, trit = 2^30, trest = 2^20, burnin = 1000, max.time = 24*60*6, maxit = 2^20, print.freq = 1000
    ))
    vect[,i]<-resm$bayes.results$p.post
    vect.mc[,i]<-resm$p.post
    masses[i]<-resm$bayes.results$s.mass
    print(masses[i])
    freqs.p[,,i] <- distrib_of_neighbourhoods
    iterats[1,i]<-mySearch$g.results[4,1]
    iterats[2,i]<-mySearch$g.results[4,2]
    print("COMPLETE ITERATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! with")
    print(iterats[2,i])
    #clear(hashStat)
    #remove(hashStat)
    #remove(statistics1)
    #remove(statistics)
  }
}
)


mySearch$bittodec(c(0,1,0,1,0,1,0,0,0,1,0,1,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1,1,1,1,1,0,1,1,0,0,1,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,1,0,0,1,0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,0,1,0,0,1,1,0)
)

mySearch$g.results[1,1]<--1000
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

# #
View((cbind(ordering$ix/100,truth[ordering$ix]/100,bias.avg.rm[ordering$ix],rmse.avg.rm[ordering$ix],bias.avg.mc[ordering$ix],rmse.avg.mc[ordering$ix])*100))

which(hashStat == 77.3217980575684)

keys<-NULL
for(i in 2^30:(2^30+50000))
{
  bit <- mySearch$dectobit.alt(i+ as.integer(rnorm(n = 1,mean = 1000,sd = 100)))
  print(bit)
  print(mySearch$bittodec(bit))
}

plot(which(statistics1[,16]>=0))
