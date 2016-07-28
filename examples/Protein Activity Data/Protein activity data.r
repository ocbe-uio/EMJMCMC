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
#library(EMJMCMC)
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

# system.time({
#
# formula1 = as.formula(paste(colnames(data.example)[89],"~ 1 +",paste0(colnames(data.example)[-89],collapse = "+")))
#
# res = runemjmcmc(formula = formula1,data = data.example,estimator =estimate.bas.lm,estimator.args =  list(data = data.example,prior = 3, g = 96 ,n=96),save.beta = F,n.models = 20000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = T,create.hash = F,pseudo.paral = F,burn.in = 100,print.freq = 100,advanced.param = list(
#                                                                                                                                                                                                                                                                                                                max.N.glob=as.integer(20),
#                                                                                                                                                                                                                                                                                                                min.N.glob=as.integer(5),
#                                                                                                                                                                                                                                                                                                                max.N=as.integer(3),
#                                                                                                                                                                                                                                                                                                                min.N=as.integer(1),
#                                                                                                                                                                                                                                                                                                                printable = F))
# })
#
# View(statistics1[1:1000,])


# create either a standard hash table (default for now)
hashStat <- hash()

# or the one based on bigmemory package N/B do not create both data objects simultaneously, since this can lead to unpredicted results

#dataframe for results; n/b +1 is required for the summary statistics

#statistics1 <- big.matrix(nrow = 2 ^(23)+1, ncol = 16,init = NA, type = "double")
#statistics <- describe(statistics1)

#dataframe for results; n/b +1 is required for the summary statistics
#hash.keys1 <- big.matrix(nrow = 2 ^(23)+1, ncol = 88,init = 0, type = "char")
#hash.keys <- describe(hash.keys1)


#create MySearch object with default parameters
mySearch = EMJMCMC2016()


# load functions as in BAS article by Clyde, Ghosh and Littman to reproduce their first example
mySearch$estimator = estimate.bas.lm
mySearch$estimator.args = list(data = data.example,prior = 3, g = 96 ,n=96)
mySearch$parallelize = lapply# if the hash provided by Decision Patterns is used parallel computing is not performed!?
#
# full enumeration is infeasible
# system.time(
#   FFF<-mySearch$full_selection(list(statid=6, totalit =32769, ub = 13600,mlikcur=-Inf,waiccur =100000))
# )
# # check that all models are enumerated during the full search procedure
# idn<-which(!is.na(statistics1[,1]))
# length(idn)
# hashStat

# define parameters of the search

mySearch$printable.opt=F
mySearch$max.cpu = as.integer(10)
mySearch$locstop.nd=FALSE
mySearch$max.cpu.glob = as.integer(10)
mySearch$max.N.glob=as.integer(20)
mySearch$min.N.glob=as.integer(5)
mySearch$max.N=as.integer(3)
mySearch$min.N=as.integer(1)
mySearch$recalc.margin = (500000)
distrib_of_proposals = c(76.91870,71.25264,87.68184,90.55921,17812.39852)
distrib_of_neighbourhoods=t(array(data = c(7.6651604,16.773326,14.541629,12.839445,12.964227,13.048343,7.165434,
                                           0.9936905,15.942490,11.040131,3.200394,15.349051,15.466632,4.676458,
                                           1.5184551,9.285762,6.125034,3.627547,13.343413,12.923767,5.318774,
                                           14.5295380,1.521960,11.804457,5.070282,6.934380,10.578945,2.455602,
                                           26.0826035,12.453729,14.340435,14.863495,10.028312,12.685017,13.806295),dim = c(7,5)))
mySearch$hash.length<-as.integer(20)
mySearch$double.hashing<-T

#Proceed for the predefined number of iterations

Niter <- 10
thining<-1
system.time({

  mliklist<-array(data = 0, dim = c(2^mySearch$hash.length, Niter))
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


    hashStat <- hash()


    mySearch$g.results[1,1]<--Inf
    mySearch$g.results[1,2]<-1
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
    resm<-mySearch$modejumping_mcmc(list(varcur=NULL,statid=-1, distrib_of_proposals =distrib_of_proposals,distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = 0.000000000001, trit = 2^30, trest = 2^20, burnin = 100, max.time = 24*60*6, maxit = 2^20, print.freq = 1000
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
    lHash<-length(hashStat)
    mliks <- values(hashStat)[which((1:(lHash * 3)) %%3 == 1)]
    mliklist[,i]<-mliks[1:2^mySearch$hash.length]


    lHash<-length(hashStat)
    mliks <- values(hashStat)[which((1:(lHash * 3)) %%3 == 1)]
    sum(exp(mliks))
    smilks100000<-sort(mliks,decreasing = T)[1:100000]
    boxplot(smilks100000,xaxt="n",ylab="log(Marginal Likelihood)",xlab="Replicates",horizontal=FALSE,pch=".",cex.lab=1.7,cex.axis=1.5,omd=c(0,0.7,0,0.7))
    smilks100000[1:10]


    write(mliklist[,i], file = paste("mliks",i,".csv"),
          ncolumns = 1,
          append = FALSE, sep = " ")
    write(vect[,i], file = paste("pp",i,".rs.csv"),
          ncolumns = 1,
          append = FALSE, sep = " ")
    write(vect.mc[,i], file = paste("pp",i,".mc.csv"),
          ncolumns = 1,
          append = FALSE, sep = " ")

    remove(hashStat)
    #clear(hashStat)
    #remove(hashStat)
    #remove(statistics1)
    #remove(statistics)
  }
}
)


print("model coverages")
mean(masses)
median(masses)
print("mean # of iterations")# even smaller on average than in BAS
mean(iterats[1,])
print("mean # of estimations")# even smaller on average than in BAS
mean(iterats[2,])

