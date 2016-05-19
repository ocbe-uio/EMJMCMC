# R CMD build EMJMCMC2016 to build a package 


install.packages("/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC_1.2.tar.gz", repos = NULL, type="source")


library(package = EMJMCMC2016)
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
library(EMJMCMC)


data(data)
dim(data)


# set up the parameters of the simulation or optimization
M<-5
size<-1
# set up the parameters of the simulation or optimization
fparam.example <- colnames(data)[c(8:10,12:17,21:24,29)]
fobserved.example <- colnames(data)[5]
data.example <- data[c(1:500,7000:7500,10000:10500),]
data.example$express<-data.example$express+rnorm(n = 1502,mean = 1000,sd = 800)
#create MySearch object with default parameters
  mySearch = EMJMCMC2016()
mySearch$parallelize = lapply
#dataframe for results; n/b +1 is required for the summary statistics
statistics1 <- read.big.matrix(filename = "/mn/anatu/ansatte-u3/aliaksah/importantresults.csv")
statistics <- describe(statistics1)

mySearch$estimator = inla
args<-list(family = "poisson",data = data.example)
args$control.compute = list(dic = TRUE, waic = TRUE, mlik = TRUE)

formula2 <- as.formula("methylated_bases ~  1 + CHG + DT1 +f(data.example$pos,model=\"ar1\")")
fm4<-do.call(inla, c(args,formula = formula2))
summary(fm4)

mySearch$latent.formula  = ""; "+f(data.example$pos,model=\"ar1\")"
mySearch$estimator.args = args
mySearch$printable.opt = F


crits<-read.big.matrix(filename = "/mn/anatu/ansatte-u3/aliaksah/importantresults.csv")[,1:3]


# estimator based on precalculated and saved into crit data
esimator<-function(formula, crits)
{
  values <- strsplit(as.character(formula)[3],split = " + ",fixed = T)
  vec<-array(0,dim = Nvars)
  for(i in 2:(length(values[[1]])-1))
  {
    iid <- which(fparam.example == values[[1]][i])
    if(length(iid)>0)
      vec[iid]<-1
  }
  id<-bittodec(vec)+1
  return(list(mlik = crits[id,1],waic = crits[id,2] , dic =  crits[id,3]))
}

esimator(formula = formula2, crits = crits)
mySearch$estimator = esimator
mySearch$estimator.args = list(crits = crits)

# this one MUST be completed before moving to the experiments
system.time(
  FFF<-mySearch$full_selection(list(statid=-1, totalit =2^14+1, ub = 10,mlikcur=-Inf,waiccur =100000))
)
# completed in   7889  for 1048576 models whilst BAS took 6954.101 seonds and thus now advantage of using C versus R is clearly seen as neglectible  (14688.209 user seconds)
# BAS completed the same job in 

# check that all models are enumerated during the full search procedure
idn<-which(!is.na(statistics1[,1]))

statistics1[length(idn),]<-NA


View(idn)
statistics1[idn,1]<- -100000
idn<-which((statistics1[,1]==-99000))
length(idn)
statistics1[idn,4:15]<- 0
statistics1[,1]<-statistics1[,1]+1000 #multiple marginal log likelihoods with some constant
statistics1[,7]<-0

mySearch$visualize_results(statistics1, "test3",1024, crit=list(mlik = T, waic = T, dic = T),draw_dist = TRUE)

# once full search is completed, get the truth for the experiment
ppp<-mySearch$post_proceed_results(statistics1 = statistics1)
truth = ppp$p.post # make sure it is equal to Truth column from the article
truth.m = ppp$m.post
truth.prob = ppp$s.mass
ordering = sort(ppp$p.post,index.return=T)
fake500 <- sum(exp(x = (sort(statistics1[,1],decreasing = T)[1:2])),na.rm = TRUE)/truth.prob
print("pi truth")
sprintf("%.10f",truth[ordering$ix])

iddx <- sort(statistics1[,1],decreasing = T,index.return=T,na.last = NA)$ix
# check that all models are enumerated during the full search procedure

statistics1[as.numeric(iddx[376:2^14]),1:15]<-NA

# once full search is completed, get the truth for the experiment
ppp.best<-mySearch$post_proceed_results(statistics1 = statistics1)
best = ppp.best$p.post # make sure it is equal to Truth column from the article
bset.m = ppp.best$m.post
best.prob = ppp.best$s.mass/truth.prob
print("pi best")
sprintf("%.10f",best[ordering$ix])
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
# # 
View((cbind(best.bias[ordering$ix],best.rmse[ordering$ix])*100))


statistics1[1:2^14,5:14]<-0

mySearch$save_results_csv(statistics1, "important results")


mySearch$printable.opt=T
mySearch$max.cpu = as.integer(1)
mySearch$locstop.nd=FALSE
mySearch$max.cpu.glob = as.integer(1)
mySearch$max.N.glob=as.integer(1)
mySearch$min.N.glob=as.integer(1)
mySearch$max.N=as.integer(1)
mySearch$min.N=as.integer(1)
mySearch$recalc.margin = as.integer(2^13)


distrib_of_proposals = c(10,0,0,0,1000)

distrib_of_neighbourhoods=t(array(data = c(7.6651604,16.773326,14.541629,12.839445,2.964227,13.048343,7.165434,
                                           0.9936905,15.942490,11.040131,3.200394,15.349051,5.466632,14.676458,
                                           1.5184551,9.285762,6.125034,3.627547,13.343413,2.923767,15.318774,
                                           14.5295380,1.521960,11.804457,5.070282,6.934380,10.578945,12.455602,
                                           1,1,1,1,1,1,1),dim = c(7,5)))

Niter <- 1
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
    #distrib_of_neighbourhoods=array(data = runif(n = 5*7,min = 0, max = 20),dim = c(5,7))
    #distrib_of_proposals = runif(n = 5,min = 0, max = 100)
    #distrib_of_proposals[5]=sum(distrib_of_proposals[1:4])*runif(n = 1,min = 50, max = 150)
    print("BEGIN ITERATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print(i)
    set.seed(10*i)
    initsol=rbinom(n = length(fparam.example),size = 1,prob = 0.5)
    inits[i] <- mySearch$bittodec(initsol)
    freqs[,i]<- distrib_of_proposals
    resm<-mySearch$modejumping_mcmc(list(varcur=initsol,statid=-1, distrib_of_proposals =distrib_of_proposals,distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = 0.000000001, trit = 50000, trest = 3000, burnin = 3, max.time = 30, maxit = 100000, print.freq =500))
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

# load functions as in BAS article by Clyde, Ghosh and Littman to reproduce their first example
mySearch$estimator = estimate.bas.lm
mySearch$estimator.args = list(data = data.example,prior = 3, g = 96 ,n=96)

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
