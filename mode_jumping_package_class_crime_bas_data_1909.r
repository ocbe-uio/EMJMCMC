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

#define your working directory, where the data files are stored
workdir<-"/results"

#prepare data
simx <- read.table(text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/simcen-x1.txt"),sep = ",")
simy <- read.table(text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/simcen-y1.txt"))
data.example <- cbind(simy,simx)
names(data.example)[1]="Y"

#fparam <- c("Const",colnames(data)[-1])
fparam.example <- colnames(data.example)[-1]
fobserved.example <- colnames(data.example)[1]

#dataframe for results; n/b +1 is required for the summary statistics
statistics1 <- big.matrix(nrow = 2 ^(length(fparam.example))+1, ncol = 15,init = NA, type = "double")
statistics <- describe(statistics1)

#create MySearch object with default parameters
mySearch = EMJMCMC2016()


# load functions as in BAS article by Clyde, Ghosh and Littman to reproduce their first example
mySearch$estimator = estimate.bas.lm
mySearch$estimator.args = list(data = data.example,prior = 3, g = 47 ,n=47)


# carry out full enumeration to discover the truth
system.time(
  FFF<-mySearch$full_selection(list(statid=6, totalit =32769, ub = 36,mlikcur=-Inf,waiccur =100000))
)
# check that all models are enumerated during the full search procedure
idn<-which(!is.na(statistics1[,1]))
length(idn)

# once full search is completed, get the truth for the experiment
ppp<-mySearch$post_proceed_results(statistics1 = statistics1)
truth = ppp$p.post # make sure it is equal to Truth column from the article
truth.m = ppp$m.post
truth.prob = ppp$s.mass
ordering = sort(ppp$p.post,index.return=T)
fake500 <- sum(exp(x = sort(statistics1[,1],decreasing = T)[1:3276]),na.rm = TRUE)/truth.prob
print("pi truth")
sprintf("%.10f",truth[ordering$ix])

#estimate best performance ever
iddx <- sort(statistics1[,1],decreasing = T,index.return=T)$ix
statistics1[as.numeric(iddx[2000:2^15]),1:15]<-NA

# once full search is completed, get the truth for the experiment
ppp.best<-mySearch$post_proceed_results(statistics1 = statistics1)
best = ppp.best$p.post # make sure it is equal to Truth column from the article
bset.m = ppp.best$m.post
best.prob = ppp.best$s.mass/truth.prob
print("pi best")
sprintf("%.10f",best[ordering$ix])


best.bias.m<-sqrt(mean((bset.m - truth.m)^2,na.rm = TRUE))*100000
best.rmse.m<-sqrt(mean((bset.m - truth.m)^2,na.rm = TRUE))*100000

best.bias<- best - truth
best.rmse<- abs(best - truth)
# #
View((cbind(best.bias[ordering$ix],best.rmse[ordering$ix])*100))



print("pi truth")
sprintf("%.10f",truth[ordering$ix])

# vecbuf<-vect
# freqs.buf<- freqs
# freqs.p.buf <- freqs.p

#reproduce the 1st experiment as in BAS article


mySearch$switch.type=as.integer(5)
mySearch$switch.type.glob=as.integer(1)
mySearch$max.N.glob=as.integer(6)
mySearch$min.N.glob=as.integer(5)
mySearch$max.N=as.integer(12)
mySearch$min.N=as.integer(5)
mySearch$recalc.margin = as.integer(2^15+1)
mySearch$max.cpu=as.integer(4)
mySearch$p.add = array(data = 0.5,dim = 15)
mySearch$printable.opt = T
mySearch$p.add = array(data = 0.5,dim = 15)
fff<-mySearch$forward_selection(list(varcur=rep(0,length(fparam.example)),mlikcur=-Inf,waiccur =Inf,locstop = FALSE,statid=-1))
bbb<-mySearch$backward_selection(list(varcur=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),mlikcur=-Inf,waiccur =Inf,locstop = FALSE,statid=-1))


cor(data.example)[1,]

resss<-0
masss<-0
mmm<-0

statistics1 <- big.matrix(nrow = 2 ^(length(fparam.example))+1, ncol = 15,init = NA, type = "double")
statistics <- describe(statistics1)
mySearch$g.results[4,1]<-0
mySearch$g.results[4,2]<-0
distrib_of_neighbourhoods=array(data = 2,dim = c(5,7))
distrib_of_neighbourhoods[,1]<-8
distrib_of_neighbourhoods[,5]<-0
distrib_of_neighbourhoods[,6]<-0

while(mmm<9900)
{

  initsol=rbinom(n = length(fparam.example),size = 1,prob = 0.5)
  mySearch$p.add = array(data = 0.5,dim = 15)
  resm<-mySearch$modejumping_mcmc(list(varcur=initsol,statid=5, distrib_of_proposals = c(0,0,50,0,5500),distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = 0.0001, trit = 3200, trest = 32000 , burnin = 30, max.time = 30, maxit = 100000, print.freq =5))
  resss<-c(resss,mySearch$g.results[4,1])
  mySearch$g.results[4,2]
  mySearch$g.results[4,1]
  resm$bayes.results$s.mass/truth.prob
  mmm<-mySearch$post_proceed_results(statistics1 = statistics1)$s.mass/truth.prob*10000
  masss<-c(masss,mmm)
  plot(resss,ylim = c(0,20000),type = "l", main = "EMJMCMC", ylab = "iterations, mass x 10^3",xlab = "test points (every 50 iterations of EMJCMCMC)" )
  lines(masss,col = 3)
  #   set.seed(length(resss)*20000)
  #   mySearch$seed=as.integer(runif(n = 1,min = length(resss),max = length(resss)*20000))
  #   mySearch$forw_backw_walk(list(steps=2,p1=0.5,reverse = FALSE))
  #   resss<-c(resss,mySearch$g.results[4,2])
  #   mmm<-mySearch$post_proceed_results(statistics1 = statistics1)$s.mass/truth.prob*10000
  #   masss<-c(masss,mmm)
  #   plot(resss,ylim = c(0,10000))
  #   points(masss,col = 4)
  #   set.seed(length(resss)*10000)
}

mySearch$switch.type=as.integer(1)
mySearch$switch.type.glob=as.integer(1)
mySearch$printable.opt = TRUE


mySearch$max.N.glob=as.integer(4)
mySearch$min.N.glob=as.integer(4)
mySearch$max.N=as.integer(2)
mySearch$min.N=as.integer(2)
mySearch$recalc.margin = as.integer(400)
distrib_of_proposals = c(76.91870,71.25264,87.68184,60.55921,17812.39852)
distrib_of_neighbourhoods=t(array(data = c(7.6651604,16.773326,14.541629,12.839445,2.964227,13.048343,7.165434,
                                           0.9936905,15.942490,11.040131,3.200394,15.349051,5.466632,14.676458,
                                           1.5184551,9.285762,6.125034,3.627547,13.343413,2.923767,15.318774,
                                           14.5295380,1.521960,11.804457,5.070282,6.934380,10.578945,12.455602,
                                           6.0826035,2.453729,14.340435,14.863495,1.028312,12.685017,13.806295),dim = c(7,5)))
distrib_of_neighbourhoods[7]=distrib_of_neighbourhoods[7]/50
Niter <- 100
thining<-1
system.time({

  vect <-array(data = 0,dim = c(15,Niter))
  vect.mc <-array(data = 0,dim = c(15,Niter))

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
    mySearch$p.add = array(data = 0.5,dim = 15)
    #distrib_of_neighbourhoods=array(data = runif(n = 5*7,min = 0, max = 20),dim = c(5,7))
    #distrib_of_proposals = runif(n = 5,min = 0, max = 100)
    #distrib_of_proposals[5]=sum(distrib_of_proposals[1:4])*runif(n = 1,min = 50, max = 150)
    print("BEGIN ITERATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print(i)
    set.seed(10*i)
    initsol=rbinom(n = length(fparam.example),size = 1,prob = 0.8)
    inits[i] <- mySearch$bittodec(initsol)
    freqs[,i]<- distrib_of_proposals
    resm<-mySearch$modejumping_mcmc(list(varcur=initsol,statid=5, distrib_of_proposals = distrib_of_proposals,distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = 0.0001, trit = 3270, trest = 32000 , burnin = 50, max.time = 30, maxit = 100000, print.freq =500))
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
print("mean # of iterations")# even smaller on average than in BAS
mean(iterats[1,])
print("mean # of estimations")# even smaller on average than in BAS
mean(iterats[2,])



# correlation between the MSE and the masses, obviously almost minus 1
cor(rmse.m,masses)
cor(rmse.m.mc,masses)

truth.buf <- array(data = 0,dim = c(15,Niter))
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
View((cbind(bias.avg.rm[ordering$ix],rmse.avg.rm[ordering$ix],bias.avg.mc[ordering$ix],rmse.avg.mc[ordering$ix])*100))


# indexs<-(sort.int(cols.b, index.return=TRUE))
# plot(y = inits[indexs$ix], x = 1:1000,col = 1, type = "p")
# points(cols.b[indexs$ix]*10000,col = 4)
# points(cols.r[indexs$ix]*50000,col = 5)
# points(cols.r[indexs$ix]*50000+cols.b[indexs$ix]*10000,col = 3)
# points(freqs.p[1,indexs$ix]*30,col= 6)
# points(freqs.p[2,indexs$ix]*30,col= 7)
# points(freqs.p[3,indexs$ix]*30,col= 8)
# points(freqs.p[4,indexs$ix]*30,col= 2)
#

#View((proceeded$p.post-truth)[c(12,14,10,8,6,7,13,11,15,5,9,2,1,3,4)])
#View((proceeded$p.post^2+truth^2-2*proceeded$p.post*truth)[c(12,14,10,8,6,7,13,11,15,5,9,2,1,3,4)])
# cor(cols.b,cols.r)
# cor(cols.b,inits)
# cor(cols.r,inits)
cor(masses,freqs.p[1,1:243])
cor(masses,freqs.p[2,1:243])
cor(masses,freqs.p[3,1:243])
cor(masses,freqs.p[4,1:243])
cor(masses,freqs.p[5,1:243])

#View(statistics1[which(!is.na(statistics1[,1]))])

# build a package
#package.skeleton(name ="EMJMCMC", code_files =paste(workdir,"mode_jumping_package_class.r",sep = "",collapse = ""))
freqs.p[,,order.deviat$ix[1:Nlim]]
freqs[,order.deviat$ix[1:Nlim]]
iterats[1,order.deviat$ix[1:Nlim]]
iterats[2,order.deviat$ix[1:Nlim]]
inits[order.deviat$ix[1:Nlim]]

statistics1 <- big.matrix(nrow = 2 ^(length(fparam.example))+1, ncol = 15,init = NA, type = "double")
statistics <- describe(statistics1)
mySearch$g.results[4,1]<-0
mySearch$g.results[4,2]<-0
mySearch$p.add = array(data = 0.5,dim = 15)


mySearch$max.N.glob=as.integer(4)
mySearch$min.N.glob=as.integer(4)
mySearch$max.N=as.integer(2)
mySearch$min.N=as.integer(2)
mySearch$recalc.margin = as.integer(400)
distrib_of_neighbourhoods=freqs.p[,,order.deviat$ix[1]]
distrib_of_proposals = freqs[,order.deviat$ix[1]]
initsol = inits[order.deviat$ix[1]]
resm<-mySearch$modejumping_mcmc(list(varcur=mySearch$dectobit(inits[order.deviat$ix[1]]),statid=5, distrib_of_proposals = distrib_of_proposals,distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = 0.0001, trit = 3273, trest = 3272 , burnin = 3, max.time = 30, maxit = 100000, print.freq =500))
mySearch$g.results[4,2]
mySearch$g.results[4,1]
resm$bayes.results$s.mass/truth.prob

mean(iterats[2,order.deviat$ix[1:Nlim]])
