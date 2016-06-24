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
#bulid the sourse file https://github.com/aliaksah/EMJMCMC2016/blob/master/R/the_mode_jumping_package2.r instead, since the latest published build EMJMCMC_1.2.tar.gz doen not have prediction functionality (coming in the package soon with other cookies ;) ) available
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

#define the working directory

workdir<-""

# get the data
M<-5
size<-1

data.example2 <- read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/Epigenetic%20Data/epigen.txt"),sep = ",",header = T)[,2:30]
data.example2$strfact<-0
data.example2$strfact[which(data.example2$strand=="+")]<-1
data.example2$maexpr<-data.example2$Mα*data.example2$express
data.example2$mgexpr<-data.example2$Mγ*data.example2$express
data.example2$mdexpr<-data.example2$Mδ*data.example2$express


workdir = ""
data.example1<-data.example2[which(data.example2$total_bases<=2),]


data.example<-data.example2[which(data.example2$total_bases>2),]

plot(data.example$total_bases, col = 4, pch = 16,ylim = c(0,30),main = "Epigenetic data sample",xlab = "nucleobase id",ylab = "reads")
points(data.example$methylated_bases,col=2, pch = 17)

#specify variables
length(which(data.example2$Mα==1))
length(which(data.example2$Mβ==1))
length(which(data.example2$Mγ==1))
length(which(data.example2$Mδ==1))
length(which(data.example2$MIKC==1))
length(which(data.example2$none==1)) # I have found now an issue with none, which does not seem to be present!?

View(cor(data.example2[c(5,8:10,12:17,21,23,24,26,29,30)])) #and also Md is highly correlated with the expression levels, so the betas
# can be explained with overparametrization

fparam.example <- colnames(data.example2)[c(8:10,12:17,21,23,26,30,31,32,33)]
fobserved.example <- colnames(data.example2)[5]


# now apply logistic regression
#create MySearch object with default parameters. N/B default estimator is INLA!
mySearch = EMJMCMC2016()
#mySearch$parallelize = lapply



# specify some INLA realted parameters
mySearch$estimator = estimate.inla.ar1
args<-list(args = list(family = "binomial",data = data.example, Ntrials = data.example$total_bases))
args$args$control.compute = list(dic = TRUE, waic = TRUE, mlik = TRUE)
mySearch$latent.formula  = "+f(data.example$pos,model=\"ar1\")"
mySearch$estimator.args = args
mySearch$printable.opt = F
mySearch$save.beta=T

estimate.inla.ar1(formula = formula2,args)


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
mySearch$recalc.margin = as.integer(2^10)
distrib_of_proposals = c(76.91870,71.25264,87.68184,60.55921,15812.39852)


# draw some paprameters of the search
pie(c(distrib_of_proposals[5],sum(distrib_of_proposals[1:4])), labels = c("NO","YES"), main="Local optimizers addressed")
pie(distrib_of_proposals[1:4], labels = c("SA1", "LMCMC", "SA2", "GREEDY"), main="Local optimizers frequencies")


#distrib_of_proposals = с(0,0,0,0,10)
distrib_of_neighbourhoods=t(array(data = c(7.6651604,16.773326,14.541629,12.839445,2.964227,13.048343,7.165434,
                                           0.9936905,15.942490,11.040131,3.200394,15.349051,5.466632,14.676458,
                                           1.5184551,9.285762,6.125034,3.627547,13.343413,2.923767,15.318774,
                                           14.5295380,1.521960,11.804457,5.070282,6.934380,10.578945,12.455602,
                                           6.0826035,2.453729,14.340435,14.863495,1.028312,12.685017,13.806295),dim = c(7,5)))

# draw some paprameters of the search
pie(distrib_of_neighbourhoods[1,], labels = c("TYPE 1", "TYPE 2", "TYPE 3", "TYPE 4", "TYPE 5", "TYPE 6", "TYPE 7"), main="Proposal frequencies in SA1")
pie(distrib_of_neighbourhoods[2,], labels = c("TYPE 1", "TYPE 2", "TYPE 3", "TYPE 4", "TYPE 5", "TYPE 6", "TYPE 7"), main="Proposal frequencies in LMCMC")
pie(distrib_of_neighbourhoods[3,], labels = c("TYPE 1", "TYPE 2", "TYPE 3", "TYPE 4", "TYPE 5", "TYPE 6", "TYPE 7"), main="Proposal frequencies in SA2")
pie(distrib_of_neighbourhoods[4,], labels = c("TYPE 1", "TYPE 2", "TYPE 3", "TYPE 4", "TYPE 5", "TYPE 6", "TYPE 7"), main="Proposal frequencies in GREEDY")
pie(distrib_of_neighbourhoods[5,], labels = c("TYPE 1", "TYPE 2", "TYPE 3", "TYPE 4", "TYPE 5", "TYPE 6", "TYPE 7"), main="Proposal frequencies in NO")


#carry the search (training out)
statistics1 <- big.matrix(nrow = 2 ^(length(fparam.example))+1, ncol =  16+length(fparam.example),init = NA, type = "double")
statistics <- describe(statistics1)
mySearch$g.results[4,1]<-0
mySearch$g.results[4,2]<-0
mySearch$p.add = array(data = 0.5,dim = length(fparam.example))
initsol=rbinom(n = length(fparam.example),size = 1,prob = 0.5)
resm<-mySearch$modejumping_mcmc(list(varcur=initsol,statid=5, distrib_of_proposals =distrib_of_proposals,distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = 0.000000001, trit = 999000, trest = 1000, burnin = 3, max.time = 30, maxit = 100000, print.freq =50))


#analyze the results from the search
ppp<-mySearch$post_proceed_results(statistics1 = statistics1)
truth = ppp$p.post # make sure it is equal to Truth column from the article
truth.m = ppp$m.post
truth.prob = ppp$s.mass
ordering = sort(ppp$p.post,index.return=T)
print("pi truth")
sprintf("%.10f",truth[ordering$ix])
sprintf(fparam.example[ordering$ix])

par(mar = c(10,4,4,2) + 4.1)

template = "INLA"

jpeg(file=paste(workdir,template,"barrm.jpg",sep = ""))
barplot(resm$bayes.results$p.post,density = 46,border="black",main = "Marginal Inclusion (RM)",ylab="Probability",names.arg =fparam.example,las=2)
dev.off()
jpeg(file=paste(workdir,template,"barmc.jpg",sep = ""))
barplot(resm$p.post,density = 46,border="black",main = "Marginal Inclusion (MC)",ylab="Probability",names.arg =fparam.example,las=2)
dev.off()

# n/b visualization iisues on Windows! To be fixed!
#mySearch$visualize_results(statistics1, template, 200, crit=list(mlik = T, waic = T, dic = T),draw_dist =T)
View(statistics1[,])

best<-fparam.example[which(!is.na(statistics1[which(statistics1[,15]==max(statistics1[,15],na.rm = T))[1],17:(16+length(fparam.example))]))]
#best<-fparam.example[-c(1,3,5,6,7,8,9,15)]

# condiotion on posterior mode in  the model space for now
data.example3<-data.example2
idtest<-which(data.example2$total_bases<=2)
idtrain<-which(data.example2$total_bases>2)
data.example3$methylated_bases[idtest]<-NA
#example of the underlying model within INLA
formula2 <- as.formula(paste0("methylated_bases ~  1 +",paste(best,collapse = "+"), " + f(data.example3$pos,model=\"ar1\")"))
args<-list(family = "binomial",data = data.example3, Ntrials = data.example3$total_bases)
fm4<-do.call(inla, c(args,formula = formula2))
summary(fm4)

fm4$summary.hyperpar$mean[1]
coef<-fm4$summary.fixed$mode
coef[1]<-coef[1]+fm4$summary.hyperpar$mean[1]

# carry out glm BIC prior selection
mySearch$parallelize = mclapply
# specify some INLA realted parameters
mySearch$estimator = estimate.glm
mySearch$estimator.args = list(data = data.example,prior = 2,family = binomial,observ="cbind(methylated_bases,unmethylated_bases)")
#carry the search (training out)
statistics1 <- big.matrix(nrow = 2 ^(length(fparam.example))+1, ncol =  16+length(fparam.example),init = NA, type = "double")
statistics <- describe(statistics1)
mySearch$g.results[4,1]<-0
mySearch$g.results[4,2]<-0
mySearch$p.add = array(data = 0.5,dim = length(fparam.example))
initsol=rbinom(n = length(fparam.example),size = 1,prob = 0.5)
resm<-mySearch$modejumping_mcmc(list(varcur=initsol,statid=5, distrib_of_proposals =distrib_of_proposals,distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = 0.000000001, trit = 999000, trest = 10000, burnin = 3, max.time = 30, maxit = 100000, print.freq =50))

#analyze the results from the search
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
#mySearch$visualize_results(statistics1, "test", 200, crit=list(mlik = T, waic = T, dic = T),draw_dist =T)
best<-fparam.example[which(!is.na(statistics1[which(statistics1[,15]==max(statistics1[,15],na.rm = T))[1],17:(16+length(fparam.example))]))]


#now check the glm results
formula3 <- as.formula(paste0("cbind(methylated_bases,unmethylated_bases)~ 1 +",paste(best,collapse = "+")))
fit1 <- glm( formula = formula3,family=binomial,data=data.example)
summary(fit1)


g<-function(x)
{
  return((x = 1/(1+exp(-x))))
}

#classify the unknown

res1<-mySearch$foreast.matrix(nvars = 16,ncases = 552,link.g = g,covariates = data.example1[1:552,c(8:10,12:17,21,23,26,30,31,32,33)])$forecast

res<-as.integer(res1>=0.5)
length(which(res>=0.5))
length(which(res<0.5))
length(res)

sum(data.example1$methylation_call)
sum(abs(res-data.example1$methylation_call))
length(which(res-data.example1$methylation_call<0))
length(which(res-data.example1$methylation_call>0))
length(which(res-data.example1$methylation_call==0))


jpeg(file=paste(workdir,template,"predict.jpg",sep = ""))
plot(data.example1$total_bases, col = 4, pch = 16,ylim = c(-11,3),main = "Epigenetic data sample",xlab = "nucleobase id",ylab = "reads")
points(data.example1$methylated_bases,col=2, pch = 17)
lines(10*data.example1$methylation_call-11,col = 5)
lines(10*res1[1,] - 11,col=8)
lines(10*g(fm4$summary.random[[1]]$mode[idtest])-11,col =1,lwd=4)
lines(10*g(fm4$summary.random[[1]]$mean[idtest])-11,col =7,lwd=2)
dev.off()
#classify the known
res1<-mySearch$foreast.matrix(nvars = 16,ncases = 950,link.g = g,covariates = data.example[1:950,c(8:10,12:17,21,23,26,30,31,32,33)])$forecast
res<-as.integer(res1>=0.5)
length(which(res>=0.5))
length(which(res<0.5))
length(res)

sum(data.example$methylation_call)
sum(abs(res-data.example$methylation_call))

length(which(res-data.example$methylation_call<0))
length(which(res-data.example$methylation_call>0))
length(which(res-data.example$methylation_call==0))

jpeg(file=paste(workdir,template,"classify.jpg",sep = ""))
plot(data.example$total_bases, col = 4, pch = 16,ylim = c(-11,30),main = "Epigenetic data sample",xlab = "nucleobase id",ylab = "reads")
points(data.example$methylated_bases,col=2, pch = 17)
lines(10*data.example$methylation_call-11,col = 5)
lines(10*res1[1,] - 11,col=8)
lines(10*g(fm4$summary.random[[1]]$mode[idtrain])-11,col =1,lwd=4)
lines(10*g(fm4$summary.random[[1]]$mean[idtrain])-11,col =7,lwd=2)
dev.off()

# scale up interesting regions:

plot(data.example$total_bases[1:50], col = 4, pch = 16,ylim = c(-10,30),main = "Epigenetic data sample",xlab = "nucleobase id",ylab = "reads")
points(data.example$methylated_bases[1:50],col=2, pch = 17)
lines(10*data.example$methylation_call[1:50]-11,col = 5)
lines(10*res1[1,][1:50] - 11,col=8)
lines(10*g(fm4$summary.random[[1]]$mode[idtrain])[1:50]-11,col =1,lwd=4)
lines(10*g(fm4$summary.random[[1]]$mean[idtrain])[1:50]-11,col =7,lwd=2)

plot(data.example$total_bases[200:300], col = 4, pch = 16,ylim = c(-10,30),main = "Epigenetic data sample",xlab = "nucleobase id",ylab = "reads")
points(data.example$methylated_bases[200:300],col=2, pch = 17)
lines(10*data.example$methylation_call[200:300]-11,col = 5)
lines(10*res1[1,][200:300] - 11,col=8)
lines(10*g(fm4$summary.random[[1]]$mode[idtrain])[200:300]-11,col =1,lwd=4)
lines(10*g(fm4$summary.random[[1]]$mean[idtrain])[200:300]-11,col =7,lwd=2)



plot(data.example$total_bases[450:500], col = 4, pch = 16,ylim = c(-10,30),main = "Epigenetic data sample",xlab = "nucleobase id",ylab = "reads")
points(data.example$methylated_bases[450:500],col=2, pch = 17)
lines(10*data.example$methylation_call[450:500]-11,col = 5)
lines(10*res1[1,][450:500] - 11,col=8)
lines(10*g(fm4$summary.random[[1]]$mode[idtrain])[450:500]-11,col =1,lwd=4)
lines(10*g(fm4$summary.random[[1]]$mean[idtrain])[450:500]-11,col =7,lwd=2)


plot(data.example$total_bases[650:750], col = 4, pch = 16,ylim = c(-10,30),main = "Epigenetic data sample",xlab = "nucleobase id",ylab = "reads")
points(data.example$methylated_bases[650:750],col=2, pch = 17)
lines(10*data.example$methylation_call[650:750]-11,col = 5)
lines(10*res1[1,][650:750] - 11,col=8)
lines(10*g(fm4$summary.random[[1]]$mode[idtrain])[650:750]-11,col =1,lwd=4)
lines(10*g(fm4$summary.random[[1]]$mean[idtrain])[650:750]-11,col =7,lwd=2)



