library(RCurl)
library(caret)
#load library EMJMCMC implementing the MJMCMC and GMJMCMC algorithms for (D)GLMMs
#our article on MJMCMC is available at https://arxiv.org/abs/1604.06398, the one on GMJMCMC is coming soon
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package2.r")
#loading the data prepared
workdir<-""
pat.names<-read.table("submission_fin.csv",sep = ",",header = T,fill=TRUE)

patients <- read.table("data_head.csv",sep = ",",header = T,fill=TRUE)[,-1]
patients2 <- as.data.frame(prcomp(rbind(read.table("trainX.csv",sep = " ",header = F,fill=TRUE),read.table("testX.csv",sep = " ",header = F,fill=TRUE)),center = T)$x[,1:60])
patients3 <- read.table("data_face.csv",sep = ",",header = T,fill=TRUE)[,-c(1,12)]
patients4 <- read.table("data_profile.csv",sep = ",",header = T,fill=TRUE)[,-c(1,12)]

#test ids
submit<-which(patients[,106]==999)
#prepare test data from our best submission
test<-cbind(patients3[submit,-112],patients4[submit,-112],patients2[-(1:(dim(patients2)[1]-length(submit))),],patients[submit,-106])
test$label<-pat.names$cancer



#prepare train data
patients<-patients[-submit,]
patients3<-patients3[-submit,]
patients4<-patients4[-submit,]
patients2<-patients2[1:(dim(patients2)[1]-length(submit)),]
nas<-which(is.na(patients[,103]))
nas3<-which(is.na(patients3[,109]))
nas4<-which(is.na(patients4[,109]))
ty<-unique(union(union(nas,nas3),nas4))
patients<-cbind(patients3[-ty,-112],patients4[-ty,-112],patients2[-ty,],patients[-ty,])

#delete patients with missing data
nass<-(is.na(patients))
sum(nass)
sum(colSums(nass))
todel<-which(rowSums(nass)>0)
patients<-patients[-todel,]
nass<-(patients==Inf | patients == -Inf)
sum(nass)
sum(colSums(nass))
todel<-which(rowSums(nass)>0)
patients<-patients[-todel,]


nns<-rep("V",dim(patients)[2])
nns<-paste0(nns,1:dim(patients)[2])
nns[dim(patients)[2]]<-"label"
names(patients)<-nns
names(test)<-nns
test$label<-pat.names$cancer


#choose the best covariates
cor.t<-cor(test)
scor<-sort(abs(cor.t[-dim(patients)[2],dim(patients)[2]]),decreasing = T)
nvar<-names(scor[1:60])
# prepare MJMCMC
data.example <- as.data.frame(patients[,c(order(abs(cor.t[-dim(patients)[2],dim(patients)[2]]),decreasing = T)[1:60],dim(patients)[2])])
set.seed(1433)
#the penilized estimator function
estimate.bas.glm.cpen <- function(formula, data, family, prior, logn,r = 0.1,yid=1)
{
  
  #only poisson and binomial families are currently adopted
  X <- model.matrix(object = formula,data = data)
  out <- bayesglm.fit(x = X, y = data[,yid], family=family,coefprior=prior)
  fmla.proc<-as.character(formula)[2:3]
  fobserved <- fmla.proc[1]
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam <-stri_split_fixed(str = fmla.proc[2],pattern = "+I",omit_empty = F)[[1]]
  sj<-(stri_count_fixed(str = fparam, pattern = "("))
  mlik = (-(out$deviance -2*log(r)*sum(sj)))/2
  
  return(list(mlik = mlik,waic = -(out$deviance + 2*out$rank) , dic =  -(out$deviance + logn*out$rank),summary.fixed =list(mean = coefficients(out))))
  
}
#MJMCMC
t<-system.time({
  
  formula1 = as.formula("label ~ 1 + V239 + V276 + V323 + V231 + V308 + V324 + V325 + 
    V315 + V317 + V316 + V359 + V320 + V372 + V274 + V277 + V268 + 
    V312 + V225 + V373 + V240 + V360 + V321 + V387 + V370 + V243 + 
    V313 + V365 + V251 + V228 + V385 + V357 + V280 + V246 + V247 + 
    V376 + V368 + V261 + V378 + V238 + V332 + V329 + V355 + V363 + 
    V296 + V386 + V338 + V379 + V267 + V237 + V294 + V307 + V384 + 
    V257 + V278 + V358 + V302 + V260 + V366 + V345 + V342")
  
  res = runemjmcmc(formula = formula1,data = data.example,presearch=T, locstop =T,estimator =estimate.bas.glm.cpen,estimator.args =  list(data = data.example,prior = aic.prior(),family = binomial(), logn = dim(data.example)[1],r=exp(-1),yid=dim(data.example)[2]),recalc_margin = 50, save.beta = T,interact = F,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=2,last.mutation=1000,mutation_rate = 100, max.tree.size = 200000, Nvars.max = 21,p.allow.replace=0.1,p.allow.tree=0.1,p.nor=0.3,p.and = 0.7),n.models = 50000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
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


Nvars<-mySearch$Nvars
linx <-mySearch$Nvars+4
lHash<-length(hashStat)
mliks <- values(hashStat)[which((1:(lHash * linx)) %% linx == 1)]
betas <- values(hashStat)[which((1:(lHash * linx)) %% linx == 4)]
for(i in 1:(Nvars-1))
{
  betas<-cbind(betas,values(hashStat)[which((1:(lHash * linx)) %% linx == (4+i))])
}
betas<-cbind(betas,values(hashStat)[which((1:(lHash * linx)) %% linx == (0))])


t<-system.time({
  
  res<-mySearch$forecast.matrix.na(link.g = g,covariates = (test[,-dim(patients)[2]]),betas = betas,mliks.in = mliks)$forecast
  
})


#res[239]=mean(res,na.rm=T)
res[which(is.na(res))]<-0.265
summary(res)

length(res)

p.pred<-res
p.pred[which(p.pred==0)]=10^(-5)
p.pred[which(p.pred==1)]=1-10^(-5)
logloss = -sum(test$label*log(p.pred)+(1-test$label)*log(1-p.pred))/(length(res))


#save the results
pat.names$cancer<-p.pred
write.csv(x = pat.names,file = "stage1_sample_submission_22.csv",row.names = F)

