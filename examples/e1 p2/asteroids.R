#compare on Asteroids
library(RCurl)

#define your working directory, where the data files are stored
workdir<-""

#prepare the test set data
simx <- read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/asteroid%20data/Recognize/NEAs.txt"),sep = ",",header = T,fill=TRUE)
simy <-  read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/asteroid%20data/Recognize/NotNeas8%2B.txt"),sep = ",",header = T,fill=TRUE)
simx$neo<-1
simy$neo<-0
data.example1 <- as.data.frame(t(cbind(t(simy),t(simx))),stringsAsFactors = T)
transform<-colnames(data.example1)[-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25)]
for(i in 1:length(transform))
{
  print(i)
  data.example1[[transform[i]]]<-as.numeric(as.character(data.example1[[transform[i]]]))
}


#prepare the training set data
simx <- read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/asteroid%20data/Teach/NeoPHA.txt"),sep = ",",header = T,fill=TRUE)
simy <-  read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/asteroid%20data/Teach/NotNeo-Type7.txt"),sep = ",",header = T,fill=TRUE)
simx$neo<-1
simy$neo<-0
data.example <- as.data.frame(t(cbind(t(simy),t(simx))),stringsAsFactors = T)
for(i in 1:length(transform))
{
  print(i)
  data.example[[transform[i]]]<-as.numeric(as.character(data.example[[transform[i]]]))
}


gc()

#GMJMCMC
system.time({
  
  formula1 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)],collapse = "+")))
  
  res = runemjmcmc(formula = formula1,data = data.example,estimator =estimate.bas.glm,estimator.args =  list(data = data.example,prior = aic.prior(),family = binomial(), logn = log(64)),recalc_margin = 50, save.beta = T,interact = T,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=2,mutation_rate = 100, max.tree.size = 200000, Nvars.max = 16,p.allow.replace=0.1,p.allow.tree=0.1,p.nor=0.3,p.and = 0.7),n.models = 5000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 100,advanced.param = list(
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

(1-sum(abs(res-data.example1$neo),na.rm = T)/20720)


#FNR
ps<-which(data.example1$neo==1)
sum(abs(res[ps]-data.example1$neo[ps]))/(sum(abs(res[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
sum(abs(res[ns]-data.example1$neo[ns]))/(sum(abs(res[ns]-data.example1$neo[ns]))+length(ns))

#MJMCMC
system.time({
  
  formula1 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)],collapse = "+")))
  
  res = runemjmcmc(formula = formula1,data = data.example,estimator =estimate.bas.glm,estimator.args =  list(data = data.example,prior = aic.prior(),family = binomial(), logn = log(64)),recalc_margin = 50, save.beta = T,interact = F,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=2,mutation_rate = 100, max.tree.size = 200000, Nvars.max = 16,p.allow.replace=0.1,p.allow.tree=0.1,p.nor=0.3,p.and = 0.7),n.models = 450,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 100,advanced.param = list(
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

(1-sum(abs(res-data.example1$neo),na.rm = T)/20720)


#FNR
ps<-which(data.example1$neo==1)
sum(abs(res[ps]-data.example1$neo[ps]))/(sum(abs(res[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
sum(abs(res[ns]-data.example1$neo[ns]))/(sum(abs(res[ns]-data.example1$neo[ns]))+length(ns))


gc()

#GMJMCMC elnet
system.time({
  
  formula1 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)],collapse = "+")))
  
  res = runemjmcmc(formula = formula1,data = data.example,estimator =estimate.elnet,estimator.args =  list(data = data.example,family = "binomial", response = "neo",alpha = 1),recalc_margin = 50, save.beta = T,interact = T,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=2,mutation_rate = 100, max.tree.size = 200000, Nvars.max = 16,p.allow.replace=0.1,p.allow.tree=0.4,p.nor=0.3,p.and = 0.7),n.models = 4000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 100,advanced.param = list(
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

(1-sum(abs(res-data.example1$neo),na.rm = T)/20720)


#FNR
ps<-which(data.example1$neo==1)
sum(abs(res[ps]-data.example1$neo[ps]))/(sum(abs(res[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
sum(abs(res[ns]-data.example1$neo[ns]))/(sum(abs(res[ns]-data.example1$neo[ns]))+length(ns))

#MJMCMC elnet
system.time({
  
  formula1 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)],collapse = "+")))
  
  res = runemjmcmc(formula = formula1,data = data.example,estimator =estimate.elnet,estimator.args =  list(data = data.example,family = "binomial", response = "neo",alpha = 1),recalc_margin = 50, save.beta = T,interact = F,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=2,mutation_rate = 100, max.tree.size = 200000, Nvars.max = 16,p.allow.replace=0.1,p.allow.tree=0.1,p.nor=0.3,p.and = 0.7),n.models = 450,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 100,advanced.param = list(
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

(1-sum(abs(res-data.example1$neo),na.rm = T)/20720)


#FNR
ps<-which(data.example1$neo==1)
sum(abs(res[ps]-data.example1$neo[ps]))/(sum(abs(res[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
sum(abs(res[ns]-data.example1$neo[ns]))/(sum(abs(res[ns]-data.example1$neo[ns]))+length(ns))


library(xgboost)
#xGboost AUC gbtree

param <- list(objective = "binary:logistic", 
              eval_metric = "auc",
              booster = "gbtree", 
              eta = 0.05,
              subsample = 0.86,
              colsample_bytree = 0.92,
              colsample_bylevel = 0.9,
              min_child_weight = 0,
              gamma = 0.005,
              max_depth = 15)

train<-as.data.frame(data.example[,-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)])
test<-as.data.frame(data.example1[,-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)])

dval<-xgb.DMatrix(data = data.matrix(train[,-1]), label = data.matrix(train[,1]),missing=NA)
watchlist<-list(dval=dval)

system.time({
m2 <- xgb.train(data = xgb.DMatrix(data = data.matrix(train[,-1]), label = data.matrix(train[,1]),missing=NA), 
                param, nrounds = 10000,
                watchlist = watchlist,
                print_every_n = 10)

})
# Predict

system.time({
dtest  <- xgb.DMatrix(data.matrix(test[,-1]),missing=NA)
})


system.time({
  out <- predict(m2, dtest)
})

out<-as.integer(out>=0.5)

print((1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))

# xgboost logLik gbtree

param <- list(objective = "binary:logistic", 
              eval_metric = "logloss",
              booster = "gbtree", 
              eta = 0.05,
              subsample = 0.86,
              colsample_bytree = 0.92,
              colsample_bylevel = 0.9,
              min_child_weight = 0,
              gamma = 0.005,
              max_depth = 15)

train<-as.data.frame(data.example[,-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)])
test<-as.data.frame(data.example1[,-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)])

dval<-xgb.DMatrix(data = data.matrix(train[,-1]), label = data.matrix(train[,1]),missing=NA)
watchlist<-list(dval=dval)

system.time({
  m2 <- xgb.train(data = xgb.DMatrix(data = data.matrix(train[,-1]), label = data.matrix(train[,1]),missing=NA), 
                  param, nrounds = 10000,
                  watchlist = watchlist,
                  print_every_n = 10)
  
})
# Predict
system.time({
  dtest  <- xgb.DMatrix(data.matrix(test[,-1]),missing=NA)
})

system.time({
out <- predict(m2, dtest)
})
out<-as.integer(out>=0.5)
print((1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))


#xGboost AUC gbtree

param <- list(objective = "binary:logistic", 
              eval_metric = "auc",
              booster = "gblinear", 
              eta = 0.05,
              subsample = 0.86,
              colsample_bytree = 0.92,
              colsample_bylevel = 0.9,
              min_child_weight = 0,
              gamma = 0.005,
              max_depth = 15)

as.h2o(test[,-1])
dval<-xgb.DMatrix(data = data.matrix(train[,-1]), label = data.matrix(train[,1]),missing=NA)
watchlist<-list(dval=dval)

system.time({
  m2 <- xgb.train(data = xgb.DMatrix(data = data.matrix(train[,-1]), label = data.matrix(train[,1]),missing=NA), 
                  param, nrounds = 10000,
                  watchlist = watchlist,
                  print_every_n = 10)
  
})
# Predict

system.time({
  dtest  <- xgb.DMatrix(data.matrix(test[,-1]),missing=NA)
})


system.time({
  out <- predict(m2, dtest)
})

out<-as.integer(out>=0.5)

print((1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))

# xgboost logLik gbtree

param <- list(objective = "binary:logistic", 
              eval_metric = "logloss",
              booster = "gblinear", 
              eta = 0.05,
              subsample = 0.86,
              colsample_bytree = 0.92,
              colsample_bylevel = 0.9,
              min_child_weight = 0,
              gamma = 0.005,
              max_depth = 15)

train<-as.data.frame(data.example[,-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)])
test<-as.data.frame(data.example1[,-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)])

dval<-xgb.DMatrix(data = data.matrix(train[,-1]), label = data.matrix(train[,1]),missing=NA)
watchlist<-list(dval=dval)

system.time({
  m2 <- xgb.train(data = xgb.DMatrix(data = data.matrix(train[,-1]), label = data.matrix(train[,1]),missing=NA), 
                  param, nrounds = 10000,
                  watchlist = watchlist,
                  print_every_n = 10)
  
})
# Predict
system.time({
  dtest  <- xgb.DMatrix(data.matrix(test[,-1]),missing=NA)
})

system.time({
  out <- predict(m2, dtest)
})
out<-as.integer(out>=0.5)
print((1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))


#GLMNET (elastic networks) # lasso a=1
library(glmnet)


system.time({
fit2 <- glmnet(as.matrix(train)[,-1], train$neo, family="binomial")
})


mmm<-as.matrix(test[,-1])
mmm[which(is.na(mmm))]<-0
system.time({
out <- predict(fit2,mmm , type = "response")[,fit2$dim[2]]
})


out<-as.integer(out>=0.5)

print((1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))

# ridge a=0

system.time({
  fit2 <- glmnet(as.matrix(train)[,-1], train$neo, family="binomial",alpha=0)
})


mmm<-as.matrix(test[,-1])
mmm[which(is.na(mmm))]<-0
system.time({
  out <- predict(fit2,mmm , type = "response")[,fit2$dim[2]]
})


out<-as.integer(out>=0.5)

print((1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))


# h2o.random forest
library(h2o)

# h2o initiate
h2o.init(nthreads=-1, max_mem_size = "6G")
h2o.removeAll()

df <- as.h2o(train)



train1 <- h2o.assign(df , "train1.hex")   
valid1 <- h2o.assign(df , "valid1.hex")
test1 <- h2o.assign(as.h2o(test[,-1]), "test1.hex")

train1[1:5,]

features = names(train1)[-1]

# in order to make the classification prediction
train1$neo <- as.factor(train1$neo)

system.time({
rf1 <- h2o.randomForest( stopping_metric = "AUC",   
                         training_frame = train1,       
                         validation_frame = valid1,     
                         x=features,                       
                         y="neo",                         
                         model_id = "rf1",   
                         ntrees = 10000,
                         stopping_rounds = 3,
                         score_each_iteration = T,
                         ignore_const_cols = T,
                         seed = 1234) 
})

system.time({
out<-h2o.predict(rf1,as.h2o(test1))[,1]
})

out<-as.data.frame(out)

out<-as.integer(as.numeric(as.character(out$predict)))


print((1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))

#h2o deeplearning

system.time({
neo.dl <- h2o.deeplearning(x = features, y = "neo",hidden=c(200,200,200,200,200,200),
                           distribution = "bernoulli", 
                           training_frame = train1,
                           validation_frame = valid1,
                           seed = 12341)
})
# now make a prediction


system.time({
  out<-h2o.predict(neo.dl,as.h2o(test1))[,1]
})

out<-as.data.frame(out)

out<-as.integer(as.numeric(as.character(out$predict)))


print((1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))


#h2o glm

system.time({
  neo.glm <- h2o.glm(x = features, y = "neo",
                             family = "binomial",
                             training_frame = train1,
                             validation_frame = valid1,
                             #lambda = 0,
                             #alpha = 0,
                             lambda_search = F,
                             seed = 12341)
})
# now make a prediction


system.time({
  out<-h2o.predict(neo.glm,as.h2o(test1))[,1]
})

out<-as.data.frame(out)

out<-as.integer(as.numeric(as.character(out$predict)))


print((1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))

#h2o naive bayes

system.time({
  neo.nb <- h2o.naiveBayes(x = features, y = "neo",
                           training_frame = train1,
                           validation_frame = valid1,
                           seed = 12341)
})
# now make a prediction


system.time({
  out<-h2o.predict(neo.nb,as.h2o(test1))[,1]
})

out<-as.data.frame(out)

out<-as.integer(as.numeric(as.character(out$predict)))


print((1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))

#h2o kmeans

system.time({
  neo.nb <- h2o.kmeans(x = c(features,"neo"),k=2,
                           training_frame = train1,
                           seed = 12341)
})


# now make a prediction


test2 <- h2o.assign(as.h2o(test), "test2.hex")

system.time({
  out<-h2o.predict(neo.dl,as.h2o(test2))[,1]
})

out<-as.data.frame(out)

out<-as.integer(as.numeric(as.character(out$predict)))


print((1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))



