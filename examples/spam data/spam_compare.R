library(RCurl)
library(glmnet)
library(xgboost)
library(h2o)
library(BAS)
library(caret)
#define your working directory, where the data files are stored
workdir<-""
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")
# define the function estimating individual DBRM models
estimate.bas.glm.cpen <- function(formula, data, link, distribution, family, prior, logn,r = 0.1,yid=1,relat=c("cosi","sigmoid","tanh","atan","sini","troot"))
{
  capture.output({out <- glm(family = family,formula = formula,data = data)})
  fmla.proc<-as.character(formula)[2:3]
  fobserved <- fmla.proc[1]
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  sj<-2*(stri_count_fixed(str = fmla.proc[2], pattern = "*"))
  sj<-sj+1*(stri_count_fixed(str = fmla.proc[2], pattern = "+"))
  for(rel in relat)
    sj<-sj+2*(stri_count_fixed(str = fmla.proc[2], pattern = rel))
  mlik = ((-out$deviance +2*log(r)*sum(sj)))/2
  return(list(mlik = mlik,waic = -(out$deviance + 2*out$rank) , dic =  -(out$deviance + logn*out$rank),summary.fixed =list(mean = coefficients(out))))
}

setwd("/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/spam data/")

data = read.table("spam.data",col.names=c(paste("x",1:57,sep=""),"X"))
data[,1:57] = scale(data[,1:57])

set.seed(1)
spam.traintest =  read.table("spam.traintest")#rbinom(n = dim(data)[1],size = 1,prob = 0.01)
train = data[spam.traintest==1,]
test = data[spam.traintest==0,]

data.example = train





results<-array(0,dim = c(11,100,5))
#GMJMCMC




# h2o initiate
h2o.init(nthreads=-1, max_mem_size = "6G")

troot<-function(x)abs(x)^(1/3)
sini<-function(x)sin(x/180*pi)
logi<-function(x)log(abs(x+0.1))
gfquar<-function(x)as.integer(x<quantile(x,probs = 0.25))
glquar<-function(x)as.integer(x>quantile(x,probs = 0.75))
gmedi<-function(x)as.integer(x>median(x))
cosi<-function(x)cos(x/180*pi)
gmean<-function(x)as.integer(x>mean(x))
gone<-function(x)as.integer(x>0)
gthird<-function(x)(abs(x)^(1/3))
gfifth<-function(x)(abs(x)^(1/5))
grelu<-function(x)(x*(x>0))
contrelu<-function(x)log(1+exp(x))

h2o.removeAll()



featgmj = hash()
featrgmj = hash()

for(ii in 1:100)
{
  print(paste("iteration ",ii))
  capture.output({withRestarts(tryCatch(capture.output({

  set.seed(ii)

  #set.seed(runif(1,1,10000))
  t<-system.time({

    formula1 = as.formula(paste(colnames(data.example)[58],"~ 1 +",paste0(colnames(data.example)[-58],collapse = "+")))
    #gen.prob =c(1,1,1,1,1)
    res = runemjmcmc(formula = formula1,data = data.example,gen.prob = c(1,1,1,1,0),estimator =estimate.bas.glm.cpen,estimator.args =  list(data = data.example,prior = aic.prior(),link = "sigmoid", distribution = "binomial" ,family = binomial(),yid=58, logn = log(155),r=exp(-0.5)),recalc_margin = 95, save.beta = T, deep.method = 3,interact = T,relations = c("cosi","sigmoid","tanh","atan","sin","contrelu"),relations.prob =c(0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=3,mutation_rate = 100,last.mutation=1000, max.tree.size = 60, Nvars.max = 100,p.allow.replace=0.5,p.allow.tree=0.4,p.nor=0.3,p.and = 0.9),n.models = 20000,unique =F,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 100,advanced.param = list(
      max.N.glob=as.integer(20),
      min.N.glob=as.integer(5),
      max.N=as.integer(3),
      min.N=as.integer(1),
      printable = F))
  })

  results[1,ii,4]<-t[3]

  ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
  cbind(ppp$p.post, mySearch$fparam)

  mySearch$g.results[,]



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

    res<-mySearch$forecast.matrix.na(link.g = g,covariates = (test),betas = betas,mliks.in = mliks)$forecast

  })

  results[1,ii,5]<-t[3]

  summary(res)

  length(res)
  res<-as.integer(res>=0.5)
  length(which(res>=0.5))
  length(which(res<0.5))
  length(res)


  #(1-sum(abs(res-test$X),na.rm = T)/length(test$X))

  results[1,ii,1]<-(1-sum(abs(res-test$X),na.rm = T)/length(test$X))
  print(results[1,ii,1])
  gc()
  #FPR
  ns<-which(test$X==0)
  results[1,ii,3]<-sum(abs(res[ns]-test$X[ns]))/(sum(abs(res[ns]-test$X[ns]))+length(ns))

  #FNR
  ps<-which(test$X==1)
  results[1,ii,2]<-sum(abs(res[ps]-test$X[ps]))/(sum(abs(res[ps]-test$X[ps]))+length(ps))

  for(i in which(ppp$p.post>0.1))
  {  if(!has.key(hash = featgmj,key =  mySearch$fparam[i]))
    featgmj[[mySearch$fparam[i]]] = as.numeric(1) else{
      featgmj[[mySearch$fparam[i]]] =as.numeric(featgmj[[mySearch$fparam[i]]]) + 1
    }
  }
  gc()
  })), abort = function(){onerr<-TRUE;out<-NULL})})
  print(results[1,ii,1])
  }

  #MJMCMC
  t<-system.time({


    formula1 = as.formula(paste(colnames(data.example)[58],"~ 1 +",paste0(colnames(data.example)[-58],collapse = "+")))
    #gen.prob =c(1,1,1,1,1)
    res = runemjmcmc(formula = formula1,data = data.example,estimator =estimate.bas.glm.cpen,estimator.args =  list(data = data.example,prior = aic.prior(),family = binomial(),yid=58, logn = log(64),r=exp(-0.5)),recalc_margin = 50, save.beta = T,interact = F,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=2,last.mutation=1000,mutation_rate = 100, max.tree.size = 200000, Nvars.max = 16,p.allow.replace=0.1,p.allow.tree=0.1,p.nor=0.3,p.and = 0.7),n.models = 450,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
      max.N.glob=as.integer(10),
      min.N.glob=as.integer(5),
      max.N=as.integer(3),
      min.N=as.integer(1),
      printable = F))
  })

  results[2,ii,4]<-t[3]

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

    res<-mySearch$forecast.matrix.na(link.g = g,covariates = (test),betas = betas,mliks.in = mliks)$forecast

  })

  results[2,ii,5]<-t[3]

  summary(res)

  length(res)
  res<-as.integer(res>=0.5)
  length(which(res>=0.5))
  length(which(res<0.5))
  length(res)
  length(which(test$X==1))


  results[2,ii,1]<-(1-sum(abs(res-test$X),na.rm = T)/length(test$X))
  print(results[1,ii,1])
  gc()

  #FNR
  ps<-which(test$X==1)
  results[2,ii,2]<-sum(abs(res[ps]-test$X[ps]))/(sum(abs(res[ps]-test$X[ps]))+length(ps))

  #FPR
  ns<-which(test$X==0)
  results[2,ii,3]<-sum(abs(res[ns]-test$X[ns]))/(sum(abs(res[ns]-test$X[ns]))+length(ns))


  gc()




  #xGboost logloss gblinear
  t<-system.time({
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


  dval<-xgb.DMatrix(data = data.matrix(train[,-58]), label = data.matrix(train[,58]),missing=NA)
  watchlist<-list(dval=dval)


  m2 <- xgb.train(data = xgb.DMatrix(data = data.matrix(train[,-58]), label = data.matrix(train[,58]),missing=NA),
                  param, nrounds = 10000,
                  watchlist = watchlist,
                  print_every_n = 10)

  })
  # Predict
  results[3,ii,4]<-t[3]
  t<-system.time({
  dtest  <- xgb.DMatrix(data.matrix(test[,-58]),missing=NA)
  })


  t<-system.time({
    out <- predict(m2, dtest)
  })
  results[3,ii,5]<-t[3]
  out<-as.integer(out>=0.5)

  print( results[3,ii,1]<-(1-sum(abs(out-test$X[1:length(out)]))/length(out)))

  #FNR
  ps<-which(test$X==1)
  results[3,ii,2]<-sum(abs(out[ps]-test$X[ps]))/(sum(abs(out[ps]-test$X[ps]))+length(ps))

  #FPR
  ns<-which(test$X==0)
  results[3,ii,3]<-sum(abs(out[ns]-test$X[ns]))/(sum(abs(out[ns]-test$X[ns]))+length(ns))

  # xgboost logLik gbtree
  t<-system.time({
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

  dval<-xgb.DMatrix(data = data.matrix(train[,-58]), label = data.matrix(train[,58]),missing=NA)
  watchlist<-list(dval=dval)


    m2 <- xgb.train(data = xgb.DMatrix(data = data.matrix(train[,-58]), label = data.matrix(train[,58]),missing=NA),
                    param, nrounds = 10000,
                    watchlist = watchlist,
                    print_every_n = 10)

  })

  results[4,ii,4]<-t[3]
  # Predict
  system.time({
    dtest  <- xgb.DMatrix(data.matrix(test[,-58]),missing=NA)
  })

  t<-system.time({
  out <- predict(m2, dtest)
  })
  out<-as.integer(out>=0.5)

  print(results[4,ii,1]<-(1-sum(abs(out-test$X[1:length(out)]))/length(out)))

  #FNR
  ps<-which(test$X==1)
  results[4,ii,2]<-sum(abs(out[ps]-test$X[ps]))/(sum(abs(out[ps]-test$X[ps]))+length(ps))

  #FPR
  ns<-which(test$X==0)
  results[4,ii,3]<-sum(abs(out[ns]-test$X[ns]))/(sum(abs(out[ns]-test$X[ns]))+length(ns))


  #GLMNET (elastic networks) # lasso a=1



  t<-system.time({
  fit2 <- glmnet(as.matrix(train)[,-58], train$X, family="binomial")
  })
  results[5,ii,4]<-t[3]

  mmm<-as.matrix(test[,-58])
  mmm[which(is.na(mmm))]<-0
  t<-system.time({
  out <- predict(fit2,mmm , type = "response")[,fit2$dim[2]]
  })
  results[5,ii,5]<-t[3]

  out<-as.integer(out>=0.5)

  print(results[5,ii,1]<-(1-sum(abs(out-test$X[1:length(out)]))/length(out)))

  #FNR
  ps<-which(test$X==1)
  results[5,ii,2]<-sum(abs(out[ps]-test$X[ps]))/(sum(abs(out[ps]-test$X[ps]))+length(ps))

  #FPR
  ns<-which(test$X==0)
  results[5,ii,3]<-sum(abs(out[ns]-test$X[ns]))/(sum(abs(out[ns]-test$X[ns]))+length(ns))

  # ridge a=0

  t<-system.time({
    fit2 <- glmnet(as.matrix(train)[,-58], train$X, family="binomial",alpha=0)
  })
  results[6,ii,4]<-t[3]

  mmm<-as.matrix(test[,-58])
  mmm[which(is.na(mmm))]<-0
  t<-system.time({
    out <- predict(fit2,mmm , type = "response")[,fit2$dim[2]]
  })

  results[6,ii,5]<-t[3]

  out<-as.integer(out>=0.5)

  print(results[6,ii,1]<-(1-sum(abs(out-test$X[1:length(out)]))/length(out)))

  #FNR
  ps<-which(test$X==1)
  results[6,ii,2]<-sum(abs(out[ps]-test$X[ps]))/(sum(abs(out[ps]-test$X[ps]))+length(ps))

  #FPR
  ns<-which(test$X==0)
  results[6,ii,3]<-sum(abs(out[ns]-test$X[ns]))/(sum(abs(out[ns]-test$X[ns]))+length(ns))

  gc()

  # h2o.random forest



  df <- as.h2o(train)



  train1 <- h2o.assign(df , "train1.hex")
  valid1 <- h2o.assign(df , "valid1.hex")
  test1 <- h2o.assign(as.h2o(test[,-58]), "test1.hex")

  train1[1:5,]

  features = names(train1)[-58]

  # in order to make the classification prediction
  train1$X <- as.factor(train1$X)

  t<-system.time({
  rf1 <- h2o.randomForest( stopping_metric = "AUC",
                           training_frame = train1,
                           validation_frame = valid1,
                           x=features,
                           y="X",
                           model_id = "rf1",
                           ntrees = 10000,
                           stopping_rounds = 3,
                           score_each_iteration = T,
                           ignore_const_cols = T,
                           seed = ii)
  })
  results[7,ii,4]<-t[3]
  t<-system.time({
  out<-h2o.predict(rf1,as.h2o(test1))[,1]
  })
  results[7,ii,5]<-t[3]
  out<-as.data.frame(out)

  out<-as.integer(as.numeric(as.character(out$predict)))


  print(results[7,ii,1]<-(1-sum(abs(out-test$X[1:length(out)]))/length(out)))

  #FNR
  ps<-which(test$X==1)
  results[7,ii,2]<-sum(abs(out[ps]-test$X[ps]))/(sum(abs(out[ps]-test$X[ps]))+length(ps))

  #FPR
  ns<-which(test$X==0)
  results[7,ii,3]<-sum(abs(out[ns]-test$X[ns]))/(sum(abs(out[ns]-test$X[ns]))+length(ns))

  #h2o deeplearning

  t<-system.time({
  neo.dl <- h2o.deeplearning(x = features, y = "X",hidden=c(200,200,200,200,200,200),
                             distribution = "bernoulli",
                             training_frame = train1,
                             validation_frame = valid1,
                             seed = ii)
  })
  # now make a prediction

  results[8,ii,4]<-t[3]
  t<-system.time({
    out<-h2o.predict(neo.dl,as.h2o(test1))[,1]
  })
  results[8,ii,5]<-t[3]
  out<-as.data.frame(out)

  out<-as.integer(as.numeric(as.character(out$predict)))


  print(results[8,ii,1]<-(1-sum(abs(out-test$X[1:length(out)]))/length(out)))

  #FNR
  ps<-which(test$X==1)
  results[8,ii,2]<-sum(abs(out[ps]-test$X[ps]))/(sum(abs(out[ps]-test$X[ps]))+length(ps))

  #FPR
  ns<-which(test$X==0)
  results[8,ii,3]<-sum(abs(out[ns]-test$X[ns]))/(sum(abs(out[ns]-test$X[ns]))+length(ns))


  #h2o glm

  t<-system.time({
    neo.glm <- h2o.glm(x = features, y = "X",
                               family = "binomial",
                               training_frame = train1,
                               validation_frame = valid1,
                               #lambda = 0,
                               #alpha = 0,
                               lambda_search = F,
                               seed = ii)
  })
  # now make a prediction
  results[9,ii,4]<-t[3]

  t<-system.time({
    out<-h2o.predict(neo.glm,as.h2o(test1))[,1]
  })
  results[9,ii,5]<-t[3]
  out<-as.data.frame(out)

  out<-as.integer(as.numeric(as.character(out$predict)))


  print(results[9,ii,1]<-(1-sum(abs(out-test$X[1:length(out)]))/length(out)))

  #FNR
  ps<-which(test$X==1)
  results[9,ii,2]<-sum(abs(out[ps]-test$X[ps]))/(sum(abs(out[ps]-test$X[ps]))+length(ps))

  #FPR
  ns<-which(test$X==0)
  results[9,ii,3]<-sum(abs(out[ns]-test$X[ns]))/(sum(abs(out[ns]-test$X[ns]))+length(ns))

  #h2o naive bayes

  t<-system.time({
    neo.nb <- h2o.naiveBayes(x = features, y = "X",
                             training_frame = train1,
                             validation_frame = valid1,
                             seed = ii)
  })
  # now make a prediction

  results[10,ii,4]<-t[3]
  t<-system.time({
    out<-h2o.predict(neo.nb,as.h2o(test1))[,1]
  })
  results[10,ii,5]<-t[3]
  out<-as.data.frame(out)

  out<-as.integer(as.numeric(as.character(out$predict)))


  print(results[10,ii,1]<-(1-sum(abs(out-test$X[1:length(out)]))/length(out)))

  #FNR
  ps<-which(test$X==1)
  results[10,ii,2]<-sum(abs(out[ps]-test$X[ps]))/(sum(abs(out[ps]-test$X[ps]))+length(ps))

  #FPR
  ns<-which(test$X==0)
  results[10,ii,3]<-sum(abs(out[ns]-test$X[ns]))/(sum(abs(out[ns]-test$X[ns]))+length(ns))

  #h2o kmeans
  #set.seed(runif(1,1,10000))
  t<-system.time({

    formula1 = as.formula(paste(colnames(data.example)[58],"~ 1 +",paste0(colnames(data.example)[-58],collapse = "+")))
    #gen.prob =c(1,1,1,1,1)
    res = runemjmcmc(formula = formula1,data = data.example,gen.prob = c(1,1,1,1,0),estimator =estimate.bas.glm.cpen,estimator.args =  list(data = data.example,prior = aic.prior(),family = binomial(),yid=58, logn = log(155),r=exp(-0.5)),recalc_margin = 95, save.beta = T,interact = T,relations = c("cosi","sigmoid","tanh","atan","sini","troot"),relations.prob =c(0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=4,mutation_rate = 100,last.mutation=1000, max.tree.size = 6, Nvars.max = 100,p.allow.replace=0.5,p.allow.tree=0.4,p.nor=0.3,p.and = 0.9),n.models = 20000,unique =F,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 100,advanced.param = list(
      max.N.glob=as.integer(20),
      min.N.glob=as.integer(5),
      max.N=as.integer(3),
      min.N=as.integer(1),
      printable = F))
  })

  results[11,ii,4]<-t[3]

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

    res<-mySearch$forecast.matrix.na(link.g = g,covariates = (test),betas = betas,mliks.in = mliks)$forecast

  })

  results[11,ii,5]<-t[3]

  summary(res)

  length(res)
  res<-as.integer(res>=0.5)
  length(which(res>=0.5))
  length(which(res<0.5))
  length(res)


  #(1-sum(abs(res-test$X),na.rm = T)/length(test$X))

  for(i in which(ppp$p.post>0.1))
  {  if(!has.key(hash = featrgmj,key =  mySearch$fparam[i]))
    featrgmj[[mySearch$fparam[i]]] = as.numeric(1) else{
      featrgmj[[mySearch$fparam[i]]] =as.numeric(featrgmj[[mySearch$fparam[i]]]) + 1
    }
  }

  results[11,ii,1]<-(1-sum(abs(res-test$X),na.rm = T)/length(test$X))
  print(results[1,ii,1])
  gc()
  #FPR
  ns<-which(test$X==0)
  results[11,ii,3]<-sum(abs(res[ns]-test$X[ns]))/(sum(abs(res[ns]-test$X[ns]))+length(ns))

  #FNR
  ps<-which(test$X==1)
  results[11,ii,2]<-sum(abs(res[ps]-test$X[ps]))/(sum(abs(res[ps]-test$X[ps]))+length(ps))
  gc()

  })), abort = function(){onerr<-TRUE;out<-NULL})})
  print(results[,ii,1])
}

ids<-NULL
for(i in 1:100)
{
  if(min(results[1,i,1])>0)
    ids<-c(ids,i)

}



ress<-results[,ids,]

summary.results<-array(data = NA,dim = c(15,15))
for(i in 1:1)
{
  for(j in 1:5)
  {
    summary.results[i,(j-1)*3+1]<-min(ress[i,,j])
    summary.results[i,(j-1)*3+2]<-median(ress[i,,j])
    summary.results[i,(j-1)*3+3]<-max(ress[i,,j])
  }
}
summary.results<-as.data.frame(summary.results)
names(summary.results)<-c("min(prec)","median(prec)","max(prec)","min(fnr)","median(fnr)","max(fnr)","min(fpr)","median(fpr)","max(fpr)","min(ltime)","median(ltime)","max(ltime)","min(ptime)","median(ptime)","max(ptime)")
rownames(summary.results)[1:11]<-c("GMJMCMC(AIC)","MJMCMC(AIC)","lXGBOOST(logLik)","tXGBOOST(logLik)","LASSO","RIDGE","RFOREST","DEEPNETS","NAIVEBAYESS","LR","KMEANS")



for(i in 1:15)
{
  plot(density(ress[i,,1],bw = "SJ"), main="Compare Kernel Density of precisions")
  polygon(density(ress[i,,1],bw = "SJ"), col="red", border="blue")

}

write.csv(x = cbind(keys(featgmj),values(featgmj)),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/neo classification/spamfeatgmj.csv")
write.csv(x = cbind(keys(featrgmj),values(featrgmj)),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/neo classification/spamfeatrgmj.csv")


