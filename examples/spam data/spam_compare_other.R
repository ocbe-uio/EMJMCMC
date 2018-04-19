library(RCurl)
library(glmnet)
library(xgboost)
library(h2o)
library(BAS)
library(caret)

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

h2o.removeAll()


for(ii in 1:100)
{
  print(paste("iteration ",ii))
  capture.output({withRestarts(tryCatch(capture.output({
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


  })), abort = function(){onerr<-TRUE;out<-NULL})})
  print(results[,ii,1])
}

ids<-NULL
for(i in 1:100)
{
  if(min(results[2,i,1])>0)
    ids<-c(ids,i)

}



ress<-results[,ids,]

summary.results<-array(data = NA,dim = c(15,15))
for(i in 1:11)
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
rownames(summary.results)[1:11]<-c("GMJMCMC(AIC)","MJMCMC(AIC)","lXGBOOST(logLik)","tXGBOOST(logLik)","LASSO","RIDGE","RFOREST","DEEPNETS","LR","NAIVEBAYESS","KMEANS")


write.csv(x = summary.results,file = "sumresothers.csv")
