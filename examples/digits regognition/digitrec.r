#compare on Asteroids
library(RCurl)
library(caret)
library("EBImage")

digits <- read.table(file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/digits regognition/train.csv",sep = ",",header = T,fill=TRUE)

digits$Y1<-as.integer(digits$label==1)
digits$Y2<-as.integer(digits$label==2)
digits$Y3<-as.integer(digits$label==3)
digits$Y4<-as.integer(digits$label==4)
digits$Y5<-as.integer(digits$label==5)
digits$Y6<-as.integer(digits$label==6)
digits$Y7<-as.integer(digits$label==7)
digits$Y8<-as.integer(digits$label==8)
digits$Y9<-as.integer(digits$label==9)
digits$Y10<-as.integer(digits$label==10)


data<-round(runif(10,1,42000))

df<-NULL
for (k in 1:42)
{


  matrix2 <-matrix(data = digits[k,2:785],nrow = 28,ncol = 28,byrow = F)

  matrix2<-normalize(matrix2)
  #image(matrix2, axes=FALSE, col=topo.colors(12))
  matrix2<-(matrix2>otsu(matrix2, range = c(0, 1), levels = 256))
  matrix2<-resize(matrix2, w = 12, h = 12)
  matrix2<-normalize(matrix2)
  matrix2<-(matrix2>otsu(matrix2, range = c(0, 1), levels = 256))
  df<-cbind(df,as.array(matrix2))
}




index <- createDataPartition(digits$label, p = 0.1, list = FALSE)

test <- digits[-index, ]
train <- digits[index, ]

data.example <- as.data.frame(train,stringsAsFactors = T)


gc()


gfquar<-function(x)as.integer(x<quantile(x,probs = 0.25))
glquar<-function(x)as.integer(x>quantile(x,probs = 0.75))
gmedi<-function(x)as.integer(x>median(x))
cosi<-function(x)cos(x/180*pi)
gmean<-function(x)as.integer(x>mean(x))
dig<-1

# h2o initiate
h2o.init(nthreads=-1, max_mem_size = "6G")

h2o.removeAll()

for(ii in 1:100)
{
  print(paste("iteration ",ii))
  capture.output({withRestarts(tryCatch(capture.output({

    set.seed(ii)

    idss<-which(abs(cor(x = train[1:785],y=train[785+dig]))>0.2)

    formula1 = as.formula(paste(colnames(train)[785+dig],"~ 1 +",paste0(colnames(train)[idss],collapse = "+")))
    #set.seed(runif(1,1,10000))
    #set.seed(runif(1,1,10000))
    t<-system.time({

      res = runemjmcmc(formula = formula1,data = data.example,estimator =estimate.bas.glm.cpen,estimator.args =  list(data = data.example,prior = aic.prior(),family = binomial(), logn = log(64),r=exp(-1),yid=785+dig),recalc_margin = 95,locstop=T,presearch=T,save.beta = T,interact = T,relations = c("","cosi","sigmoid","tanh","atan","erf","gmean","gmedi","gfquar","glquar"),relations.prob =c(0.8,0.1,0.1,0.1,0.1,0.1,0.1,0.5,0.1,0.1),interact.param=list(allow_offsprings=2,mutation_rate = 100,last.mutation=1000, max.tree.size = 4, Nvars.max =50,p.allow.replace=0.1,p.allow.tree=0.5,p.nor=0.3,p.and = 0.7),n.models = 15000,unique =F,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
        max.N.glob=as.integer(10),
        min.N.glob=as.integer(5),
        max.N=as.integer(3),
        min.N=as.integer(1),
        printable = F))
    })


    results[1,ii,4]<-t[3]

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

      res<-mySearch$forecast.matrix.na(link.g = g,covariates = (test[,idss]),betas = betas,mliks.in = mliks)$forecast

    })
    data<-runif(1,1,42000)

    for (k in data)
    {
      row <- NULL
      for (n in 1:784)
        row[n] <- digits[k,n]

      matrix1 <- matrix(row,28,28,byrow=FALSE)
      matrix2 <- matrix(rep(0,784),28,28)

      for (i in 1:28)
        for (j in 1:28)
          matrix2[i,28-j+1] <- matrix1[i,j]

      image(matrix2, axes=FALSE, col=topo.colors(12))
      matrix2[which(matrix2>10)]<-255
      matrix2[which(matrix2<10)]<-0
      image(matrix2, axes=FALSE, col=topo.colors(100))
    }

    index <- createDataPartition(digits$label, p = 0.05, list = FALSE)

    test <- digits[-index, ]
    train <- digits[index, ]

    results[1,ii,5]<-t[3]

    summary(res)

    length(res)
    res<-as.integer(res>=0.5)
    length(which(res>=0.5))
    length(which(res<0.5))
    length(res)
    length(which(test$label==1))

    #(1-sum(abs(res-data.example1$neo),na.rm = T)/20702)

    results[1,ii,1]<-(1-sum(abs(res-test$Y1),na.rm = T)/length(res))
    print(results[1,ii,1])
    gc()

    #FNR
    ps<-which(test$label==1)
    results[1,ii,2]<-sum(abs(res[ps]-test$label[ps]))/(sum(abs(res[ps]-test$label[ps]))+length(ps))

    #FPR
    ns<-which(test$label==0)
    results[1,ii,3]<-sum(abs(res[ns]-test$label[ns]))/(sum(abs(res[ns]-test$label[ns]))+length(ns))

    gc()
  })), abort = function(){onerr<-TRUE;out<-NULL})})
  print(results[1,ii,1])

}

#MJMCMC
t<-system.time({

  formula1 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)],collapse = "+")))

  res = runemjmcmc(formula = formula1,data = data.example,estimator =estimate.bas.glm,estimator.args =  list(data = data.example,prior = aic.prior(),family = binomial(), logn = log(64)),recalc_margin = 50, save.beta = T,interact = F,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=2,last.mutation=1000,mutation_rate = 100, max.tree.size = 200000, Nvars.max = 16,p.allow.replace=0.1,p.allow.tree=0.1,p.nor=0.3,p.and = 0.7),n.models = 450,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
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

  res<-mySearch$forecast.matrix.na(link.g = g,covariates = (data.example1[1:20702,-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)]),betas = betas,mliks.in = mliks)$forecast

})

results[2,ii,5]<-t[3]

summary(res)

length(res)
res<-as.integer(res>=0.5)
length(which(res>=0.5))
length(which(res<0.5))
length(res)
length(which(data.example1$neo==1))

results[2,ii,1]<-(1-sum(abs(res-data.example1$neo),na.rm = T)/20702)


#FNR
ps<-which(data.example1$neo==1)
results[2,ii,2]<-sum(abs(res[ps]-data.example1$neo[ps]))/(sum(abs(res[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
results[2,ii,3]<-sum(abs(res[ns]-data.example1$neo[ns]))/(sum(abs(res[ns]-data.example1$neo[ns]))+length(ns))


gc()

#GMJMCMC elnet
t<-system.time({

  formula1 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)],collapse = "+")))

  res = runemjmcmc(formula = formula1,data = data.example,estimator =estimate.elnet,estimator.args =  list(data = data.example,family = "binomial", response = "neo",alpha = 1),recalc_margin = 99, save.beta = T,interact = T,relations = c("","cosi","sigmoid","tanh","atan","erf"),relations.prob =c(0.8,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=2,mutation_rate = 100,last.mutation=500, max.tree.size = 4, Nvars.max = 15,p.allow.replace=0.4,p.allow.tree=0.4,p.nor=0.3,p.and = 0.7),n.models = 7000,unique = F,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
    max.N.glob=as.integer(10),
    min.N.glob=as.integer(5),
    max.N=as.integer(3),
    min.N=as.integer(1),
    printable = F))
})

results[3,ii,4]<-t[3]

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

  res<-mySearch$forecast.matrix.na(link.g = g,covariates = (data.example1[1:20702,-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)]),betas = betas,mliks.in = mliks)$forecast

})

results[3,ii,5]<-t[3]

summary(res)

length(res)
res<-as.integer(res>=0.5)
length(which(res>=0.5))
length(which(res<0.5))
length(res)
length(which(data.example1$neo==1))

results[3,ii,1]<-(1-sum(abs(res-data.example1$neo),na.rm = T)/20702)


#FNR
ps<-which(data.example1$neo==1)
results[3,ii,2]<-sum(abs(res[ps]-data.example1$neo[ps]))/(sum(abs(res[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
results[3,ii,3]<-sum(abs(res[ns]-data.example1$neo[ns]))/(sum(abs(res[ns]-data.example1$neo[ns]))+length(ns))


#MJMCMC elnet
t<-system.time({

  formula1 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)],collapse = "+")))

  res = runemjmcmc(formula = formula1,data = data.example,estimator =estimate.elnet,estimator.args =  list(data = data.example,family = "binomial", response = "neo",alpha = 1),recalc_margin = 50, save.beta = T,interact = F,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=2,mutation_rate = 100, max.tree.size = 200000,last.mutation=500, Nvars.max = 16,p.allow.replace=0.1,p.allow.tree=0.1,p.nor=0.3,p.and = 0.7),n.models = 450,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
    max.N.glob=as.integer(10),
    min.N.glob=as.integer(5),
    max.N=as.integer(3),
    min.N=as.integer(1),
    printable = F))
})

results[4,ii,4]<-t[3]
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

  res<-mySearch$forecast.matrix.na(link.g = g,covariates = (data.example1[1:20702,-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)]),betas = betas,mliks.in = mliks)$forecast

})
results[4,ii,5]<-t[3]
summary(res)

length(res)
res<-as.integer(res>=0.5)
length(which(res>=0.5))
length(which(res<0.5))
length(res)
length(which(data.example1$neo==1))

results[4,ii,1]<-(1-sum(abs(res-data.example1$neo),na.rm = T)/20702)


#FNR
ps<-which(data.example1$neo==1)
results[4,ii,2]<-sum(abs(res[ps]-data.example1$neo[ps]))/(sum(abs(res[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
results[4,ii,3]<-sum(abs(res[ns]-data.example1$neo[ns]))/(sum(abs(res[ns]-data.example1$neo[ns]))+length(ns))



#xGboost AUC gbtree
t<-system.time({
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


  m2 <- xgb.train(data = xgb.DMatrix(data = data.matrix(train[,-1]), label = data.matrix(train[,1]),missing=NA),
                  param, nrounds = 10000,
                  watchlist = watchlist,
                  print_every_n = 10)

})
# Predict
results[5,ii,4]<-t[3]
t<-system.time({
  dtest  <- xgb.DMatrix(data.matrix(test[,-1]),missing=NA)
})


t<-system.time({
  out <- predict(m2, dtest)
})
results[5,ii,5]<-t[3]
out<-as.integer(out>=0.5)

print( results[5,ii,1]<-(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
results[5,ii,2]<-sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
results[5,ii,3]<-sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))

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

  train<-as.data.frame(data.example[,-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)])
  test<-as.data.frame(data.example1[,-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)])

  dval<-xgb.DMatrix(data = data.matrix(train[,-1]), label = data.matrix(train[,1]),missing=NA)
  watchlist<-list(dval=dval)


  m2 <- xgb.train(data = xgb.DMatrix(data = data.matrix(train[,-1]), label = data.matrix(train[,1]),missing=NA),
                  param, nrounds = 10000,
                  watchlist = watchlist,
                  print_every_n = 10)

})

results[6,ii,4]<-t[3]
# Predict
system.time({
  dtest  <- xgb.DMatrix(data.matrix(test[,-1]),missing=NA)
})

t<-system.time({
  out <- predict(m2, dtest)
})
out<-as.integer(out>=0.5)
print(results[6,ii,1]<-(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))
results[6,ii,5]<-t[3]
#FNR
ps<-which(data.example1$neo==1)
results[6,ii,2]<-sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
results[6,ii,3]<-sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))


#xGboost AUC gbtree
t<- system.time({
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


  m2 <- xgb.train(data = xgb.DMatrix(data = data.matrix(train[,-1]), label = data.matrix(train[,1]),missing=NA),
                  param, nrounds = 10000,
                  watchlist = watchlist,
                  print_every_n = 10)

})
# Predict
results[7,ii,4]<-t[3]
t<-system.time({
  dtest  <- xgb.DMatrix(data.matrix(test[,-1]),missing=NA)
})


t<-system.time({
  out <- predict(m2, dtest)
})
results[7,ii,5]<-t[3]
out<-as.integer(out>=0.5)

print(results[7,ii,1]<-(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
results[7,ii,2]<-sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
results[7,ii,3]<-sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))

# xgboost logLik gbtree
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

  train<-as.data.frame(data.example[,-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)])
  test<-as.data.frame(data.example1[,-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)])

  dval<-xgb.DMatrix(data = data.matrix(train[,-1]), label = data.matrix(train[,1]),missing=NA)
  watchlist<-list(dval=dval)


  m2 <- xgb.train(data = xgb.DMatrix(data = data.matrix(train[,-1]), label = data.matrix(train[,1]),missing=NA),
                  param, nrounds = 10000,
                  watchlist = watchlist,
                  print_every_n = 10)

})
results[8,ii,4]<-t[3]
# Predict
t<-system.time({
  dtest  <- xgb.DMatrix(data.matrix(test[,-1]),missing=NA)
})

t<-system.time({
  out <- predict(m2, dtest)
})

results[8,ii,5]<-t[3]
out<-as.integer(out>=0.5)
print(results[8,ii,1]<-(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
results[8,ii,2]<-sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
results[8,ii,3]<-sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))


#GLMNET (elastic networks) # lasso a=1



t<-system.time({
  fit2 <- glmnet(as.matrix(train)[,-1], train$neo, family="binomial")
})
results[9,ii,4]<-t[3]

mmm<-as.matrix(test[,-1])
mmm[which(is.na(mmm))]<-0
t<-system.time({
  out <- predict(fit2,mmm , type = "response")[,fit2$dim[2]]
})
results[9,ii,5]<-t[3]

out<-as.integer(out>=0.5)

print(results[9,ii,1]<-(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
results[9,ii,2]<-sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
results[9,ii,3]<-sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))

# ridge a=0

t<-system.time({
  fit2 <- glmnet(as.matrix(train)[,-1], train$neo, family="binomial",alpha=0)
})
results[10,ii,4]<-t[3]

mmm<-as.matrix(test[,-1])
mmm[which(is.na(mmm))]<-0
t<-system.time({
  out <- predict(fit2,mmm , type = "response")[,fit2$dim[2]]
})

results[10,ii,5]<-t[3]

out<-as.integer(out>=0.5)

print(results[10,ii,1]<-(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
results[10,ii,2]<-sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
results[10,ii,3]<-sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))

gc()

# h2o.random forest



df <- as.h2o(train)



train1 <- h2o.assign(df , "train1.hex")
valid1 <- h2o.assign(df , "valid1.hex")
test1 <- h2o.assign(as.h2o(test[,-1]), "test1.hex")

train1[1:5,]

features = names(train1)[-1]

# in order to make the classification prediction
train1$neo <- as.factor(train1$neo)

t<-system.time({
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
results[11,ii,4]<-t[3]
t<-system.time({
  out<-h2o.predict(rf1,as.h2o(test1))[,1]
})
results[11,ii,5]<-t[3]
out<-as.data.frame(out)

out<-as.integer(as.numeric(as.character(out$predict)))


print(results[11,ii,1]<-(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
results[11,ii,2]<-sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
results[11,ii,3]<-sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))

#h2o deeplearning

t<-system.time({
  neo.dl <- h2o.deeplearning(x = features, y = "neo",hidden=c(200,200,200,200,200,200),
                             distribution = "bernoulli",
                             training_frame = train1,
                             validation_frame = valid1,
                             seed = 12341)
})
# now make a prediction

results[12,ii,4]<-t[3]
t<-system.time({
  out<-h2o.predict(neo.dl,as.h2o(test1))[,1]
})
results[12,ii,5]<-t[3]
out<-as.data.frame(out)

out<-as.integer(as.numeric(as.character(out$predict)))


print(results[12,ii,1]<-(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
results[12,ii,2]<-sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
results[12,ii,3]<-sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))


#h2o glm

t<-system.time({
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
results[13,ii,4]<-t[3]

t<-system.time({
  out<-h2o.predict(neo.glm,as.h2o(test1))[,1]
})
results[13,ii,5]<-t[3]
out<-as.data.frame(out)

out<-as.integer(as.numeric(as.character(out$predict)))


print(results[13,ii,1]<-(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
results[13,ii,2]<-sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
results[13,ii,3]<-sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))

#h2o naive bayes

t<-system.time({
  neo.nb <- h2o.naiveBayes(x = features, y = "neo",
                           training_frame = train1,
                           validation_frame = valid1,
                           seed = 12341)
})
# now make a prediction

results[14,ii,4]<-t[3]
t<-system.time({
  out<-h2o.predict(neo.nb,as.h2o(test1))[,1]
})
results[14,ii,5]<-t[3]
out<-as.data.frame(out)

out<-as.integer(as.numeric(as.character(out$predict)))


print(results[14,ii,1]<-(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
results[14,ii,2]<-sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
results[14,ii,3]<-sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))

#h2o kmeans

t<-system.time({
  neo.nb <- h2o.kmeans(x = c(features,"neo"),k=2,
                       training_frame = train1,
                       seed = 12341)
})
results[15,ii,4]<-t[3]

# now make a prediction


test2 <- h2o.assign(as.h2o(test), "test2.hex")

t<-system.time({
  out<-h2o.predict(neo.dl,as.h2o(test2))[,1]
})
results[15,ii,5]<-t[3]
out<-as.data.frame(out)

out<-as.integer(as.numeric(as.character(out$predict)))


print(results[15,ii,1]<-(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

#FNR
ps<-which(data.example1$neo==1)
results[15,ii,2]<-sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

#FPR
ns<-which(data.example1$neo==0)
results[15,ii,3]<-sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))

gc()
})), abort = function(){onerr<-TRUE;out<-NULL})})
}

ids<-NULL
for(i in 1:100)
{
  if(min(results[,i,1])>0)
    ids<-c(ids,i)

}

ress<-results[,ids[1:100],]

summary.results<-array(data = NA,dim = c(15,15))
for(i in 1:15)
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
rownames(summary.results)<-c("GMJMCMC(AIC)","MJMCMC(AIC)","GMJMCMC(LASSO)","MJMCMC(LASSO)","tXGBOOST(AUC)","tXGBOOST(logLik)","lXGBOOST(AUC)","lXGBOOST(logLik)","LASSO","RIDGE","RFOREST","DEEPNETS","NAIVEBAYESS","LR","KMEANS")


for(i in 1:15)
{
  plot(density(ress[i,,1],bw = "SJ"), main="Compare Kernel Density of precisions")
  polygon(density(ress[i,,1],bw = "SJ"), col="red", border="blue")

}

write.csv(x = round(summary.results,4),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/e1 p2/asteroids3.csv")

