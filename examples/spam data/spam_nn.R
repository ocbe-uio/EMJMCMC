#Neural network on spam data

library(MASS)
library(nnet)
library(class)
setwd("/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/spam data/")

spam = read.table("spam.data",col.names=c(paste("x",1:57,sep=""),"y"))
spam$y = as.factor(spam$y)
spam[,1:57] = scale(spam[,1:57])
spam.traintest = read.table("spam.traintest")
spam.train = spam[spam.traintest==1,]
spam.test = spam[spam.traintest==0,]

#Neural network with one hidden layer
#Note: The following call is VERY time-consuming
m = 10
spam.nnet = nnet(y~.,data=spam.train,size=m,decay=0.1,
                 MaxNWts=10000,maxit=300)
pred = predict(spam.nnet,spam.test,type="class")
1-mean(abs(as.numeric(pred)-as.numeric(as.character(spam.test$y))))
table(spam.test$y,pred)
err.nnet=mean(spam.test$y!=pred)

#Deep learning with three hidden layers
library(RSNNS)
spamTargets <- decodeClassLabels(spam$y)
spamTargets.train = spamTargets[spam.traintest==1,]
zipTargets.test = spamTargets[spam.traintest==0,]
maxit = c(500)
err.dnet = rep(NA,length(maxit))
lab = c(0,1)
for(i in 1:length(maxit))
{
  spam.dnet = mlp(as.matrix(spam.train[,1:57]),spamTargets.train, size = c(50),
                 learnFuncParams=c(0.3),maxit=maxit[i])
  pred.dnet = lab[apply(predict(spam.dnet,spam.test[,1:57]),1,which.is.max)]
  table(spam.test$y,pred.dnet)
  err.dnet[i] = mean(spam.test$y!=pred.dnet)
}

fit.glm = glm(y~.,data=spam.train,family=binomial())
stepAIC(fit.glm)
m = 10
spam.nnet = nnet(y~x1 + x5 + x6 + x7 + x8 + x9 + x15 + x16 + x20 +
                   x21 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + x31 + x33 +
                   x35 + x36 + x41 + x42 + x44 + x45 + x46 + x48 + x49 + x52 +
                   x53 + x55 + x57,data=spam.train,size=m,decay=0.1,
                 MaxNWts=10000,maxit=300)
pred = predict(spam.nnet,spam.test,type="class")
table(spam.test$y,pred)
err.nnet=mean(spam.test$y!=pred)
1-mean(abs(as.numeric(pred)-as.numeric(as.character(spam.test$y))))
library(devtools)
source_url('https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r')

#Make plot
par(mar=numeric(4),mfrow=c(1,1),family='serif')
plot.nnet(spam.dnet)
