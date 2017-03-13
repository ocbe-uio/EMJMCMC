library(LogicReg)
library(scrime)
library(hash)
library(parallel)

M=100
resa<-array(data = 0,dim = c(16,M*3))
NM= 1000
post.popul <- array(0,M)
max.popul <- array(0,M)
set.seed(040590)
fblr()

y <- X1$Y1
# begin with the MCLR
system.time({
  lmod <- logreg(
    resp=y,bin=X1[,-51],
    type=3, select = 7,
    ntrees=5,
    nleaves =5*5,
    seed = 040590,
    mc.control=logreg.mc.control(nburn=5000, niter=100000, hyperpars=log(2),output = 3),
    tree.control=logreg.tree.control(treesize=2
                                     ,opers=2)
  )
})

lmod$size
iter<-sum(lmod$size[,2])
# marginal posteriors of covariates in over models
paste(which(lmod$single > 0.5*iter),paste("X"))
paste(which(lmod$double > 0.5*iter),paste("XX"))
#paste(which(lmod$triple > 0.5*iter),paste("XXX"))





library(LogicReg)
library(scrime)
library(hash)
library(parallel)
library(inline)
includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')


rhash<-hash()
clear(rhash)

runpar<-function(vect)
{
  set.seed(as.integer(vect))
  print(vect)
  
  X1<- as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = 0.3),dim = c(1000,50)))
  Y1=-0.7+1*((1-X1$V1)*(X1$V4)) + 1*(X1$V8*X1$V11)+1*(X1$V5*X1$V9)
  X1$Y1<-round(1.0/(1.0+exp(-Y1)))
  
  
  #setwd(dir = paste("mclr",vect))
  y <- X1$Y1
  # begin with the MCLR
  system.time({
    lmod <- logreg(
      resp=y,bin=X1[,-51],
      type=3, select = 7,
      ntrees=5,
      nleaves =5*5,
      seed = vect,
      mc.control=logreg.mc.control(nburn=5000, niter=100000, hyperpars=log(2),output = 3),
      tree.control=logreg.tree.control(treesize=2
                                       ,opers=2)
    )
  })
  
  
  lmod$size
  iter<-sum(lmod$size[,2])
  # marginal posteriors of covariates in over models
  n1<-paste(which(lmod$single > 0.5*iter),paste("X"))
  n2<-paste(which(lmod$double > 0.5*iter),paste("XX"))
  n3<-paste(which(lmod$triple > 0.5*iter),paste("XXX"))
  print("end")
  return(c(n1,n2,n3))
}

M<-100
N<-32
vect<-(1:M)

results<-mclapply(X = vect,FUN = runpar,mc.preschedule = T,mc.cores = 1)
wait()


for(names.all in results)
{
  for(name in names.all)
  {
    if(name %in% keys(rhash))
    {
      rhash[[name]]<-(as.integer(rhash[[name]])+1)
    }
    else
    {
      rhash[[name]]<-as.integer(1)
    }
    
  }
}

write.csv(file = "res1MCLR.csv",x = cbind(keys(rhash),values(rhash)))