library(LogicReg)
library(scrime)

M=100
resa<-array(data = 0,dim = c(16,M*3))
NM= 1000
post.popul <- array(0,M)
max.popul <- array(0,M)
set.seed(040590)

X1<- as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))
Y1<-rnorm(n = 1000,mean = 1+0.7*(X1$V1*X1$V4) + 0.8896846*(X1$V8*X1$V11)+1.434573*(X1$V5*X1$V9),sd = 1)
X1$Y1<-round(1.0/(1.0+exp(-Y1)))


y <- X1$Y1
# begin with the MCLR
system.time({
  lmod <- logreg(
    resp=y,bin=X1[,-51],
    type=3, select = 7,
    ntrees=5,
    nleaves =5*5,
    seed = 040590,
    mc.control=logreg.mc.control(nburn=500, niter=100000, hyperpars=log(2),output = 3),
    tree.control=logreg.tree.control(treesize=5
                                     ,opers=2)
  )
})

lmod$size
iter<-sum(lmod$size[,2])
# marginal posteriors of covariates in over models
lmod$single[1]/iter
lmod$single[4]/iter
lmod$single[8]/iter
lmod$single[11]/iter
lmod$single[5]/iter
lmod$single[9]/iter
max(lmod$single)/iter
min(lmod$single)/iter
which(lmod$single == max(lmod$single))
order(lmod$single)
# pairwise probability of the pairs of covariates in a tree over models
lmod$double[4,1]/iter
lmod$double[11,8]/iter
lmod$double[9,5]/iter
max(lmod$double)/iter
which(lmod$double == max(lmod$double))
#lmod$triple the best tupple of them is
max(lmod$triple)/iter

#clearly does not work for this simple data

#proceed with FBLR
snp <- as.matrix(X1[,1:50])
bin <- snp2bin(snp)
int <- apply(bin,1,function(x) (x[1] == 1 & x[3] == 0)*1)
y <- as.integer(X1$Y1)
# normally more iterations should be used
system.time({
fblr(y, bin, niter=100000, nburn=1000, int.level = 5, kmax = 15, geo = 1)
})

bbb=analyse.models("fblr_mcmc.txt",size.freq = F,moco = c(20,20,20),int.freq = F,kmax = 15, int.level = 5,bin.names = names(X1)[1:50])


