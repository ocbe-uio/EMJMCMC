#compare on Asteroids
library(RCurl)
library(caret)

library(inline)
includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')

source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package2.r")


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


parall.gmj <<- mclapply

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

xxx<-array()
digits[,2:785][which(digits[,2:785]<10)]<-0

index <- createDataPartition(digits$label, p = 0.05, list = FALSE)

test <- digits[-index, ]
train <- digits[index, ]

data.example <- as.data.frame(train,stringsAsFactors = T)

digits<-digits[,-which(colSums(digits)==0)]

which(colSums(digits)==0)

g<-function(x)
{
  return((x = 1/(1+exp(-x))))
}

data<-runif(1,1,42000)


index <- createDataPartition(digits$label, p = 0.005, list = FALSE)

test <- digits[-index, ]
train <- digits[index, ]

data.example <- as.data.frame(train,stringsAsFactors = T)


runpar<-function(vect)
{
  
  set.seed(as.integer(vect[23]))
  do.call(runemjmcmc, vect[1:22])
  
  ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
  ppp$p.post
  
 
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
    
    res<-mySearch$forecast.matrix.na(link.g = g,covariates = (vect$test),betas = betas,mliks.in = mliks)$forecast
    
  })
  
  rm(betas)
  rm(mliks)
  clear(hashStat)
  rm(hashStat)
  

  return(list(p.post =  ppp$p.post, fparam = mySearch$fparam, res = res))
}



gc()


gfquar<-function(x)as.integer(x<quantile(x,probs = 0.25,na.rm = T))
glquar<-function(x)as.integer(x>quantile(x,probs = 0.75,na.rm = T))
gmedi<-function(x)as.integer(x>median(x))
cosi<-function(x)cos(x/180*pi)
gmean<-function(x)as.integer(x>mean(x))

#idss<-which(abs(cor(x = train[2:785],y=train[785+dig]))>0.3)

M<-10

formula1 = "x"

vect<-list(formula = formula1,data = data.example,estimator =estimate.bas.glm.cpen,estimator.args =  list(data = data.example,prior = aic.prior(),family = binomial(), logn = log(64),r=exp(-0.5),yid=785),recalc_margin = 95,locstop=T,presearch=T,save.beta = T,interact = T,relations = c("","cosi","sigmoid","tanh","atan","erf","gmean","gmedi","gfquar","glquar"),relations.prob =c(0.8,0.1,0.1,0.1,0.1,0.1,0.1,0.5,0.1,0.1),interact.param=list(allow_offsprings=2,mutation_rate = 100,last.mutation=2500, max.tree.size = 20, Nvars.max =50,p.allow.replace=0.1,p.allow.tree=0.5,p.nor=0.3,p.and = 0.7),n.models = 30000,unique =F,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
    max.N.glob=as.integer(10),
    min.N.glob=as.integer(5),
    max.N=as.integer(3),
    min.N=as.integer(1),
    printable = F))

length(vect)

params <- list(vect)[rep(1,M)]


for(dig in 1:M)
{
  params[[dig]]$cpu<-dig
  idss<-which(abs(cor(x = train[2:785],y=train[785+dig]))>0.2)
  params[[dig]]$formula<-as.formula(paste(colnames(train)[785+dig],"~ 1 +",paste0(colnames(train)[idss],collapse = "+")))
  params[[dig]]$estimator.args <- list(data = data.example,prior = aic.prior(),family = binomial(), logn = log(64),r=exp(-1),yid=785+dig)
  params[[dig]]$test<-test
}
gc()

length(params[[1]])

results<-parall.gmj(X = params,FUN = runpar,mc.preschedule = T, mc.cores = M)
#print(results)
wait()

cbind(results[[1]]$fparam,results[[1]]$p.post)

for(i in 1:M)
{
  
  res<-results[[i]]$res
  results[[i]]$prec<-(1-sum(abs(res-test[,785+i]),na.rm = T)/length(res))
  
  #FNR
  ps<-which(test[,785+i]==1)
  results[[i]]$fnr<-sum(abs(res[ps]-test[,785+i][ps]))/(sum(abs(res[ps]-test[,785+i][ps]))+length(ps))
  
  #FPR
  ns<-which(test[,785+i]==0)
  results[[i]]$fpr<-sum(abs(res[ns]-test[,785+i][ns]))/(sum(abs(res[ns]-test[,785+i][ns]))+length(ns))
  
  gc()
  print(results[[i]]$prec)
}

which(is.na(digits))

