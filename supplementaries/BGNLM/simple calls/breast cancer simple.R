#inference
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")

#***********************IMPORTANT******************************************************
# if a multithreaded backend openBLAS for matrix multiplications
# is installed on your machine, please force it to use 1 thread explicitly
library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)
#***********************IMPORTANT******************************************************


#prediction
estimate.bas.glm.cpen = function(formula, data, link, distribution, family, prior, logn,r = 0.1,yid=1,relat=c("gauss","tanh","atan","sin"))
{
  capture.output({out = glm(family = family,formula = formula,data = data)})
  fmla.proc=as.character(formula)[2:3]
  fobserved = fmla.proc[1]
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  sj=2*(stri_count_fixed(str = fmla.proc[2], pattern = "*"))
  sj=sj+1*(stri_count_fixed(str = fmla.proc[2], pattern = "+"))
  for(rel in relat)
    sj=sj+2*(stri_count_fixed(str = fmla.proc[2], pattern = rel))
  mlik = ((-out$deviance +2*log(r)*sum(sj)))/2
  return(list(mlik = mlik,waic = -(out$deviance + 2*out$rank) , dic =  -(out$deviance + logn*out$rank),summary.fixed =list(mean = coefficients(out))))
}

#load some of the possible nonlinearities of interest
troot=function(x)abs(x)^(1/3)
sini=function(x)sin(x/180*pi)
logi=function(x)log(abs(x+0.1))
gfquar=function(x)as.integer(x<quantile(x,probs = 0.25))
glquar=function(x)as.integer(x>quantile(x,probs = 0.75))
gmedi=function(x)as.integer(x>median(x))
cosi=function(x)cos(x/180*pi)
gmean=function(x)as.integer(x>mean(x))
gone=function(x)as.integer(x>0)
gthird=function(x)(abs(x)^(1/3))
gfifth=function(x)(abs(x)^(1/5))
grelu=function(x)(x*(x>0))
contrelu=function(x)log(1+exp(x))
gauss=function(x)exp(-x*x)
compmax = 21

#read in the train and test data sets
test = read.csv("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/DBRM%20supplementaries/breast%20cancer/test.csv",header = T,sep=",")[,-1]
train = read.csv("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/DBRM%20supplementaries/breast%20cancer/train.csv",header = T,sep=",")[,-1]

#transform the train data set to a data.example data.frame that EMJMCMC class will internally use
data.example = as.data.frame(train)

#specify the link function that will be used in the prediction phase
g=function(x)
{
  return((x = 1/(1+exp(-x))))
}

formula1 = as.formula(paste(colnames(data.example)[31],"~ 1 +",paste0(colnames(data.example)[-31],collapse = "+")))

#define the number or cpus
M = 32

res = pinferunemjmcmc(n.cores = M, report.level =  0.5 , num.mod.best = NM,simplify = T, predict = T,test.data = as.data.frame(test),link.function = g, runemjmcmc.params =list(formula = formula1,data = data.example,gen.prob = c(1,1,1,1,0),estimator =estimate.bas.glm.cpen,estimator.args =  list(data = data.example,prior = aic.prior(),family = binomial(),yid=31, logn = log(143),r=exp(-0.5)),recalc_margin = 95, save.beta = T,interact = T,relations = c("gauss","tanh","atan","sin"),relations.prob =c(0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=4,mutation_rate = 100,last.mutation=1000, max.tree.size = 6, Nvars.max = 20,p.allow.replace=0.5,p.allow.tree=0.4,p.nor=0.3,p.and = 0.9),n.models = 7000,unique =T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
  max.N.glob=as.integer(10),
  min.N.glob=as.integer(5),
  max.N=as.integer(3),
  min.N=as.integer(1),
  printable = F)))

print(pROC::auc(predictor = res$predictions,response = test$X))

for(i in 1:M){
  print(pROC::auc(predictor =  res$threads.stats[[i]]$preds,response = test$X))
  print( res$threads.stats[[i]]$post.populi)
}

for(jjjj in 1:10)
{
  resw = as.integer(res$predictions>=0.1*jjjj)
  prec=(1-sum(abs(resw-test$X),na.rm = T)/length(resw))
  print(prec)
  #FNR
  ps=which(test$X==1)
  fnr=sum(abs(resw[ps]-test$X[ps]))/(sum(abs(resw[ps]-test$X[ps]))+length(ps))
  
  #FPR
  ns=which(test$X==0)
  fpr=sum(abs(resw[ns]-test$X[ns]))/(sum(abs(resw[ns]-test$X[ns]))+length(ns))
  
}