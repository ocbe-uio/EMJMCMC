#read the most recent stable version of the package
#source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")
#install.packages("https://github.com/aliaksah/EMJMCMC2016/blob/master/EMJMCMC_1.4.2_R_x86_64-pc-linux-gnu.tar.gz?raw=true", 
#                 repos = NULL, type="source")
# load the EMJMCMC package
library(EMJMCMC)
set.seed(040590)
#make sure that you are using Mac Os or Linux (mclapply is currently not supported for Windows unless some mclapply hack function for Windows is preloaded in your R session)



## Construct a binary correlation matrix for M = 50 cariables
M = 50
m = clusterGeneration::rcorrmatrix(M,alphad=2.5) 
#print the highest 10 correlations in the data
print(unique(sort(abs(m),decreasing = T))[1:10])
#print the lowest 10 correlations in the data
print(unique(sort(abs(m),decreasing = F))[1:10])
#simulate 1000 binary variables with a given correlation's structure
X = bindata::rmvbin(1000, margprob = rep(0.5,M), bincorr = m)


Xp = data.frame(X)

Xp$temper = rpois(1000,lambda = 34)
Y=rnorm(n = 1000,mean = 1+0.7*(Xp$X1*Xp$X4) + 0.8896846*(Xp$X8*Xp$X11)+1.434573*(Xp$X5*Xp$X9) + 0.01*Xp$temper,sd = 1)
Xp$Y=Y


#library(aRxiv)
set.seed(040590)
et = arima.sim(list(order=c(1,0,0), ar=0.8), n=1000)
Xp$Y = Xp$Y + et



#define the function estimating parameters of a given Gaussian logic regression with Jeffrey's prior
estimate.logic.inla = function(formula,family = "gaussian", latent = "+f(pos,model=\"ar1\")", data, n = 1000, m = 50 , r = 1,k.max=15, bias = -5000)
{
  
  
  print(formula)
  if(is.null(formula))
  { 
    print("formula missing")
    return(list(mlik =  -10000 + rnorm(1,0,1),waic =10000+ rnorm(1,0,1) , dic =  10000+ rnorm(1,0,1),summary.fixed =list(mean = 1)))
  }
  
  fmla.proc=as.character(formula)[2:3]
  fobserved = fmla.proc[1]
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam =stri_split_fixed(str = fmla.proc[2],pattern = "+",omit_empty = F)[[1]]
  
  p = length(fparam)
  if(p>k.max)
  {
    return(list(mlik = -10000+ rnorm(1,0,1),waic = 10000+ rnorm(1,0,1) , dic =  10000+ rnorm(1,0,1),summary.fixed =list(mean = 0)))
  }
  
  out<-NULL
  capture.output({tryCatch(capture.output({
    
    formula = as.formula(paste0(fmla.proc[1],"~",fmla.proc[2],latent))
    
    out <-inla(family = family,silent = 2L,data = data,formula = formula,control.compute = list(dic = TRUE, waic = TRUE, mlik = TRUE))
  }))})
  
  # use dic and aic as bic and aic correspondinly
  
  coef<-out$summary.fixed$mode
  #coef[1]<-coef[1]+out$summary.hyperpar$mode[1]
  
  sj=(stri_count_fixed(str = fparam, pattern = "&"))
  sj=sj+(stri_count_fixed(str = fparam, pattern = "|"))
  sj=sj+1
  Jprior = prod(factorial(sj)/((m^sj)*2^(2*sj-2)))
  
  if(is.null(out))
    return(list(mlik = -10000+log(r)*(sj),waic =  10000 , dic = 10000, summary.fixed =list(mean = NULL)))
  
  if(length(out$waic[1]$waic)==0)
    return(list(mlik = -10000+log(r)*(sj),waic =  10000 , dic = 10000, summary.fixed =list(mean = NULL)))
  
  mlik = (out$mlik[1]+log(Jprior) + p*log(r)+n) + bias
  if(mlik==-Inf||is.na(mlik))
    mlik = -10000 + rnorm(1,0,1)
  if(mlik==Inf)
    mlik = 10000 + rnorm(1,0,1)
  return(list(mlik = mlik,waic =  out$waic[1]$waic , dic = out$dic[1]$dic, summary.fixed =list(mean = coef)))
}

#estimate.logic.inla(data = Xp,formula = Y ~ temper+ 1 + X1 + X2 + X3 + X4 + X5)

#estimate.logic.inla(data = Xp,formula = Y ~ I(temper) + 1 + I(X1&X4) + I(X8&X11) + I(X5&X9),latent = "+f(pos,model=\"ar1\")")

Xp$pos = 1:length(Xp$X1)
#Xp$pos1 = 1:length(Xp$X1)
#Xp$pos2 = 1:length(Xp$X1)
#Xp$pos3 = 1:length(Xp$X1)

#estimate.logic.inla(data = Xp,formula = Y ~ temp+ 1 + X1 + X2 + X3 + X4 + X5 + f(Xp$pos,model="ar1") + f(data.example$pos1,model="rw1")+f(data.example$pos2,model="iid"))
# specify the initial formula
formula1 = as.formula(paste(colnames(Xp)[52],"~ 1 +",paste0(colnames(Xp)[-c(51,52,53)],collapse = "+")))
# specify the link function
g = function(x) x
data.example = Xp
res.mixed = pinferunemjmcmc(n.cores = 1, report.level =  0.2, num.mod.best = 1000,simplify = T,predict = F,test.data = NULL,link.function = g,runemjmcmc.params = list(formula = formula1,latnames = "",data = data.example,estimator = estimate.logic.inla,estimator.args = list(data = data.example, n = 1000, r = 1, m =stri_count_fixed(as.character(formula1)[3],"+"),k.max =10),recalc_margin = 249, save.beta = F,interact = T,outgraphs=F,relations=c("sin","cos"),relations.prob =c(0.5,0.5),interact.param=list(allow_offsprings=1,mutation_rate = 250,last.mutation=10000, max.tree.size = 3, Nvars.max =15,p.allow.replace=0.9,p.allow.tree=0.01,p.nor=0.01,p.and = 0.9),n.models = 10000,unique =T,max.cpu = 4,max.cpu.glob = 4,create.table = F,presearch = T, create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
  max.N.glob=as.integer(10),
  min.N.glob=as.integer(5),
  max.N=as.integer(3),
  min.N=as.integer(1),
  printable = T)))

res.mixed$threads.stats[[1]]$fparam
res.mixed$threads.stats[[1]]$p.post

res.mixed$allposteriors

