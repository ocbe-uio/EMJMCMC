#this script performs the experiment for the single runs of GMJMCMC/RGMJMCMC algoorithm on breast cancer data.

#read in the package most recent version
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")


#***********************IMPORTANT******************************************************
# if a multithreaded backend openBLAS for matrix multiplications
# is installed on your machine, please force it to use 1 thread explicitly
# you might want to experiment with the combinations of blas_set_num_threads and ncores
library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)
#***********************IMPORTANT******************************************************
#specify the estimator function returning p(Y|m)p(m), model selection criteria and the vector of the modes for the beta coefficients

#specify the estimator function returning p(Y|m)p(m), model selection criteria and the vector of the modes for the beta coefficients
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

#read in the train and test data sets
test = read.csv("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/BGNLM/breast%20cancer/test.csv",header = T,sep=",")[,-1]
train = read.csv("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/BGNLM/breast%20cancer/train.csv",header = T,sep=",")[,-1]


#transform and attach the train data set to a data.example data.frame that EMJMCMC class will internally use
data.example = as.data.frame(train)

attach(data.example)
set.seed(10)

#specify the function for parallel computations
parall.gmj <<- mclapply


#perfrom garbage collection
gc()



#specify the link function that will be used in the prediction phase
g=function(x)
{
  return((x = 1/(1+exp(-x))))
}

#a function that is run to perform analyses on each of the addressed threads
runpar=function(vect)
{

  set.seed(as.integer(vect$seed))
  do.call(runemjmcmc, vect[1:vect$runlen])

  ppp=mySearch$post_proceed_results_hash(hashStat = hashStat)
  ppp$p.post


  Nvars=mySearch$Nvars
  linx =mySearch$Nvars+4
  lHash=length(hashStat)
  mliks = values(hashStat)[which((1:(lHash * linx)) %% linx == 1)]
  betas = values(hashStat)[which((1:(lHash * linx)) %% linx == 4)]
  cterm=max(values(hashStat)[1,],na.rm = T)
  post.populi=sum(exp(values(hashStat)[1,][1:1000]-cterm),na.rm = T)
  for(i in 1:(Nvars-1))
  {
    betas=cbind(betas,values(hashStat)[which((1:(lHash * linx)) %% linx == (4+i))])
  }
  betas=cbind(betas,values(hashStat)[which((1:(lHash * linx)) %% linx == (0))])

  t=system.time({

    res=mySearch$forecast.matrix.na(link.g = g, covariates = (vect$test),betas = betas,mliks.in = mliks)$forecast

  })

  rm(betas)
  rm(mliks)
  clear(hashStat)
  rm(hashStat)
  return(list(p.post =  ppp$p.post, fparam = mySearch$fparam, res = res, post.populi = post.populi, cterm = cterm))
}


#perform garbage collection
gc()
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

M=32
max.num = 128

#specify data structures for the results
result = array(0,dim = c(max.num,4))
featgmj = hash()

#run the algorithms
for(run in 1:9)
{


  #specify the initial formula
  set.seed(run)
  formula1 = as.formula(paste(colnames(data.example)[31],"~ 1 +",paste0(colnames(data.example)[-31],collapse = "+")))


  #specify tuning parameters of the algorithm for exploring BGNLM of interest
  #notice that allow_offsprings=3 corresponds to the GMJMCMC runs and
  #allow_offsprings=4 -to the RGMJMCMC runs
  vect= list(formula = formula1,data = data.example,gen.prob = c(1,1,1,1,0),estimator =estimate.bas.glm.cpen,estimator.args =  list(data = data.example,prior = aic.prior(),link = "sigmoid", distribution = "binomial" ,family = binomial(),yid=31, logn = log(143),r=exp(-0.5)),recalc_margin = 95, save.beta = T, deep.method = run%%4,interact = ifelse(run<9,T,F),relations = c("gauss","tanh","atan","sin"),relations.prob =c(0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=ifelse(run<5,3,4),mutation_rate = 100,last.mutation=1000, max.tree.size = 6, Nvars.max = 20,p.allow.replace=0.5,p.allow.tree=0.4,p.nor=0.3,p.and = 0.9),n.models = 7000,unique =T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
    max.N.glob=as.integer(10),
    min.N.glob=as.integer(5),
    max.N=as.integer(3),
    min.N=as.integer(1),
    printable = F))


  len = length(vect)

  params = list(vect)[rep(1,max.num)]

  #specify additional information to be used on different threads (e.g. seeds)
  #n/b here we will run 100 runs of the single threaded algorithms in parallel
  #in order to minimize the time untill the results are all ready
  for(jj in 1:max.num)
  {
    params[[jj]]$runlen = len
    params[[jj]]$cpu=jj
    params[[jj]]$seed = jj*run
    params[[jj]]$test=test
  }
  #perform garbage collection
  gc()

  #run the threads with the specified parameters on different seeds
  results=parall.gmj(X = params,FUN = runpar,mc.preschedule = F, mc.cores = M,mc.cleanup = T)

  #check which threads had non-zero exit status
  nulls = NULL
  not.nulls =  NULL
  for(k in 1:max.num)
  {
    if(length(results[[k]])==1||length(results[[k]]$cterm)==0)
    {
      nulls=c(nulls,k)
      next
    }
    else
    {
      not.nulls = c(not.nulls,k)
    }

  }

  #get 100 non.null results to proceed
  #normally there should not be any nulls in general
  if(length(not.nulls)>100)
    not.nulls = not.nulls[1:100]


  #perform classifications for the individuals runs
  #and calculate corresponding performance metrics
  for(i in not.nulls)
  {

    res = as.integer(results[[i]]$res>0.5)
    prec=(1-sum(abs(res-test$X),na.rm = T)/length(res))

    #FNR
    ps = which(test$X==1)
    fnr = sum(abs(res[ps]-test$X[ps]))/(sum(abs(res[ps]-test$X[ps]))+length(ps))

    #FPR
    ns = which(test$X==0)
    fpr = sum(abs(res[ns]-test$X[ns]))/(sum(abs(res[ns]-test$X[ns]))+length(ns))

    crit=results[[i]]$cterm

    result[i,1]=prec
    result[i,2]=fnr
    result[i,3]=fpr
    result[i,4]=crit

    #summarize the features found in all of the runs
    for(j in which(results[[i]]$p.post>0.1))
    {  if(!has.key(hash = featgmj,key =  results[[i]]$fparam[j]))
      featgmj[[results[[i]]$fparam[j]]] = as.numeric(1) else{
        featgmj[[results[[i]]$fparam[j]]] =as.numeric(featgmj[[results[[i]]$fparam[j]]]) + 1
      }
    }

  }

  #write results down
  write.csv(file = paste0("results",run,".csv"),x=result)

  ids = NULL
  for(i in 1:100)
  {
    if(result[i,1]>0)
      ids=c(ids,i)

  }



  #make the joint summary of the runs, including min, max and medians of the performance metrics
  ress=result[ids,]
  summary.results=array(data = NA,dim = c(1,12))
  for(j in 1:4)
  {
    summary.results[1,(j-1)*3+2]=min(ress[,j],na.rm = T)
    summary.results[1,(j-1)*3+1]=median(ress[,j],na.rm = T)
    summary.results[1,(j-1)*3+3]=max(ress[,j],na.rm = T)
  }
  summary.results=format(round(summary.results, digits = 4),4)
  summary.results=as.data.frame(summary.results)

  names(summary.results)=c("median(prec)","min(prec)","max(prec)","median(fnr)","min(fnr)","max(fnr)","median(fpr)","min(fpr)","max(fpr)","median(crit)","min(crit)","max(crit)")
  rownames(summary.results)[1]=paste0("BGNLM",run)

  #write the final reults into the files
  write.csv(file =paste0("summary1",run,".csv"),x=paste0(summary.results[1,1]," (",summary.results[1,2],",",summary.results[1,3],")&",summary.results[1,4]," (",summary.results[1,5],",",summary.results[1,6],")&",summary.results[1,7]," (",summary.results[1,8],",",summary.results[1,9],")&",summary.results[1,10]," (",summary.results[1,11],",",summary.results[1,12],")"),row.names = F,col.names = F)
  write.csv(x = cbind(keys(featgmj),values(featgmj)),file = paste0("spamfeatgmj",run,".csv"))

  #perorm garbage collection
  clear(featgmj)
  gc()
}


