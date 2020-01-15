#this script performs the experiment for the parallel runs of GMJMCMC/RGMJMCMC algoorithm on spam data.

#read in the package most recent version
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")

#specify the estimator function returning p(Y|m)p(m), model selection criteria and the vector of the modes for the beta coefficients

estimate.bas.glm.cpen = function(formula, data, family, prior, logn,r = 0.1,yid=1,relat =c("sigmoid","glquar","gmedi","gfquar","tanh","atan","sin"))
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

#define a hashtable, where the found features are stored
featgmj = hash()


#specify the function for parallel computations
parall.gmj <<- mclapply

setwd("/nr/samba/user/ahu/EMJMCMC2016/supplementaries/BGNLM/spam/")

#read in the train and test data sets
data = read.table("spam.data",col.names=c(paste("x",1:57,sep=""),"X"))
data[,1:57] = scale(data[,1:57])

set.seed(1)
spam.traintest =  read.table("spam.traintest")#rbinom(n = dim(data)[1],size = 1,prob = 0.01)
train = data[spam.traintest==1,]
test = data[spam.traintest==0,]

#transform the train data set to a data.example data.frame that EMJMCMC class will internally use
data.example = train

#perfrom garbage collection
gc()

#define an array for storing final results from the individual 100 runs
total = array(0,dim = c(100,10,3))

#specify the link function that will be used in the prediction phase
g=function(x)
{
  return((x = 1/(1+exp(-x))))
}

#a function that is run to perform analyses on each of the addressed threads
runpar=function(vect)
{
  #specify the seed
  set.seed(as.integer(vect$seed))
  #run emjmcmc on the given vector of tuning parameters
  do.call(runemjmcmc, vect[1:vect$runlen])

  #postproceed results and obtain marginal inclusiona probabilities and posterior marginal model proabilities
  #from the set of explored in a given run models
  ppp=mySearch$post_proceed_results_hash(hashStat = hashStat)
  ppp$p.post

  #get the modes of beta coefficients for the explored models
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

  #make predictions for the test data set
  t=system.time({

    res=mySearch$forecast.matrix.na(link.g = g, covariates = (vect$test),betas = betas,mliks.in = mliks)$forecast

  })

  #remove the non-required further variables
  rm(betas)
  rm(mliks)
  clear(hashStat)
  rm(hashStat)
  #return the variables of intererst for this example, including
  #marginal inclusion probabilities (p.post), the set of corresponding features (fparam)
  #prediction results (res), the sum of parginal likelihoods times the priors (post.populi)
  #maximal marginal likelihood times the prior in the set of explored models
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
compmax = 101

run =4

#specify the number of threads used in the runs
M=32
for(j in 1:100)
{

  set.seed(j)
  #specify the initial formula
  formula1 = as.formula(paste(colnames(data.example)[58],"~ 1 +",paste0(colnames(data.example)[-58],collapse = "+")))
  #specify tuning parameters of the algorithm for exploring BGNLM of interest
  #notice that allow_offsprings=3 corresponds to the GMJMCMC runs and
  #allow_offsprings=4 -to the RGMJMCMC runs
  vect= list(formula = formula1,data = data.example,gen.prob = c(1,1,1,1,0),estimator =estimate.bas.glm.cpen,estimator.args =  list(data = data.example,prior = aic.prior() ,family = binomial(),yid=58, logn = log(1536),r=exp(-0.5)),recalc_margin = 95, save.beta = T,interact = ifelse(run<9,T,F),relations = c("gauss","tanh","atan","sin","glquar","gmedi","gfquar"),relations.prob =c(0.1,0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=ifelse(run<5,3,4),mutation_rate = 100,last.mutation=1000, max.tree.size = 60, Nvars.max = 100,p.allow.replace=0.5,p.allow.tree=0.4,p.nor=0.3,p.and = 0.9),n.models = 20000,unique =F,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 100,advanced.param = list(
    max.N.glob=as.integer(20),
    min.N.glob=as.integer(5),
    max.N=as.integer(3),
    min.N=as.integer(1),
    printable = F))
  len = length(vect)

  params = list(vect)[rep(1,M)]

  #specify additional information to be used on different threads (e.g. seeds)
  for(jj in 1:M)
  {
    params[[jj]]$runlen = len
    params[[jj]]$cpu=jj
    params[[jj]]$seed = jj*j
    params[[jj]]$test=test
  }
  #perform garbage collection
  gc()

  length(params[[1]])
  #explore BGNLM on M threads in parallel using the GMJMCMC/RGMJMCMC algorithm
  results=parall.gmj(X = params,FUN = runpar,mc.preschedule = F, mc.cores = M,mc.cleanup = T)

  #prepare the data structures for final analysis of the runs
  resa=array(data = 0,dim = c(compmax,M*3))
  post.popul = array(0,M)
  max.popul = array(0,M)
  nulls=NULL

  #check which threads had non-zero exit status
  for(k in 1:M)
  {
    if(length(results[[k]])==1||length(results[[k]]$cterm)==0)
    {
      nulls=c(nulls,k)
      next
    }
    else
    {
      not.null = k
    }

  }

  #for all of the successful runs collect the results into the data structures
  for(k in 1:M)
  {
    if(k %in% nulls)
    {
      results[[k]]=results[[not.null]]
    }
    max.popul[k]=results[[k]]$cterm
    post.popul[k]=results[[k]]$post.populi
    if(length(resa[,k*3-2])==(length(results[[k]]$fparam)+1))
    {
      resa[,k*3-2]=c(results[[k]]$fparam,"Post.Gen.Max")
      resa[,k*3-1]=c(results[[k]]$p.post,results[[k]]$cterm)
      resa[,k*3]=rep(post.popul[k],length(results[[k]]$p.post)+1)
    }else
    {
      resa[,k*3-2]=rep(results[[k]]$fparam[1],length(resa[,k*3-2]))
      resa[,k*3-1]=rep(0,length(resa[,k*3-1]))
      resa[,k*3]=rep(-10^9,length(resa[,k*3]))
    }

  }

  #renormalize estimates of the marginal inclusion probabilities
  #based on all of the runs
  ml.max=max(max.popul)
  post.popul=post.popul*exp(-ml.max+max.popul)
  p.gen.post=post.popul/sum(post.popul)

  #perform BMA of the redictions across the runs
  res1 = results[[1]]$res*p.gen.post[1]
  for(i in 2:M)
  {

    res1=res1+results[[i]]$res*p.gen.post[i]

  }
  #calculate the performance metrics
  for(jjjj in 1:10)
  {
    res = as.integer(res1>=0.1*jjjj)
    prec=(1-sum(abs(res-test$X),na.rm = T)/length(res))

    #FNR
    ps=which(test$X==1)
    fnr=sum(abs(res[ps]-test$X[ps]))/(sum(abs(res[ps]-test$X[ps]))+length(ps))

    #FPR
    ns=which(test$X==0)
    fpr=sum(abs(res[ns]-test$X[ns]))/(sum(abs(res[ns]-test$X[ns]))+length(ns))

    total[j,jjjj,1]=prec
    total[j,jjjj,2]=fnr
    total[j,jjjj,3]=fpr



  }

  #obtain posteriors for the explored features across all of the threads in a given run
  hfinal=hash()
  for(ii in 1:M)
  {
    resa[,ii*3]=p.gen.post[ii]*as.numeric(resa[,ii*3-1])
    resa[length(resa[,ii*3]),ii*3]=p.gen.post[ii]
    if(p.gen.post[ii]>0)
    {
      for(jj in 1:(length(resa[,ii*3])-1))
      {
        if(resa[jj,ii*3]>0)
        {
          #print(paste0(ii,"  and ",jj))
          if(as.integer(has.key(hash = hfinal,key =resa[jj,ii*3-2]))==0)
            hfinal[[resa[jj,ii*3-2]]]=as.numeric(resa[jj,ii*3])
          else
            hfinal[[resa[jj,ii*3-2]]]=hfinal[[resa[jj,ii*3-2]]]+as.numeric(resa[jj,ii*3])
        }

      }
    }
  }

  #write the features from a given run into the file
  posteriors=values(hfinal)
  write.csv(x =posteriors,row.names = F,file = paste0("parres/posteriorsspam_",j,".csv"))
  print(posteriors)
  clear(hfinal)
  rm(hfinal)
  rm(resa)
  rm(post.popul)
  rm(max.popul)
  posteriors=as.data.frame(posteriors)
  posteriors=data.frame(X=row.names(posteriors),x=posteriors$posteriors)
  posteriors$X=as.character(posteriors$X)
  write.csv(x =posteriors,row.names = F,file = paste0("parres/posteriorsspam_",j,".csv"))
  print(paste0("end simulation ",j))

  #print the run's metrics and clean the results


  print(max(total[j,,1]))

  write.csv(file =paste0("parres/resultsrun_",j,".csv"),x=total[j,5,])


  rm(results)
  gc()
}

#check all of the files
results=array(NA,dim = c(2,100,5))
for(j in 1:100)
{
  if(!file.exists(paste0("parres/resultsrun_",j,".csv")))
    next
  tmp = read.csv(paste0("parres/resultsrun_",j,".csv"))
  
  results[1,j,1]=  tmp[1,2]
  results[1,j,2]=  tmp[2,2]
  results[1,j,3] =   tmp[3,2]
  
  
}
#make the joint summary of the runs, including min, max and medians of the performance metrics
summary.results=array(data = NA,dim = c(1,15))


for(j in 1:5)
{
  summary.results[1,(j-1)*3+1]=min(results[1,,j],na.rm = T)
  summary.results[1,(j-1)*3+2]=median(results[1,,j],na.rm = T)
  summary.results[1,(j-1)*3+3]=max(results[1,,j],na.rm = T)
}
summary.results=as.data.frame(summary.results)

names(summary.results)=c("median(prec)","min(prec)","max(prec)","median(fnr)","min(fnr)","max(fnr)","median(fpr)","min(fpr)","max(fpr)")
rownames(summary.results)[1]=paste0("BGNLM",run)


#write the final results into the files
write.csv(file =paste0("parres/summary1",run,".csv"),x=paste0(summary.results[1,1]," (",summary.results[1,2],",",summary.results[1,3],")&",summary.results[1,4]," (",summary.results[1,5],",",summary.results[1,6],")&",summary.results[1,7]," (",summary.results[1,8],",",summary.results[1,9],")"),row.names = F,col.names = F)
write.csv(file ="parres/results.csv",x=total)
write.csv(x = cbind(keys(featgmj),values(featgmj)),file = paste0("parres/spamfeatgmj",1,".csv"))
