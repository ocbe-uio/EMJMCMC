#read the appropriate version of the package
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package2.r")

#define some function for cleaning up after parallel computations
library(inline)
includes = '#include <sys/wait.h>'
code = 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait = cfunction(body=code, includes=includes, convention='.C')

#define the function estimating parameters of a given Gaussian logic regression with Jeffrey's prior
estimate.logic.lm = function(formula, data, n, m, r = 1)
{
  out = lm(formula = formula,data = data)
  p = out$rank
  fmla.proc=as.character(formula)[2:3]
  fobserved = fmla.proc[1]
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam =stri_split_fixed(str = fmla.proc[2],pattern = "+",omit_empty = F)[[1]]
  sj=(stri_count_fixed(str = fparam, pattern = "&"))
  sj=sj+(stri_count_fixed(str = fparam, pattern = "|"))
  sj=sj+1
  Jprior = prod(factorial(sj)/((m^sj)*2^(2*sj-2)))
  #tn=sum(stri_count_fixed(str = fmla.proc[2], pattern = "I("))
  if(length(which(is.na(coef(out))))>0)
    mlik =-10000
  else
    mlik = (-BIC(out)+2*log(Jprior) + 2*p*log(r)+n)/2
  if(mlik==-Inf)
    mlik = -10000
  return(list(mlik = mlik,waic = AIC(out)-n , dic =  BIC(out)-n,summary.fixed =list(mean = coef(out))))
}

#define the function for perfroming parallel computations
parall.gmj = mclapply

#define the function simplifying logical expressions at the end of the search
simplifyposteriors=function(X,posteriors,th=0.0001,thf=0.5)
{
  posteriors=posteriors[-which(posteriors[,2]<th),]
  rhash=hash()
  for(i in 1:length(posteriors[,1]))
  {
    expr=posteriors[i,1]
    print(expr)
    res=model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
    res[,1]=res[,1]-res[,2]
    ress=c(stri_flatten(res[,1],collapse = ""),stri_flatten(res[,2],collapse = ""),posteriors[i,2],expr)
    if(!(ress[1] %in% values(rhash)||(ress[2] %in% values(rhash))))
      rhash[[ress[1]]]=ress
    else
    {
      if(ress[1] %in% keys(rhash))
      {
        rhash[[ress[1]]][3]= (as.numeric(rhash[[ress[1]]][3]) + as.numeric(ress[3]))
        if(stri_length(rhash[[ress[1]]][4])>stri_length(expr))
          rhash[[ress[1]]][4]=expr
      }
      else
      {
        rhash[[ress[2]]][3]= (as.numeric(rhash[[ress[2]]][3]) + as.numeric(ress[3]))
        if(stri_length(rhash[[ress[2]]][4])>stri_length(expr))
          rhash[[ress[2]]][4]=expr
      }
    }
    
  }
  res=as.data.frame(t(values(rhash)[c(3,4),]))
  res$V1=as.numeric(as.character(res$V1))
  res=res[which(res$V1>thf),]
  res=res[order(res$V1, decreasing = T),]
  clear(rhash)
  rm(rhash)
  res[which(res[,1]>1),1]=1
  colnames(res)=c("posterior","tree")
  return(res)
}


#define number of simulations
MM = 100
#define number of threads to be used
M = 32
#define the size of the simulated samples
NM= 1000
#define \k_{max} + 1 from the paper
compmax = 41
#define treshold for preinclusion of the tree into the analysis
th=(10)^(-5)
#define a final treshold on the posterior marginal probability for reporting a tree  
thf=0.05

#define a function performing the map step for a given thread
runpar=function(vect)
{
  
  set.seed(as.integer(vect[22]))
  do.call(runemjmcmc, vect[1:21])
  vals=values(hashStat)
  fparam=mySearch$fparam
  cterm=max(vals[1,],na.rm = T)
  ppp=mySearch$post_proceed_results_hash(hashStat = hashStat)
  post.populi=sum(exp(values(hashStat)[1,][1:NM]-cterm),na.rm = T)
  clear(hashStat)
  rm(hashStat)
  rm(vals)
  return(list(post.populi = post.populi, p.post =  ppp$p.post, cterm = cterm, fparam = fparam))
}


#perform MM runs of GMJMCMC on M threads each
for(j in 1:MM)
{
  tryCatch({
    
    #prepare the data for simulation j
    set.seed(j)
    X4= as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))
    Y4=rnorm(n = 1000,mean = 1+7*(X4$V4*X4$V17*X4$V30*X4$V10)+7*(as.integer((X4$V50*X4$V19+X4$V13*X4$V11)>0)) + 9*(X4$V37*X4$V20*X4$V12)+ 7*(X4$V1*X4$V27*X4$V3)
             +3.5*(X4$V9*X4$V2) + 6.6*(X4$V21*X4$V18) + 1.5*X4$V7 + 1.5*X4$V8,sd = 1)
    X4$Y4=Y4
    
    #specify the initial formula
    formula1 = as.formula(paste(colnames(X4)[51],"~ 1 +",paste0(colnames(X4)[-c(51)],collapse = "+")))
    data.example = as.data.frame(X4)
    #specify tuning parameters of the algorithm for exploring the space of Bayesian logic regressions of interest
    #notice that allow_offsprings=1 corresponds to the GMJMCMC algorithm for Bayesian logic regression
    vect=list(formula = formula1,outgraphs=F,data = X4,estimator = estimate.logic.lm,estimator.args =  list(data = data.example,n = 1000, m = 50),recalc_margin = 249, save.beta = F,interact = T,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=1,mutation_rate = 250,last.mutation = 10000, max.tree.size = 4, Nvars.max =40,p.allow.replace=0.7,p.allow.tree=0.2,p.nor=0,p.and = 0.9),n.models = 20000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,print.freq = 1000,advanced.param = list(
      max.N.glob=as.integer(10),
      min.N.glob=as.integer(5),
      max.N=as.integer(3),
      min.N=as.integer(1),
      printable = F))
    
    params = list(vect)[rep(1,M)]
    #specify additional information to be used on different threads (e.g. seeds)
    for(i in 1:M)
    {
      params[[i]]$cpu=i
      params[[i]]$simul="scenario_3_"
      params[[i]]$simid=j
    }
    #perform garbage collection
    gc()
    
    #explore Byesian logic regression on M threads in parallel using the GMJMCMC algorithm
    print(paste0("begin simulation ",j))
    results=parall.gmj(X = params,FUN = runpar,mc.preschedule = T, mc.cores = M)
    #clean up
    gc()
    wait()
    
    #prepare the data structures for final analysis of the runs
    resa=array(data = 0,dim = c(compmax,M*3))
    post.popul = array(0,M)
    max.popul = array(0,M)
    nulls=NULL
    not.null=1
    #check which threads had non-zero exit status
    for(k in 1:M)
    {
      if(length(results[[k]])<=1||length(results[[k]])==0)
      {
        nulls=c(nulls,k)
        next
      }
      else
      {
        not.null = k
      }
      
    }
    
    #for all of the successful runs collect the results into the corresponding data structures
    for(k in 1:M)
    {
      if(k %in% nulls)
      {
        results[[k]]=results[[not.null]]
      }
      max.popul[k]=results[[k]]$cterm
      post.popul[k]=results[[k]]$post.populi
      resa[,k*3-2]=c(results[[k]]$fparam,"Post.Gen.Max")
      resa[,k*3-1]=c(results[[k]]$p.post,results[[k]]$cterm)
      resa[,k*3]=rep(post.popul[k],length(results[[k]]$p.post)+1)
      
    }
    #delete the unused further variables and perfrom garbage collection
    rm(results)
    gc()
    #renormalize estimates of the marginal inclusion probabilities
    #based on all of the runs
    ml.max=max(max.popul)
    post.popul=post.popul*exp(-ml.max+max.popul)
    p.gen.post=post.popul/sum(post.popul)
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
            if(as.integer(has.key(hash = hfinal,key =resa[jj,ii*3-2]))==0)
              hfinal[[resa[jj,ii*3-2]]]=as.numeric(resa[jj,ii*3])
            else
              hfinal[[resa[jj,ii*3-2]]]=hfinal[[resa[jj,ii*3-2]]]+as.numeric(resa[jj,ii*3])
          }
          
        }
      }
    }
    
    posteriors=values(hfinal)
    #delete the unused further variables
    clear(hfinal)
    rm(hfinal)
    rm(resa)
    rm(post.popul)
    rm(max.popul)
    #simplify the found trees and their posteriors
    posteriors=as.data.frame(posteriors)
    posteriors=data.frame(X=row.names(posteriors),x=posteriors$posteriors)
    posteriors$X=as.character(posteriors$X)
    tryCatch({
      res1=simplifyposteriors(X = X4,posteriors = posteriors, th,thf)
      write.csv(x =res1,row.names = F,file = paste0("post3etaOld_",j,".csv"))
    },error = function(err){
      print("error")
      write.csv(x =posteriors,row.names = F,file = paste0("posteriors3etaOld_",j,".csv"))
    },finally = {
      
      print(paste0("end simulation ",j))
      
    }) 
    #delete the unused further variables and perfrom garbage collection
    rm(X4)
    rm(data.example)
    rm(vect)
    rm(params)
    gc()
    print(paste0("end simulation ",j))
  },error = function(err){
    print("error")
    j=j-1
    print(paste0("repeat  simulation ",j))
  },finally = {
    
    print(paste0("end simulation ",j))
    
  }) 
  
}
