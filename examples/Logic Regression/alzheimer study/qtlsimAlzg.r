source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package2.r")

library(inline)
includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')

parall.gmj <<- mclapply

estimate.logic.glm <- function(formula, data, family =  binomial(), n, m, r = 1)
{
  X <- model.matrix(object = formula,data = data)
  out <- bayesglm.fit(x = X, y = data[,1], family=family,coefprior=aic.prior())
  p <- out$rank
  fmla.proc<-as.character(formula)[2:3]
  fobserved <- fmla.proc[1]
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam <-stri_split_fixed(str = fmla.proc[2],pattern = "+",omit_empty = F)[[1]]
  sj<-(stri_count_fixed(str = fparam, pattern = "&"))
  sj<-sj+(stri_count_fixed(str = fparam, pattern = "|"))
  sj<-sj+1
  Jprior <- sum(log(factorial(sj)/((m^sj)*2^(2*sj-2))))
  mlik = (-(out$deviance + log(n)*(out$rank)) + 2*(Jprior))/2+n
  if(mlik==-Inf)
    mlik = -10000
  return(list(mlik = mlik,waic = -(out$deviance + 2*out$rank) , dic =  -(out$deviance + log(n)*out$rank),summary.fixed =list(mean = coefficients(out))))
}



MM = 1
M = 64
NM= 1000
compmax = 26
th<-(10)^(-5)
thf<-0.05



paral<-function(X,FUN)
{
  return(mclapply(X = X,FUN = FUN,mc.preschedule = T, mc.cores = 32))
}

runpar<-function(vect)
{
  
  set.seed(as.integer(vect[24]))
  do.call(runemjmcmc, vect[1:23])
  vals<-values(hashStat)
  fparam<-mySearch$fparam
  cterm<-max(vals[1,],na.rm = T)
  ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
  post.populi<-sum(exp(values(hashStat)[1,][1:NM]-cterm),na.rm = T)
  clear(hashStat)
  rm(hashStat)
  rm(vals)
  return(list(post.populi = post.populi, p.post =  ppp$p.post, cterm = cterm, fparam = fparam))
}




simplifyposteriors<-function(X,posteriors,th=0.0001,thf=0.5)
{
  posteriors<-posteriors[-which(posteriors[,2]<th),]
  rhash<-hash()
  for(i in 1:length(posteriors[,1]))
  {
    expr<-posteriors[i,1]
    print(expr)
    res<-model.matrix(data=X,object = as.formula(paste0("Class1~",expr)))
    res[,1]<-res[,1]-res[,2]
    ress<-c(stri_flatten(res[,1],collapse = ""),stri_flatten(res[,2],collapse = ""),posteriors[i,2],expr)
    if(!(ress[1] %in% values(rhash)||(ress[2] %in% values(rhash))))
      rhash[[ress[1]]]<-ress
    else
    {
      if(ress[1] %in% keys(rhash))
      {
        rhash[[ress[1]]][3]<- (as.numeric(rhash[[ress[1]]][3]) + as.numeric(ress[3]))
        if(stri_length(rhash[[ress[1]]][4])>stri_length(expr))
          rhash[[ress[1]]][4]<-expr
      }
      else
      {
        rhash[[ress[2]]][3]<- (as.numeric(rhash[[ress[2]]][3]) + as.numeric(ress[3]))
        if(stri_length(rhash[[ress[2]]][4])>stri_length(expr))
          rhash[[ress[2]]][4]<-expr
      }
    }
    
  }
  res<-as.data.frame(t(values(rhash)[c(3,4),]))
  res$V1<-as.numeric(as.character(res$V1))
  res<-res[which(res$V1>thf),]
  res<-res[order(res$V1, decreasing = T),]
  clear(rhash)
  rm(rhash)
  res[which(res[,1]>1),1]<-1
  colnames(res)<-c("posterior","tree")
  return(res)
}




j=1
for(j in 1:MM)
{
  
  set.seed(j)
 
 X<-read.csv("qtlalzg.csv")[,-1]


  formula1 = as.formula(paste(colnames(X)[1],"~ 1 +",paste0(colnames(X)[-1],collapse = "+")))
  data.example = as.data.frame(X)
  
  print(formula1)

  #estimate.logic.lm(formula = formula1,data = data.example,n = 2000,m = 217)$mlik

  vect<-list(formula = formula1,data = X,outgraphs=F,estimator = estimate.logic.glm,locstop=F,presearch=T ,estimator.args =  list(data = data.example,n = 380, m = 34),recalc_margin = 249, save.beta = F,interact = T,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=1,mutation_rate = 250,last.mutation = 10000, max.tree.size = 4, Nvars.max =50,p.allow.replace=0.9,p.allow.tree=0.2,p.nor=0.1,p.and = 0.7),n.models = 35000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,print.freq = 1000,advanced.param = list(
    max.N.glob=as.integer(10),
    min.N.glob=as.integer(5),
    max.N=as.integer(3),
    min.N=as.integer(1),
    printable = F))
  
  params <- list(vect)[rep(1,M)]
  
  for(i in 1:M)
  {
    params[[i]]$cpu<-i
    params[[i]]$simul<-"scenario_R_"
    params[[i]]$simid<-j
  }
  gc()

  print(paste0("begin simulation ",j))
  results<-parall.gmj(X = params,FUN = runpar,mc.preschedule = T, mc.cores = M)
  #print(results)
  wait()


resa<-array(data = 0,dim = c(compmax,M*3))
post.popul <- array(0,M)
max.popul <- array(0,M)
nulls<-NULL
not.null<-1
for(k in 1:M)
{
  if(length(results[[k]])<4)
  {
    nulls<-c(nulls,k)
    next
  }
  else
  {
    not.null <- k
  }
  
}


for(kk in 1:M)
{
  if(kk %in% nulls)
  {
    k<-not.null
  }
  else
  {
    k <- kk
  }
  max.popul[k]<-results[[k]]$cterm
  post.popul[k]<-results[[k]]$post.populi
  resa[,k*3-2]<-c(results[[k]]$fparam,"Post.Gen.Max")
  resa[,k*3-1]<-c(results[[k]]$p.post,results[[k]]$cterm)
  resa[,k*3]<-rep(post.popul[k],length(results[[k]]$p.post)+1)
  
}

  gc()
  rm(results)

  ml.max<-max(max.popul)
  post.popul<-post.popul*exp(-ml.max+max.popul)
  p.gen.post<-post.popul/sum(post.popul)
  hfinal<-hash()
  for(ii in 1:M)
  {
    resa[,ii*3]<-p.gen.post[ii]*as.numeric(resa[,ii*3-1])
    resa[length(resa[,ii*3]),ii*3]<-p.gen.post[ii]
    if(p.gen.post[ii]>0)
    {
      for(jj in 1:(length(resa[,ii*3])-1))
      {
        if(resa[jj,ii*3]>0)
        {
          #print(paste0(ii,"  and ",jj))
          if(as.integer(has.key(hash = hfinal,key =resa[jj,ii*3-2]))==0)
            hfinal[[resa[jj,ii*3-2]]]<-as.numeric(resa[jj,ii*3])
          else
            hfinal[[resa[jj,ii*3-2]]]<-hfinal[[resa[jj,ii*3-2]]]+as.numeric(resa[jj,ii*3])
        }
        
      }
    }
  }
  
  posteriors<-values(hfinal)
  clear(hfinal)
  rm(hfinal)
  rm(resa)
  rm(post.popul)
  rm(max.popul)
  posteriors<-as.data.frame(posteriors)
  posteriors<-data.frame(X=row.names(posteriors),x=posteriors$posteriors)
  posteriors$X<-as.character(posteriors$X)
  tryCatch({
    res1<-simplifyposteriors(X = X,posteriors = posteriors, th = th,thf = thf)
    write.csv(x =res1,row.names = F,file = paste0("postRealAlzg_",j,".csv"))
  },error = function(err){
    print("error")
    write.csv(x =posteriors,row.names = F,file = paste0("posteriorsRealAlzg_",j,".csv"))
  },finally = {
    
    print(paste0("end simulation ",j))
    
  }) 
  rm(X)
  rm(data.example)
  rm(vect)
  rm(params)
  gc()
  print(paste0("end simulation ",j))
  
}
