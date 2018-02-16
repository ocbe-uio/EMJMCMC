source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package2.r")

library(inline)
includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')

parall.gmj <<- mclapply

MM = 1
M = 2
NM= 1000
compmax = 26
th<-(10)^(-5)
thf<-0.05


estimate.logic.lm <- function(formula, data, n, m, r = 0.001)
{
  r=1/m
  out <- lm(formula = formula,data = data)
  p <- out$rank
  fmla.proc<-as.character(formula)[2:3]
  fobserved <- fmla.proc[1]
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam <-stri_split_fixed(str = fmla.proc[2],pattern = "+",omit_empty = F)[[1]]
  sj<-(stri_count_fixed(str = fparam, pattern = "&"))
  sj<-sj+(stri_count_fixed(str = fparam, pattern = "|"))
  sl<-length(which(sj>0))
  #mlik = (-(out$deviance + log(n)*(out$rank+sum(sj)-sl)))/2
  mlik = (-BIC(out)+2*log(r)*(sum(sj)+out$rank))/2
  if(mlik==-Inf)
    mlik = -10000
  return(list(mlik = mlik,waic = AIC(out)-n , dic =  BIC(out)-n,summary.fixed =list(mean = coef(out))))
}


paral<-function(X,FUN)
{
  return(mclapply(X = X,FUN = FUN,mc.preschedule = F, mc.cores = 32))
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



for(j in 1:MM)
{

  set.seed(j)

  X<-read.csv("XS.csv",header = F)[,-1]
  XL<-length(names(X))
  names(X)[XL]<-"Y"
  i1<-which(X[,-XL]>=0)
  arr<-array(data = 0,dim = dim(X[,-XL])[1]*dim(X[,-XL])[2])
  arr[i1]<-1
  X[,-XL]<-array(data = arr,dim = c(dim(X[,-XL])[1],dim(X[,-XL])[2]))
  
  X<-X[-which(is.na(X$Blue)),]
  #X$Blue<-as.numeric(as.character(X$Blue))
  idss<-which(abs(cor(x = X,y=X$Y))>0.0)

  formula1 = as.formula(paste(colnames(X)[XL],"~ 1 +",paste0(colnames(X)[idss][-length(idss)],collapse = "+")))
  data.example = as.data.frame(X)

  print(formula1)

  #estimate.logic.lm(formula = formula1,data = data.example,n = 2000,m = 217)$mlik

  vect<-list(formula = formula1,data = X,outgraphs=F,estimator = estimate.logic.lm,locstop=T,presearch=T ,estimator.args =  list(data = data.example,n = 234, m = 221),recalc_margin = 249, save.beta = F,interact = T,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=1,mutation_rate = 250,last.mutation = 10000, max.tree.size = 4, Nvars.max =45,p.allow.replace=0.9,p.allow.tree=0.02,p.nor=0.1,p.and = 0.7),n.models = 35000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,print.freq = 1000,advanced.param = list(
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
  for(k in 1:M)
  {
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
    res1<-simplifyposteriors(X = X,posteriors = posteriors, th = th,thf = thf,y="T1G11~")
    write.csv(x =res1,row.names = F,file = paste0("postReal_",j,".csv"))
  },error = function(err){
    print("error")
    write.csv(x =posteriors,row.names = F,file = paste0("posteriorsReal_",j,".csv"))
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
