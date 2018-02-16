# ssh -X -Y -l aliaksah abel.uio.no
# scp -r  /usit/abel/u1/aliaksah/simulations/scenario1  aliaksah@pittheus.uio.no://mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/simulations
# cat slurm-16055503.out
# squeue -u aliaksah

source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package2.r")

parall.gmj <<- mclapply

library(inline)
includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')



estimate.logic.lm.tCCH <- function(formula = NULL, data, n=1000, m=50, r = 1, p.a = 1, p.b = 2, p.r = 1.5, p.s = 0, p.v=-1, p.k = 1,k.max=21)
{
  if(is.na(formula)||is.null(formula))
    return(list(mlik =  -10000,waic =10000 , dic =  10000,summary.fixed =list(mean = 1)))
  fmla.proc<-as.character(formula)[2:3]
  fobserved <- fmla.proc[1]
  if(fmla.proc[2]=="-1")
    return(list(mlik =  -10000,waic =10000 , dic =  10000,summary.fixed =list(mean = 1)))
  out <- lm(formula = formula,data = data)
  p <- out$rank
  if(p>k.max)
  {
    return(list(mlik = -10000,waic = 10000 , dic =  10000,summary.fixed =list(mean = 0)))
  }
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam <-stri_split_fixed(str = fmla.proc[2],pattern = "+",omit_empty = F)[[1]]
  sj<-(stri_count_fixed(str = fparam, pattern = "&"))
  sj<-sj+(stri_count_fixed(str = fparam, pattern = "|"))
  sj<-sj+1
  Jprior <- prod(factorial(sj)/((m^sj)*2^(2*sj-2)))
  #tn<-sum(stri_count_fixed(str = fmla.proc[2], pattern = "I("))
  p.v = (n+1)/(p+1)
  R.2 = summary(out)$r.squared

  mlik = (-0.5*p*log(p.v) -0.5*(n-1)*log(1-(1-1/p.v)*R.2) + log(beta((p.a+p)/2,p.b/2)) - log(beta(p.a/2,p.b/2)) + log(phi1(p.b/2,(n-1)/2,(p.a+p.b+p)/2,p.s/2/p.v,R.2/(p.v-(p.v-1)*R.2))) - hypergeometric1F1(p.b/2,(p.a+p.b)/2,p.s/2/p.v,log = T)+log(Jprior) + p*log(r)+n)
  if(mlik==-Inf||is.na(mlik)||is.nan(mlik))
    mlik = -10000
  #if(mlik== -Inf)
  #  mlik = 10000
  #print(mlik)
  return(list(mlik = mlik,waic = AIC(out)-n , dic =  BIC(out)-n,summary.fixed =list(mean = coef(out))))
}


estimate.logic.lm <- function(formula= NULL, data, n, m, r = 1,k.max=21)
{
  if(is.na(formula)||is.null(formula))
    return(list(mlik =  -10000,waic =10000 , dic =  10000,summary.fixed =list(mean = 1)))
  out <- lm(formula = formula,data = data)
  p <- out$rank
  if(p>k.max)
  {
    return(list(mlik = -10000,waic = 10000 , dic =  10000,summary.fixed =list(mean = 0)))
  }
  fmla.proc<-as.character(formula)[2:3]
  fobserved <- fmla.proc[1]
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam <-stri_split_fixed(str = fmla.proc[2],pattern = "+",omit_empty = F)[[1]]
  sj<-(stri_count_fixed(str = fparam, pattern = "&"))
  sj<-sj+(stri_count_fixed(str = fparam, pattern = "|"))
  sj<-sj+1
  Jprior <- prod(factorial(sj)/((m^sj)*2^(2*sj-2)))
  #tn<-sum(stri_count_fixed(str = fmla.proc[2], pattern = "I("))
  mlik = (-BIC(out)+2*log(Jprior) + 2*p*log(r)+n)/2
  if(mlik==-Inf)
    mlik = -10000
  return(list(mlik = mlik,waic = AIC(out)-n , dic =  BIC(out)-n,summary.fixed =list(mean = coef(out))))
}

simplifyposteriors<-function(X,posteriors,th=0.0001,thf=0.5)
{
  posteriors<-posteriors[-which(posteriors[,2]<th),]
  rhash<-hash()
  for(i in 1:length(posteriors[,1]))
  {
    expr<-posteriors[i,1]
    print(expr)
    res<-model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
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


MM = 130
M = 32
NM= 1000
compmax = 21
th<-(10)^(-5)
thf<-0.05

paral<-function(X,FUN)
{
  return(mclapply(X = X,FUN = FUN,mc.preschedule = F, mc.cores = 32))
}

runpar<-function(vect)
{

  set.seed(as.integer(vect[22]))
  do.call(runemjmcmc, vect[1:21])
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


for(j in 960:916)
{
  tryCatch({

    j.id = j%%60

    jj=j
    compmax = 21
    if(j.id%in%1:10)
    {
      type1=T
      type2=F
      type3=F
      typeJ=T
      jj=j.id
      id.sym = "T1J"
    }else if(j.id%in%11:20)
    {
      type1=T
      type2=F
      type3=F
      typeJ=F
      jj=j.id-10
      id.sym = "T1G"
    }else if(j.id%in%21:30)
    {
      type1=F
      type2=T
      type3=F
      typeJ=T
      jj=j.id-20
      id.sym = "T2J"
    }else if(j.id%in%31:40)
    {
      type1=F
      type2=T
      type3=F
      typeJ=F
      jj=j.id-30
      id.sym = "T2G"
    }else if(j.id%in%41:50)
    {
      next
      type1=F
      type2=F
      type3=T
      typeJ=T
      jj=j.id-40
      compmax = jj*15+1
      id.sym = "T3J"
    }else if(j.id%in%51:60)
    {
      next
      type1=F
      type2=F
      type3=T
      typeJ=F
      jj=j.id-50
      compmax = jj*15+1
      id.sym = "T3G"
    }


    set.seed(jj)
    X2<- as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))
    Y2<-rnorm(n = 1000,mean = 1+7*(X2$V4*X2$V17*X2$V30*X2$V10) + 9*(X2$V7*X2$V20*X2$V12)+ 3.5*(X2$V9*X2$V2)+1.5*X2$V37,sd = 1)
    X2$Y2<-Y2

    formula1 = as.formula(paste(colnames(X2)[51],"~ 1 +",paste0(colnames(X2)[-c(51)],collapse = "+")))
    data.example = as.data.frame(X2)

    vect<-list(formula = formula1,data = X2,estimator = estimate.logic.lm.tCCH,estimator.args =  list(data = data.example,n = 1000, m = 50),recalc_margin = 250, save.beta = F,interact = T,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=1,mutation_rate = 300,last.mutation = 5000, max.tree.size = 4, Nvars.max =(compmax-1),p.allow.replace=0.9,p.allow.tree=0.2,p.nor=0,p.and = 0.9),n.models = 10000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,outgraphs=F,print.freq = 1000,advanced.param = list(
      max.N.glob=as.integer(10),
      min.N.glob=as.integer(5),
      max.N=as.integer(3),
      min.N=as.integer(1),
      printable = F))

    params <- list(vect)[rep(1,32)]

    for(i in 1:M)
    {
      if(type1)
      {
        Y2<-rnorm(n = 1000,mean = 1+jj*(X2$V4*X2$V17*X2$V30*X2$V10) + 9*(X2$V7*X2$V20*X2$V12)+ 3.5*(X2$V9*X2$V2)+1.5*X2$V37,sd = 1)
        X2$Y2<-Y2
        data.example = as.data.frame(X2)
        params[[i]]$data = data.example
        params[[i]]$estimator.args =  list(data = data.example,n = 1000, m = 50)
        params[[i]]$formula = formula1
        if(typeJ)
        {
          params[[i]]$estimator = estimate.logic.lm
        }
      }else if(type2)
      {
        NN = jj*100
        X2<- as.data.frame(array(data = rbinom(n = 50*NN,size = 1,prob = runif(n = 50*NN,0,1)),dim = c(NN,50)))
        Y2<-rnorm(n = NN,mean = 1+7*(X2$V4*X2$V17*X2$V30*X2$V10) + 9*(X2$V7*X2$V20*X2$V12)+ 3.5*(X2$V9*X2$V2)+1.5*X2$V37,sd = 1)
        X2$Y2<-Y2
        data.example = as.data.frame(X2)
        params[[i]]$data = data.example
        params[[i]]$estimator.args =  list(data = data.example,n = NN, m = 50)
        params[[i]]$formula = formula1
        if(typeJ)
        {
          params[[i]]$estimator = estimate.logic.lm
        }
      }else if(type3)
      {
        params[[i]]$interact.param=list(allow_offsprings=1,mutation_rate = 300,last.mutation = 5000, max.tree.size = 4, Nvars.max =jj*15,p.allow.replace=0.9,p.allow.tree=0.2,p.nor=0,p.and = 0.9)
        if(typeJ)
        {
          params[[i]]$estimator = estimate.logic.lm
        }

      }

      params[[i]]$cpu<-i
      params[[i]]$simul<-"scenario_2_"
      params[[i]]$simid<-jj
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
      if(length(results[[k]])<=1||length(results[[k]]$cterm)==0)
      {
        nulls<-c(nulls,k)
        next
      }
      else
      {
        not.null <- k
      }

    }


    for(k in 1:M)
    {
      if(k %in% nulls)
      {
        results[[k]]<-results[[not.null]]
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
        for(jjj in 1:(length(resa[,ii*3])-1))
        {
          if(resa[jjj,ii*3]>0)
          {
            #print(paste0(ii,"  and ",jj))
            if(as.integer(has.key(hash = hfinal,key =resa[jjj,ii*3-2]))==0)
              hfinal[[resa[jjj,ii*3-2]]]<-as.numeric(resa[jjj,ii*3])
            else
              hfinal[[resa[jjj,ii*3-2]]]<-hfinal[[resa[jjj,ii*3-2]]]+as.numeric(resa[jjj,ii*3])
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
      res1<-simplifyposteriors(X = X2,posteriors = posteriors, th,thf)
      write.csv(x =res1,row.names = F,file = paste0("post2eta_",id.sym,"_",j,".csv"))
    },error = function(err){
      print("error")
      write.csv(x =posteriors,row.names = F,file = paste0("posteriors2eta_",id.sym,"_",j,".csv"))
    },finally = {

      print(paste0("end simulation ",j))

    })
    rm(X2)
    rm(data.example)
    rm(vect)
    rm(params)
    gc()
    print(paste0("end simulation ",j))
  },error = function(err){
    print("error")
    print(err)
    print(paste0("repeat  simulation ",j))
  },finally = {

    print(paste0("end simulation ",j))

  })


}


