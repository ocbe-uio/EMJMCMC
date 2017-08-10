# ssh -X -Y -l aliaksah abel.uio.no
# scp -r  /usit/abel/u1/aliaksah/simulations/scenario1  aliaksah@pittheus.uio.no://mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/simulations
# cat slurm-16078690.out
# squeue -u aliaksah
#

source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")


library(inline)
includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')



estimate.dlm <- function(formula = NA, data, n, m, r = 1,sigmas = c("cosi","sigmoid","tanh","atan","sini","troot"))
{
  if(is.na(formula))
  {
    print("FORMULA MISSING")
    return(NULL)
  }
  out <- lm(formula = formula,data = data)
  p <- out$rank
  fmla.proc<-as.character(formula)[2:3]
  fobserved <- fmla.proc[1]
  sj<-(stri_count_fixed(str = fmla.proc[2], pattern = "*"))
  sj<-sj+(stri_count_fixed(str = fmla.proc[2], pattern = "+"))
  sj<-sj+sum(stri_count_fixed(str = fmla.proc[2], pattern = sigmas))
  sj<-sj-p+1
  #Jprior <- prod(factorial(sj)/((m^sj)*2^(2*sj-2)))
  #tn<-sum(stri_count_fixed(str = fmla.proc[2], pattern = "I("))
  mlik = (sj<=5)*((-BIC(out) - (m-1)*p*log(n) - m*sj*log(n))/2) + (sj>5)*(-10000)
  if(is.na(mlik))
    mlik = -10000
  if(mlik==-Inf)
    mlik = -10000
  #print(sj)
  return(list(mlik = mlik,waic = AIC(out)-n , dic =  BIC(out)-n,summary.fixed =list(mean = coef(out))))
}



parall.gmj <<- mclapply



simplifyposteriors<-function(X,posteriors,th=0.0001,thf=0.3)
{
  posteriors<-posteriors[-which(posteriors[,2]<th),]
  rhash<-hash()
  for(i in 1:length(posteriors[,1]))
  {
    expr<-posteriors[i,1]
    print(expr)
    res<-model.matrix(data=X,object = as.formula(paste0("RadiusJpt~",expr)))
    res[,1]<-res[,1]-res[,2]
    ress<-c(stri_flatten(res[,1],collapse = ""),stri_flatten(res[,2],collapse = ""),posteriors[i,2],expr)
    if(!((ress[2] %in% values(rhash))))
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

cosi<-function(x)cos(x/180*pi)
sini<-function(x)sin(x/180*pi)
expi<-function(x)
{
  r<-exp(x)
  if(r==Inf)
    return(10000000)
  else
    return(r)
}

InvX<-function(x)
{
  if(x==0)
    return(10000000)
  else
    return(1/x)

}
troot<-function(x)abs(x)^(1/3)

MM = 100
M = 4
NM= 1000
compmax = 16
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
  gc()
  return(list(post.populi = post.populi, p.post =  ppp$p.post, cterm = cterm, fparam = fparam))
}


for(j in 1:100)
{
  tryCatch({

    set.seed(j)

    X<-read.csv("/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/exaplanets/exa1.csv")

    formula1 = as.formula(paste(colnames(X)[5],"~ 1 +",paste0(colnames(X)[-5],collapse = "+")))
    data.example = as.data.frame(X)

    print(formula1)

    #wait()

    vect<-list(formula = formula1,data = data.example,estimator =estimate.dlm,estimator.args =  list(data = data.example,n = 223, m = 2),recalc_margin = 249, save.beta = F,interact = T,outgraphs=F,relations=c("cosi","sigmoid","tanh","atan","sini","troot"),relations.prob =c(0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=3,mutation_rate = 250,last.mutation=10000, max.tree.size = 5, Nvars.max =15,p.allow.replace=0.9,p.allow.tree=0.01,p.nor=0.9,p.and = 0.9),n.models = 100000,unique =F,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 100,advanced.param = list(
      max.N.glob=as.integer(10),
      min.N.glob=as.integer(5),
      max.N=as.integer(3),
      min.N=as.integer(1),
      printable = F))

    aaa=do.call(runemjmcmc,vect[1:21])
    aaa$p.post

    estimate.dlm (data = data.example,formula = SemiMajorAxisAU ~ 1 + I(troot((HostStarMassSlrMass)*PeriodDays*PeriodDays)) ,n = 223, m = 2)



    estimate.dlm(data = data.example,formula =  as.formula(paste(colnames(X)[5],"~ 1 +",paste0(mySearch$fparam[which(aaa$p.post>0.99)],collapse = "+"))),n = 223, m = 2)

    estimate.dlm(data = data.example,formula =  as.formula(paste(colnames(X)[5],"~ 1 +",paste0(mySearch$fparam[10],collapse = "+"))),n = 223, m = 2)


    params <- list(vect)[rep(1,M)]

    for(i in 1:M)
    {
      params[[i]]$cpu<-i*j
      params[[i]]$simul<-"scenario_JM_"
      params[[i]]$simid<-j
    }
    gc()
    print(paste0("begin simulation ",j))
    results<-parall.gmj(X = params,FUN = runpar,mc.preschedule = F, mc.cores = M)
    #print(results)

    wait()

    resa<-array(data = 0,dim = c(compmax,M*3))
    post.popul <- array(0,M)
    max.popul <- array(0,M)
    nulls<-NULL
    not.null<-1
    for(k in 1:M)
    {
      if(length(results[[k]])==0)
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

    print(posteriors)
    clear(hfinal)
    rm(hfinal)
    rm(resa)
    rm(post.popul)
    rm(max.popul)
    posteriors<-as.data.frame(posteriors)
    posteriors<-data.frame(X=row.names(posteriors),x=posteriors$posteriors)
    posteriors$X<-as.character(posteriors$X)
    tryCatch({
      res1<-simplifyposteriors(X = X,posteriors = posteriors, th,thf)
      write.csv(x =res1,row.names = F,file = paste0("postJA32_",j,".csv"))
    },error = function(err){
      print("error")
      write.csv(x =posteriors,row.names = F,file = paste0("posteriorsJA32_",j,".csv"))
    },finally = {

      print(paste0("end simulation ",j))

    })
    rm(X)
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
    rm(X)
    rm(data.example)
    rm(vect)
    rm(params)
    gc()
  })

}
