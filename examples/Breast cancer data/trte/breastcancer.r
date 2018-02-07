

source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")


estimate.bas.glm.cpen <- function(formula, data, family, prior, logn,r = 0.1,yid=1,relat =c("logi","sigmoid","tanh","atan","erf"))
{
  
  #only poisson and binomial families are currently adopted
  capture.output({out <- glm(family = family,formula = formula,data = data)})
  fmla.proc<-as.character(formula)[2:3]
  fobserved <- fmla.proc[1]
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  #fparam <-stri_split_fixed(str = fmla.proc[2],pattern = "+I",omit_empty = F)[[1]]
  #sj<-(stri_count_fixed(str = fparam, pattern = "("))
  sj<-2*(stri_count_fixed(str = fmla.proc[2], pattern = "*"))
  sj<-sj+1*(stri_count_fixed(str = fmla.proc[2], pattern = "+"))
  for(rel in relat)
    sj<-sj+2*(stri_count_fixed(str = fmla.proc[2], pattern = rel))
  #sj<-sj+1
  
  mlik = ((-out$deviance +2*log(r)*sum(sj)))/2
  
  #print(sj)
  #print(sum(sj))
  return(list(mlik = mlik,waic = -(out$deviance + 2*out$rank) , dic =  -(out$deviance + logn*out$rank),summary.fixed =list(mean = coefficients(out))))
  
}

set.seed(10)


parall.gmj <<- mclapply


test <- read.csv("test.csv",header = T,sep=",")[,-1]
train <- read.csv("train.csv",header = T,sep=",")[,-1]

data.example <- as.data.frame(train)

gc()


total = array(0,dim = c(100,10,3))


g<-function(x)
{
  return((x = 1/(1+exp(-x))))
}


runpar<-function(vect)
{
  
  set.seed(as.integer(vect$seed))
  do.call(runemjmcmc, vect[1:vect$runlen])
  
  ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
  ppp$p.post
  
 
  Nvars<-mySearch$Nvars
  linx <-mySearch$Nvars+4
  lHash<-length(hashStat)
  mliks <- values(hashStat)[which((1:(lHash * linx)) %% linx == 1)]
  betas <- values(hashStat)[which((1:(lHash * linx)) %% linx == 4)]
  cterm<-max(values(hashStat)[1,],na.rm = T)
  post.populi<-sum(exp(values(hashStat)[1,][1:1000]-cterm),na.rm = T)
  for(i in 1:(Nvars-1))
  {
    betas<-cbind(betas,values(hashStat)[which((1:(lHash * linx)) %% linx == (4+i))])
  }
  betas<-cbind(betas,values(hashStat)[which((1:(lHash * linx)) %% linx == (0))])
  
  t<-system.time({
    
    res<-mySearch$forecast.matrix.na(link.g = g, covariates = (vect$test),betas = betas,mliks.in = mliks)$forecast
    
  })
  
  rm(betas)
  rm(mliks)
  clear(hashStat)
  rm(hashStat)
  return(list(p.post =  ppp$p.post, fparam = mySearch$fparam, res = res, post.populi = post.populi, cterm = cterm))
}



gc()

gone<-function(x)as.integer(x>1)
gthird<-function(x)as.integer(abs(x)^(1/3))
gfifth<-function(x)as.integer(abs(x)^(1/5))
grelu<-function(x)as.integer(x*(x>0))

M<-32
for(j in 19:100)
{
  print(j)
  set.seed(j)
  formula1 = as.formula(paste(colnames(data.example)[31],"~ 1 +",paste0(colnames(data.example)[-31],collapse = "+")))
  
  vect<-list(formula = formula1,data = data.example,gen.prob = c(1,1,1,1,1),estimator =estimate.bas.glm.cpen,estimator.args =  list(data = data.example,prior = aic.prior(),family = binomial(),yid=31, logn = log(143),r=exp(-0.5)),recalc_margin = 95, save.beta = T,interact = T,relations = c("sigmoid","grelu","gone","gthird","gfifth"),relations.prob =c(0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=3,mutation_rate = 100,last.mutation=1000, max.tree.size = 6, Nvars.max = 20,p.allow.replace=0.5,p.allow.tree=0.4,p.nor=0.3,p.and = 0.9),n.models = 7000,unique =T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
    max.N.glob=as.integer(10),
    min.N.glob=as.integer(5),
    max.N=as.integer(3),
    min.N=as.integer(1),
    printable = F))
  
  len = length(vect)
  
  params <- list(vect)[rep(1,M)]
  
  for(jj in 1:M)
  {
    params[[jj]]$runlen = len
    params[[jj]]$cpu<-jj
    params[[jj]]$seed = jj*j
    params[[jj]]$test<-test
  }
  gc()
  
  length(params[[1]])
  
  results<-parall.gmj(X = params,FUN = runpar,mc.preschedule = F, mc.cores = M,mc.cleanup = T)

  post.popul <- array(0,M)
  max.popul <- array(0,M)
  nulls = NULL
  for(k in 1:M)
  {
    if(length(results[[k]])==1||length(results[[k]]$cterm)==0)
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
  }
  
  ml.max<-max(max.popul)
  post.popul<-post.popul*exp(-ml.max+max.popul)
  p.gen.post<-post.popul/sum(post.popul)

  res1 = results[[1]]$res*p.gen.post[1]
  for(i in 2:M)
  {
    
    res1<-res1+results[[i]]$res*p.gen.post[i]
  
  }  
  
  for(jjjj in 1:10)
  {
    res = as.integer(res1>=0.1*jjjj)
    prec<-(1-sum(abs(res-test$X),na.rm = T)/length(res))
    
    #FNR
    ps<-which(test$X==1)
    fnr<-sum(abs(res[ps]-test$X[ps]))/(sum(abs(res[ps]-test$X[ps]))+length(ps))
    
    #FPR
    ns<-which(test$X==0)
    fpr<-sum(abs(res[ns]-test$X[ns]))/(sum(abs(res[ns]-test$X[ns]))+length(ns))
      
    total[j,jjjj,1]=prec
    total[j,jjjj,2]=fnr
    total[j,jjjj,3]=fpr
    
    
   
  }
  write.csv(file ="results.csv",x=total)
  print(max(total[j,,1]))
  rm(results)
  gc()
}

plot(x = 1:100,y = total[,1,1],type = "l", col = 0,ylim =c(0.935,0.9813))
for(i in 2:9)
  lines(total[,i,1],type = "l", col = i)


medis = array(0,dim = c(10,3))
mins = array(0,dim = c(10,3)) 
maxs = array(0,dim = c(10,3)) 
for(i in 1:10)
  for(j in 1:3)
  {
    medis[i,j] = median(total[,i,j])
    mins[i,j] = min(total[,i,j])
    maxs[i,j] = max(total[,i,j])
  }


final = cbind(0.1*1:10,cbind((medis),(mins),(maxs))[,c(4,1,7,5,2,8,6,3,9)])


