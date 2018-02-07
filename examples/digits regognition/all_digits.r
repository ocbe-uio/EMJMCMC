#compare on Asteroids
library(RCurl)
library(caret)

source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")



estimate.glm.cpen <- function(formula, data, family, logn,r = 0.1,relat =c("gone","gthird","sigmoid","tanh","atan","erf","gfifth","grelu"))
{
  
  capture.output({out <- glm(family = family,formula = formula,data = data)})
  fmla.proc<-as.character(formula)[2:3]
  fobserved <- fmla.proc[1]
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  sj<-2*(stri_count_fixed(str = fmla.proc[2], pattern = "*"))
  sj<-sj+1*(stri_count_fixed(str = fmla.proc[2], pattern = "+"))
  for(rel in relat)
    sj<-sj+2*(stri_count_fixed(str = fmla.proc[2], pattern = rel))
  mlik = ((-out$deviance +2*log(r)*sum(sj)))/2
  return(list(mlik = mlik,waic = -(out$deviance + 2*out$rank) , dic =  -(out$deviance + logn*out$rank),summary.fixed =list(mean = coefficients(out))))
  
}


parall.gmj <<- mclapply

digits <- read.table(file = "train.csv",sep = ",",header = T,fill=TRUE)
digits$Y1<-as.integer(digits$label==1)
digits$Y2<-as.integer(digits$label==2)
digits$Y3<-as.integer(digits$label==3)
digits$Y4<-as.integer(digits$label==4)
digits$Y5<-as.integer(digits$label==5)
digits$Y6<-as.integer(digits$label==6)
digits$Y7<-as.integer(digits$label==7)
digits$Y8<-as.integer(digits$label==8)
digits$Y9<-as.integer(digits$label==9)
digits$Y10<-as.integer(digits$label==10)


digits<-digits[,-which(colSums(digits)==0)]


g<-function(x)
{
  return((x = 1/(1+exp(-x))))
}


index <- createDataPartition(digits$label, p = 0.01, list = FALSE)

test <- digits[-index, ]
train <- digits[index, ]

data.example <- as.data.frame(train,stringsAsFactors = T)


runpar<-function(vect)
{
  
  set.seed(as.integer(vect[23]))
  do.call(runemjmcmc, vect[1:22])
  
  ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
  ppp$p.post
  
 
  Nvars<-mySearch$Nvars
  linx <-mySearch$Nvars+4
  lHash<-length(hashStat)
  mliks <- values(hashStat)[which((1:(lHash * linx)) %% linx == 1)]
  betas <- values(hashStat)[which((1:(lHash * linx)) %% linx == 4)]
  cterm<-max(vals[1,],na.rm = T)
  post.populi<-sum(exp(values(hashStat)[1,][1:1000]-cterm),na.rm = T)
  for(i in 1:(Nvars-1))
  {
    betas<-cbind(betas,values(hashStat)[which((1:(lHash * linx)) %% linx == (4+i))])
  }
  betas<-cbind(betas,values(hashStat)[which((1:(lHash * linx)) %% linx == (0))])
  
  t<-system.time({
    
    res<-mySearch$forecast.matrix.na(link.g = g, cterm = cterm, post.populi = post.populi, covariates = (vect$test),betas = betas,mliks.in = mliks)$forecast
    
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

total = array(0,dim = c(10,3))

M<-32
for(dig in 4:1)
{
  idss<-which(abs(cor(x = train[1:785],y=train[785+dig]))>0.1)
  formula1 = as.formula(paste(colnames(train)[785+dig],"~ 1 +",paste0(colnames(train)[idss][-1],collapse = "+")))
  
  vect<-list(formula = formula1,data = data.example,estimator =estimate.glm.cpen,estimator.args =  list(data = data.example,family = binomial(), logn = log(dim(train)[1]),r=exp(-0.5)),recalc_margin = 95,locstop=T,presearch=F,save.beta = T,interact = T,relations = c("gone","gthird","sigmoid","tanh","atan","erf","gfifth","grelu"),relations.prob =c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=2,mutation_rate = 100,last.mutation=2500, max.tree.size = 20, Nvars.max =70,p.allow.replace=0.3,p.allow.tree=0.5,p.nor=0.3,p.and = 0.7),n.models = 30000,unique =F,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
      max.N.glob=as.integer(10),
      min.N.glob=as.integer(5),
      max.N=as.integer(3),
      min.N=as.integer(1),
      printable = F))
  
  length(vect)
  
  params <- list(vect)[rep(1,M)]
  
  for(jj in 1:M)
  {
    params[[jj]]$cpu<-jj
    params[[jj]]$test<-test
  }
  gc()
  
  length(params[[1]])
  
  results<-parall.gmj(X = params,FUN = runpar,mc.preschedule = F, mc.cores = M,mc.cleanup = T)

  
  cbind(results[[1]]$fparam,results[[1]]$p.post)
  
  post.popul <- array(0,M)
  max.popul <- array(0,M)
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

  for(i in 1:M)
  {
    
    res<-res+results[[i]]$res*p.gen.post[i]
  
  }  
  
  res = as.integer(res>0.1)
  prec<-(1-sum(abs(res-test[,785+dig]),na.rm = T)/length(res))
  
  #FNR
  ps<-which(test[,785+dig]==1)
  fnr<-sum(abs(res[ps]-test[,785+dig][ps]))/(sum(abs(res[ps]-test[,785+dig][ps]))+length(ps))
  
  #FPR
  ns<-which(test[,785+dig]==0)
  fpr<-sum(abs(res[ns]-test[,785+i][ns]))/(sum(abs(res[ns]-test[,785+i][ns]))+length(ns))
    
  total[dig,1]=prec
  total[dig,2]=fnr
  total[dig,3]=fpr
  
  rm(results)
  gc()
  print(prec)
  
  

}

