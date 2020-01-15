#inference
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")

setwd("/nr/samba/user/ahu/EMJMCMC2016/supplementaries/BGNLM/newspop/")


data.example = read.csv("OnlineNewsPopularity.csv",header = T)[,-c(1,2)]



set.seed(040590)
teid =  sample.int(size =30000,n = 39644,replace = F)

test = data.example[teid,]
data.example = data.example[-teid,]


#specify the initial formula
formula1 = as.formula(paste(colnames(test)[59],"~ 1 +",paste0(colnames(test)[-59],collapse = "+")))

#a set of nonlinearities that will be used in the DBRM model
sini=function(x)sin(x/180*pi)
expi=function(x)exp(-abs(x))
logi =function(x)log(abs(x)+1)
troot=function(x)abs(x)^(1/3)
to25=function(x)abs(x)^(2.5)
to35=function(x)abs(x)^(3.5)
#relu=function(x)max(0,x)

#specify the estimator function returning p(Y|m)p(m), model selection criteria and the vector of the modes for the beta coefficients
estimate.gamma.cpen = function(formula, data,r = 1.0/9644.0,logn=log(9644.0),relat=c("to25","expi","logi","to35","troot","sigmoid"))
{
  fparam=NULL
  fmla.proc=as.character(formula)[2:3]
  fobserved = fmla.proc[1]
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam =stri_split_fixed(str = fmla.proc[2],pattern = "+I",omit_empty = F)[[1]]
  sj=(stri_count_fixed(str = fparam, pattern = "*"))
  sj=sj+(stri_count_fixed(str = fparam, pattern = "+"))
  for(rel in relat)
    sj=sj+(stri_count_fixed(str = fparam, pattern = rel))
  sj=sj+1
  tryCatch(capture.output({
    out = glm(formula = formula,data = data, family = gaussian)
    mlik = (-(out$deviance -2*log(r)*sum(sj)))/2
    waic = -(out$deviance + 2*out$rank)
    dic =  -(out$deviance + logn*out$rank)
    summary.fixed =list(mean = coefficients(out))
    
  }, error = function(err) {
    print(err)
    mlik = -10000
    waic = -10000
    dic =  -10000
    summary.fixed =list(mean = array(0,dim=length(fparam)))
  }))
  return(list(mlik = mlik,waic = waic , dic = dic,summary.fixed =summary.fixed))
  
}
#define the number or cpus
M = 32
#define the size of the simulated samples
NM= 1000
#define \k_{max} + 1 from the paper
compmax = 16
#define treshold for preinclusion of the tree into the analysis
th=(10)^(-5)
#define a final treshold on the posterior marginal probability for reporting a tree
thf=0.05
#specify tuning parameters of the algorithm for exploring DBRM of interest
#notice that allow_offsprings=3 corresponds to the GMJMCMC runs and
#
g = function(x) x
results=array(0,dim = c(2,100,5))

for(j in 1:70)
{
  
  #specify the initial formula
  set.seed(j)

  res1 = pinferunemjmcmc(n.cores = M, report.level =  0.2, num.mod.best = NM ,simplify = T,predict = T,test.data = as.data.frame(test),link.function = g,runemjmcmc.params = list(formula = formula1,data = data.example,estimator = estimate.gamma.cpen,estimator.args =  list(data = data.example),recalc_margin = 249, save.beta = T,interact = T,outgraphs=F,relations=c("to25","expi","logi","to35","troot","sigmoid"),relations.prob =c(0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=3,mutation_rate = 250,last.mutation=60000, max.tree.size = 5, Nvars.max =75,p.allow.replace=0.9,p.allow.tree=0.01,p.nor=0.9,p.and = 0.9),n.models = 75000,unique =T,max.cpu = 4,max.cpu.glob = 4,create.table = F,presearch = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
    max.N.glob=as.integer(10),
    min.N.glob=as.integer(5),
    max.N=as.integer(3),
    min.N=as.integer(1),
    printable = F)))
  
  print(res1$feat.stat)
  

  for(j in 1:M)
  {
    if(length(res1$threads.stats[[j]])==0)
      next
    results[2,j,1]=  sqrt(mean((res1$threads.stats[[j]]$preds - test$shares)^2))
    results[2,j,2]=   (mean(abs(res1$threads.stats[[j]]$preds - test$shares)))
    results[2,j,3] =   cor(res1$threads.stats[[j]]$preds,test$shares)
    write.csv(x =res1$feat.stat,row.names = F,file = paste0("posteriorshellf_",j,".csv"))
    print(paste0("end simulation ",j))
    
    #print the run's metrics and clean the results
    write.csv(file =paste0("resultsrun_",j,".csv"),x= results[,j,])
    #rm(results)
  }
  
  results[2,j,1]=  sqrt(mean((res1$predictions - test$shares)^2))
  results[2,j,2]=   (mean(abs(res1$predictions - test$shares)))
  results[2,j,3] =   cor(res1$predictions,test$shares)
 
  write.csv(x =res1$feat.stat,row.names = F,file = paste0("posteriorshellf_",j,".csv"))
  print(paste0("end simulation ",j))
  
  #print the run's metrics and clean the results
  write.csv(file =paste0("resultsrun_",j,".csv"),x= results[,j,])
  #rm(results)
  
  
  print(sqrt(mean((res1$predictions - test$shares)^2)))
  
  gc()
}


for(j in 1:100)
{
  
  tmp = read.csv(paste0("resultsrun_",j,".csv"))
  
  results[1,j,1]=  tmp[1,2]
  results[1,j,2]=  tmp[1,3]
  results[1,j,3] =   tmp[1,4]
  
  if(tmp[2,2] == 0)
    print(j)
  
  results[2,j,1]=  tmp[2,2]
  results[2,j,2]=   tmp[2,3]
  results[2,j,3] =   tmp[2,4]
  
}

#make the joint summary of the runs, including min, max and medians of the performance metrics
summary.results=array(data = NA,dim = c(2,15))

for(i in 1:2){
for(j in 1:5)
{
  summary.results[i,(j-1)*3+1]=min(results[i,,j])
  summary.results[i,(j-1)*3+2]=median(results[i,,j])
  summary.results[i,(j-1)*3+3]=max(results[i,,j])
}
}
summary.results=as.data.frame(summary.results)


#featgmj = hash()

simplifyposteriors<-function(X,post,th=0.0001,thf=0.1,y = "shares")
{
  posteriors = (cbind(as.character(post[,1]),as.numeric(as.character(post[,2]))))
  
  rhash<-hash()
  for(i in 1:length(posteriors[,1]))
  {
    
    
    expr<-posteriors[i,1]
    print(expr)
    res<-model.matrix(data=X,object = as.formula(paste0(y,"~",expr)))
    ress<-c(stri_flatten(round(sum(res[,2]),digits = 4),collapse = ""),stri_flatten(res[,2],collapse = ""),posteriors[i,2],expr)
    if(!((ress[1] %in% keys(rhash))))
      rhash[[ress[1]]]<-ress
    else
    {
      if(ress[1] %in% keys(rhash))
      {new
        rhash[[ress[1]]][3]<- (as.numeric(rhash[[ress[1]]][3]) + as.numeric(ress[3]))
        if(stri_length(rhash[[ress[1]]][4])>stri_length(expr))
          rhash[[ress[1]]][4]<-expr
      }
    }
    
  }
  res<-as.data.frame(t(values(rhash)[c(3,4),]))
  res$V1<-as.numeric(as.character(res$V1))
  res<-res[which(res$V1>thf),]
  res<-res[order(res$V1, decreasing = T),]
  clear(rhash)
  rm(rhash)
  #res[which(res[,1]>1),1]<-1
  colnames(res)<-c("posterior","tree")
  
  row.names(res) = 1:length(res$posterior)
  return(res)
}




featgmj = hash()

for(j in 1:100)
{
  tmpf = read.csv(paste0("posteriorshell_",j,".csv"))
  #tmp = simplifyposteriors(X = data.example,post =tmpf,y = "shares")
  for(feat in as.character(tmpf$V1))
  {
    if(!has.key(hash = featgmj,key =  feat ))
    {
      featgmj[[feat]] = as.numeric(1)
    } else{
      
      featgmj[[feat]] =as.numeric(featgmj[[feat]]) + 1
    }
  }
}

tmp = simplifyposteriors(X = data.example,post =as.data.frame(cbind(keys(featgmj),as.numeric(values(featgmj)))),y = "shares")

write.csv(x =tmp,row.names = T,file = "abalonefeat.csv")
#print(paste0("end simulation ",j))


