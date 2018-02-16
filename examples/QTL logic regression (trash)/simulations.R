#scenario 1

M=100
resa<-array(data = 0,dim = c(16,M*3))
NM= 1000
post.popul <- array(0,M)
max.popul <- array(0,M)
set.seed(040590)
X1<- as.data.frame(array(data = rbinom(n = 50*100,size = 1,prob = runif(n = 50*100,0,1)),dim = c(100,50)))
Y1<-rnorm(n = 100,mean = 1+0.7*(X1$V1*X1$V4) + 0.8896846*(X1$V8*X1$V11)+1.434573*(X1$V5*X1$V9),sd = 1)
X1$Y1<-Y1

for(ii in 1:M)
{
  set.seed(ii)
  print(ii)
  #GMJMCMC
  system.time({

    formula1 = as.formula(paste(colnames(X1)[51],"~ 1 +",paste0(colnames(X1)[-c(51)],collapse = "+")))

    res = runemjmcmc(formula = formula1,data = X1,estimator = estimate.logic.lm,estimator.args =  list(data = data.example,n = 100, m = 50),recalc_margin = 250, save.beta = F,interact = T,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=1,mutation_rate = 300,last.mutation = 5000, max.tree.size = 4, Nvars.max =15,p.allow.replace=0.9,p.allow.tree=0.2,p.nor=0,p.and = 0.9),n.models = 10000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,print.freq = 100,advanced.param = list(
      max.N.glob=as.integer(10),
      min.N.glob=as.integer(5),
      max.N=as.integer(3),
      min.N=as.integer(1),
      printable = F))
  })



  cterm<-max(values(hashStat)[1,])
  max.popul[ii]<-cterm
  ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
  mySearch$g.results[,]
  post.popul[ii]<-sum(exp(values(hashStat)[1,][1:NM]-cterm),na.rm = T)
  resa[,ii*3-2]<-c(mySearch$fparam,"Post.Gen.Max")
  resa[,ii*3-1]<-c(ppp$p.post,cterm)
  resa[,ii*3]<-rep(post.popul[ii],length(ppp$p.post)+1)
  gc()
  print(cterm)

}
ml.max<-max(max.popul)
post.popul<-post.popul*exp(-ml.max+max.popul)
p.gen.post<-(post.popul)/(sum(post.popul))
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
        print(paste0(ii,"  and ",jj))
        if(as.integer(has.key(hash = hfinal,key =resa[jj,ii*3-2]))==0)
        {
          hfinal[[resa[jj,ii*3-2]]]<-as.numeric(resa[jj,ii*3])
        }else
        {
          hfinal[[resa[jj,ii*3-2]]]<-hfinal[[resa[jj,ii*3-2]]]+as.numeric(resa[jj,ii*3])
        }
      }

    }
  }
}
write.csv(x = values(hfinal),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/sim1post.csv")
clear(hfinal)
rm(hfinal)
write.csv(x = (resa),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/sim1.csv")
write.csv(x =X1,file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/X1.csv")


# now trying out the logic regression
# 
# library(LogicReg)
# 
# system.time({
# lmod <- logreg(
#   resp=X1$Y1,bin=X1[,-51], 
#   type=2, select = 7, 
#   ntrees=5,
#   nleaves =5*5,
#   mc.control=logreg.mc.control(nburn=500, niter=250000, hyperpars=log(2),output = 3,update = -1),
#   tree.control=logreg.tree.control(treesize=5
#                                    ,opers=2)
#   )
# })
# 
# lmod$size
# iter<-sum(lmod$size[,2])
# # marginal posteriors of covariates in a tree
# lmod$single[1]/iter
# lmod$single[4]/iter
# lmod$single[8]/iter
# lmod$single[11]/iter
# lmod$single[5]/iter
# lmod$single[9]/iter
# max(lmod$single)/iter
# min(lmod$single)/iter
# which(lmod$single == max(lmod$single))
# order(lmod$single)
# plot(lmod$single/iter)
# # pairwise probability of the pairs of covariates in a tree
# lmod$double[4,1]/iter 
# lmod$double[11,8]/iter 
# lmod$double[9,5]/iter
# max(lmod$double)/iter
# min(lmod$double)/iter
# plot(y = lmod$double/iter)
# which(lmod$double == max(lmod$double))
# #lmod$triple # tupples of them
# max(lmod$triple)/iter

#scenario 2

M=100
resa<-array(data = 0,dim = c(16,M*3))
NM= 1000
post.popul <- array(0,M)
max.popul <- array(0,M)
set.seed(040590)
X2<- as.data.frame(array(data = rbinom(n = 50*100,size = 1,prob = runif(n = 50*100,0,1)),dim = c(100,50)))
Y2<-rnorm(n = 100,mean = 1+3.5*(X2$V1*X2$V4) + 6.6*(X2$V8*X2$V11)+ 8.9*(X2$V5*X2$V9),sd = 1)
X2$Y2<-Y2
for(ii in 1:M)
{
  set.seed(ii)
  print(ii)
  #GMJMCMC
  system.time({
    
    formula1 = as.formula(paste(colnames(X2)[51],"~ 1 +",paste0(colnames(X2)[-c(51)],collapse = "+")))
    
    res = runemjmcmc(formula = formula1,data = X2,estimator = estimate.logic.lm,estimator.args =  list(data = data.example,n = 100, m = 50),recalc_margin = 250, save.beta = F,interact = T,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=1,mutation_rate = 300,last.mutation = 5000, max.tree.size = 4, Nvars.max =15,p.allow.replace=0.9,p.allow.tree=0.2,p.nor=0,p.and = 0.9),n.models = 10000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,print.freq = 100,advanced.param = list(
      max.N.glob=as.integer(10),
      min.N.glob=as.integer(5),
      max.N=as.integer(3),
      min.N=as.integer(1),
      printable = F))
  })
  
  
  
  cterm<-max(values(hashStat)[1,])
  max.popul[ii]<-cterm
  ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
  mySearch$g.results[,]
  post.popul[ii]<-sum(exp(values(hashStat)[1,][1:NM]-cterm),na.rm = T)
  resa[,ii*3-2]<-c(mySearch$fparam,"Post.Gen.Max")
  resa[,ii*3-1]<-c(ppp$p.post,cterm)
  resa[,ii*3]<-rep(post.popul[ii],length(ppp$p.post)+1)
  gc()
  print(cterm)
  
}

ml.max<-max(max.popul)
post.popul<-post.popul*exp(-ml.max+max.popul)
p.gen.post<-(post.popul)/(sum(post.popul))
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
        print(paste0(ii,"  and ",jj))
        if(as.integer(has.key(hash = hfinal,key =resa[jj,ii*3-2]))==0)
        {
          hfinal[[resa[jj,ii*3-2]]]<-as.numeric(resa[jj,ii*3])
        }else
        {
          hfinal[[resa[jj,ii*3-2]]]<-hfinal[[resa[jj,ii*3-2]]]+as.numeric(resa[jj,ii*3])
        }
      }
      
    }
  }
}
write.csv(x = values(hfinal),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/sim2post.csv")
clear(hfinal)
rm(hfinal)
write.csv(x = (resa),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/sim2.csv")
write.csv(x =X2,file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/X2.csv")


#scenario 3

M=100
resa<-array(data = 0,dim = c(21,M*3))
NM= 1000
post.popul <- array(0,M)
max.popul <- array(0,M)
set.seed(040590)
X3<- as.data.frame(array(data = rbinom(n = 50*100,size = 1,prob = runif(n = 50*100,0,1)),dim = c(100,50)))
Y3<-rnorm(n = 100,mean = 1+7*(X3$V4*X3$V17*X3$V30*X3$V10) + 9*(X3$V7*X3$V20*X3$V12)+ 3.5*(X3$V9*X3$V2)+1.5*X3$V37,sd = 1)
X3$Y3<-Y3

for(ii in 1:M)
{

set.seed(ii)
print(ii)

#GMJMCMC
system.time({
  
  formula1 = as.formula(paste(colnames(X3)[51],"~ 1 +",paste0(colnames(X3)[-c(51)],collapse = "+")))
  
  res = runemjmcmc(formula = formula1,data = X3,estimator = estimate.logic.lm,estimator.args =  list(data = data.example,n = 100, m = 50),recalc_margin = 250, save.beta = F,interact = T,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=1,mutation_rate = 300,last.mutation = 5000, max.tree.size = 4, Nvars.max =20,p.allow.replace=0.9,p.allow.tree=0.2,p.nor=0,p.and = 0.9),n.models = 10000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,print.freq = 100,advanced.param = list(
    max.N.glob=as.integer(10),
    min.N.glob=as.integer(5),
    max.N=as.integer(3),
    min.N=as.integer(1),
    printable = F))
})



cterm<-max(values(hashStat)[1,])
max.popul[ii]<-cterm
ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
mySearch$g.results[,]
post.popul[ii]<-sum(exp(values(hashStat)[1,][1:NM]-cterm),na.rm = T)
resa[,ii*3-2]<-c(mySearch$fparam,"Post.Gen.Max")
resa[,ii*3-1]<-c(ppp$p.post,cterm)
resa[,ii*3]<-rep(post.popul[ii],length(ppp$p.post)+1)
gc()
print(cterm)
}

ml.max<-max(max.popul)
post.popul<-post.popul*exp(-ml.max+max.popul)
p.gen.post<-(post.popul)/(sum(post.popul))
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
        print(paste0(ii,"  and ",jj))
        if(as.integer(has.key(hash = hfinal,key =resa[jj,ii*3-2]))==0)
        {
          hfinal[[resa[jj,ii*3-2]]]<-as.numeric(resa[jj,ii*3])
        }else
        {
          hfinal[[resa[jj,ii*3-2]]]<-hfinal[[resa[jj,ii*3-2]]]+as.numeric(resa[jj,ii*3])
        }
      }

    }
  }
}
write.csv(x = values(hfinal),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/sim3post.csv")
clear(hfinal)
rm(hfinal)
write.csv(x = (resa),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/sim3.csv")
write.csv(x =X3,file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/X3.csv")

#scenario 4
M=100
resa<-array(data = 0,dim = c(26,M*3))
NM= 1000
post.popul <- array(0,M)
max.popul <- array(0,M)
set.seed(040590)
X4<- as.data.frame(array(data = rbinom(n = 50*100,size = 1,prob = runif(n = 50*100,0,1)),dim = c(100,50)))
Y4<-rnorm(n = 100,mean = 1+7*(X4$V4*X4$V17*X4$V30*X4$V10)+7*(as.integer((X4$V50*X4$V19+X4$V13*X4$V11)/2)) + 9*(X4$V37*X4$V20*X4$V12)+ 7*(X4$V1*X4$V27*X4$V3)
          +3.5*(X4$V9*X4$V2) + 6.6*(X4$V21*X4$V18) + 1.5*X4$V7 + 1.5*X4$V8,sd = 1)
X4$Y4<-Y4
for(ii in 1:M)
{
  set.seed(ii)
  print(ii)
  
  #GMJMCMC
  system.time({
    
    formula1 = as.formula(paste(colnames(X4)[51],"~ 1 +",paste0(colnames(X4)[-c(51)],collapse = "+")))
    
    res = runemjmcmc(formula = formula1,data = X4,estimator = estimate.logic.lm,estimator.args =  list(data = data.example,n = 100, m = 50),recalc_margin = 249, save.beta = F,interact = T,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=1,mutation_rate = 250,last.mutation = 5000, max.tree.size = 4, Nvars.max =25,p.allow.replace=0.9,p.allow.tree=0.2,p.nor=0,p.and = 0.9),n.models = 10000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,print.freq = 100,advanced.param = list(
      max.N.glob=as.integer(10),
      min.N.glob=as.integer(5),
      max.N=as.integer(3),
      min.N=as.integer(1),
      printable = F))
  })
  
  
  
  cterm<-max(values(hashStat)[1,])
  max.popul[ii]<-cterm
  ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
  mySearch$g.results[,]
  post.popul[ii]<-sum(exp(values(hashStat)[1,][1:NM]-cterm),na.rm = T)
  resa[,ii*3-2]<-c(mySearch$fparam,"Post.Gen.Max")
  resa[,ii*3-1]<-c(ppp$p.post,cterm)
  resa[,ii*3]<-rep(post.popul[ii],length(ppp$p.post)+1)
  gc()
  print(cterm)
}

ml.max<-max(max.popul)
post.popul<-post.popul*exp(-ml.max+max.popul)
p.gen.post<-(post.popul)/(sum(post.popul))
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
        print(paste0(ii,"  and ",jj))
        if(as.integer(has.key(hash = hfinal,key =resa[jj,ii*3-2]))==0)
        {
          hfinal[[resa[jj,ii*3-2]]]<-as.numeric(resa[jj,ii*3])
        }else
        {
          hfinal[[resa[jj,ii*3-2]]]<-hfinal[[resa[jj,ii*3-2]]]+as.numeric(resa[jj,ii*3])
        }
      }
      
    }
  }
}
write.csv(x = values(hfinal),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/sim4post.csv")
clear(hfinal)
rm(hfinal)
write.csv(x = (resa),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/sim4.csv")
write.csv(x =X4,file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/X4.csv")






#scenario 5
M=100
resa<-array(data = 0,dim = c(31,M*3))
NM= 1000
post.popul <- array(0,M)
max.popul <- array(0,M)
set.seed(040590)
X5<- as.data.frame(array(data = rbinom(n = 50*100,size = 1,prob = runif(n = 50*100,0,1)),dim = c(100,50)))
Y5<-rnorm(n = 100,mean = 1+7*(X5$V20*X5$V19*X5$V16*X5$V13)+9.8*(X5$V12*X5$V8*X5$V5) + 9*(X5$V2*X5$V28*X5$V22)+ 8.9*(X5$V10*X5$V30)
          +3.5*(X5$V48*X5$V1) + 6.6*(X5$V50*X5$V3) + 1.5*X5$V32 + 1.5*X5$V7 + 1.5*X5$V24 + 1.5*X5$V11 + 1.5*X5$V4,sd = 1)
X5$Y5<-Y5
for(ii in 1:M)
{
  set.seed(ii)
  print(ii)
  
  #GMJMCMC
  system.time({
    
    formula1 = as.formula(paste(colnames(X5)[51],"~ 1 +",paste0(colnames(X5)[-c(51)],collapse = "+")))
    
    res = runemjmcmc(formula = formula1,data = X5,estimator = estimate.logic.lm,estimator.args =  list(data = data.example,n = 100, m = 50),recalc_margin = 250, save.beta = F,interact = T,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=1,mutation_rate = 300,last.mutation = 1000, max.tree.size = 4, Nvars.max =30,p.allow.replace=0.9,p.allow.tree=0.2,p.nor=0,p.and = 0.9),n.models = 12500,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,print.freq = 100,advanced.param = list(
      max.N.glob=as.integer(10),
      min.N.glob=as.integer(5),
      max.N=as.integer(3),
      min.N=as.integer(1),
      printable = F))
  })
  
  
  
  cterm<-max(values(hashStat)[1,])
  max.popul[ii]<-cterm
  ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
  mySearch$g.results[,]
  post.popul[ii]<-sum(exp(values(hashStat)[1,][1:NM]-cterm),na.rm = T)
  resa[,ii*3-2]<-c(mySearch$fparam,"Post.Gen.Max")
  resa[,ii*3-1]<-c(ppp$p.post,cterm)
  resa[,ii*3]<-rep(post.popul[ii],length(ppp$p.post)+1)
  gc()
  print(cterm)
}

ml.max<-max(max.popul)
post.popul<-post.popul*exp(-ml.max+max.popul)
p.gen.post<-(post.popul)/(sum(post.popul))
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
        print(paste0(ii,"  and ",jj))
        if(as.integer(has.key(hash = hfinal,key =resa[jj,ii*3-2]))==0)
        {
          hfinal[[resa[jj,ii*3-2]]]<-as.numeric(resa[jj,ii*3])
        }else
        {
          hfinal[[resa[jj,ii*3-2]]]<-hfinal[[resa[jj,ii*3-2]]]+as.numeric(resa[jj,ii*3])
        }
      }
      
    }
  }
}
write.csv(x = values(hfinal),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/sim5post.csv")
clear(hfinal)
rm(hfinal)
write.csv(x = (resa),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/sim5.csv")

write.csv(x =X5,file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/X5.csv")

