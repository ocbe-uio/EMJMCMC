#read the most recent stable version of the package
#source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")
#install.packages("https://github.com/aliaksah/EMJMCMC2016/blob/master/EMJMCMC_1.4.2_R_x86_64-pc-linux-gnu.tar.gz?raw=true", 
#                 repos = NULL, type="source")
# load the EMJMCMC package
library(EMJMCMC)
set.seed(040590)
#make sure that you are using Mac Os or Linux (mclapply is currently not supported for Windows unless some mclapply hack function for Windows is preloaded in your R session)



## Construct a binary correlation matrix for M = 50 cariables
M = 50
m = clusterGeneration::rcorrmatrix(M,alphad=2.5) 
#print the highest 10 correlations in the data
print(unique(sort(abs(m),decreasing = T))[1:10])
#print the lowest 10 correlations in the data
print(unique(sort(abs(m),decreasing = F))[1:10])
#simulate 1000 binary variables with a given correlation's structure
X = bindata::rmvbin(1000, margprob = rep(0.5,M), bincorr = m)

print(unique(sort(abs(cor(X)),decreasing = T))[1:10])

melted_cormat = reshape2::melt(cor(as.data.frame(X)))
ggplot2::ggplot(data = melted_cormat,
                ggplot2::aes(x=Var1, y=Var2, fill=value)) +
  ggplot2::geom_tile() +
  ggplot2::labs(fill = "Corr") +
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.title.y =  ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 legend.title = ggplot2::element_text( size = 26),
                 legend.text = ggplot2::element_text( size = 17)
  )

#simulate Gaussian responses from a model with two-way interactions
#which is considered in S.4 of the paper
X = data.frame(X)
Y=rnorm(n = 1000,mean = 1+1.43*(X$X5*X$X9)+ 0.89*(X$X8*X$X11)+0.7*(X$X1*X$X4),sd = 1)
X$Y=Y

#specify the initial formula
formula1 = as.formula(paste(colnames(X)[M+1],"~ 1 +",paste0(colnames(X)[-c(M+1)],collapse = "+")))
df = as.data.frame(X)

# run the inference with robust g prior
res4G = LogicRegr(formula = formula1,data = df ,family = "Gaussian",prior = "G",report.level = 0.5,d = 10,cmax = 2,kmax = 5,p.and = 0.9,p.not = 0.1,p.surv = 0.2,ncores = 5)
print(res4G$feat.stat)

# run the inference with Jeffrey's prior
res4J = LogicRegr(formula = formula1,data = df,family = "Gaussian",prior = "J",report.level = 0.5,d = 10,cmax = 2,kmax = 5,p.and = 0.9,p.not = 0.1,p.surv = 0.2,ncores = 5)
print(res4J$feat.stat)



Xp = X
Xp$age = rpois(1000,lambda = 34)
Y=rnorm(n = 1000,mean = 1+0.7*(Xp$X1*Xp$X4) + 0.8896846*(Xp$X8*Xp$X11)+1.434573*(Xp$X5*Xp$X9) + 2*Xp$age,sd = 1)
#Xp$age = NULL
Xp$Y=Y

teid =  sample.int(size =100,n = 1000,replace = F)
test = Xp[teid,]
train = Xp[-teid,]

sum(test$age)

# specify the initial formula
formula1 = as.formula(paste("Y ~ 1 +",paste0(colnames(test)[-c(51,52)],collapse = "+")))
# specify the link function
g = function(x) x
res.alt = pinferunemjmcmc(n.cores = 30, report.level =  0.2, num.mod.best = 100,simplify = T,predict = T,test.data = test,link.function = g,runemjmcmc.params = list(formula = formula1,latnames = c("I(age)"),data = train,estimator = estimate.logic.lm.jef,estimator.args = list(data = train, n = dim(train)[1], m =stri_count_fixed(as.character(formula1)[3],"+"),k.max = 15),recalc_margin = 249, save.beta = T,interact = T,outgraphs=F,relations=c("sin","cos"),relations.prob =c(0.5,0.5),interact.param=list(allow_offsprings=1,mutation_rate = 250,last.mutation=10000, max.tree.size = 5, Nvars.max =15,p.allow.replace=0.9,p.allow.tree=0.01,p.nor=0.01,p.and = 0.9),n.models = 10000,unique =T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
  max.N.glob=10,
  min.N.glob=5,
  max.N=3,
  min.N=1,
  printable = F)))

print(res.alt$feat.stat)
print(sqrt(mean((res.alt$predictions-test$Y)^2)))
print(mean(abs((res.alt$predictions-test$Y))))

res.alt = pinferunemjmcmc(n.cores = 30, report.level =  0.2, num.mod.best = 100,simplify = T,predict = T,test.data = as.data.frame(test),link.function = g,runemjmcmc.params = list(formula = formula1,latnames = c("I(age)"),data = as.data.frame(train),estimator = estimate.logic.lm.jef,estimator.args = list(data = as.data.frame(train), n = dim(as.data.frame(train))[1], m =stri_count_fixed(as.character(formula1)[3],"+"),k.max = 15),recalc_margin = 249, save.beta = T,interact = T,outgraphs=F,relations=c("to25","expi","logi","to35","troot","sigmoid"),relations.prob =c(0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=1,mutation_rate = 250,last.mutation=10000, max.tree.size = 5, Nvars.max =15,p.allow.replace=0.9,p.allow.tree=0.01,p.nor=0.01,p.and = 0.9),n.models = 10000,unique =T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
                            max.N.glob=as.integer(10),
                            min.N.glob=as.integer(5),
                            max.N=as.integer(3),
                            min.N=as.integer(1),
                            printable = F)))





library(HDeconometrics)
ridge = ic.glmnet(x = train[,-51],y=train$Y,family = "gaussian",alpha = 0)
predict.ridge = predict(ridge$glmnet,newx = as.matrix(test[,-51]),type = "response")[,which(ridge$glmnet$lambda == ridge$lambda)]
print(sqrt(mean((predict.ridge-test$Y)^2)))
print(mean(abs((predict.ridge-test$Y))))

tmean = 1+2*test$age+0.7*(test$X1*test$X4) + 0.8896846*(test$X8*test$X11)+1.434573*(test$X5*test$X9)
print(sqrt(mean((tmean -test$Y)^2)))
print(mean(abs((tmean -test$Y))))


plot(test$Y,type = "p",pch = 16)
points(tmean, col = "green", pch = 19)
points(res.alt$predictions,col="red",pch = 3)
points(predict.ridge,col = "blue",pch = 4)


#library(aRxiv)
set.seed(040590)
et = arima.sim(list(order=c(1,0,0), ar=0.7), n=1000)
Xp$Y = Xp$Y + 9*et



#define the function estimating parameters of a given Gaussian logic regression with Jeffrey's prior
estimate.logic.inla = function(formula,family = "gaussian", data, n = 10000, m = 50 , r = 1,k.max=15, bias = 0)
{
  
  
  print(formula)
  if(is.na(formula)||is.null(formula))
    return(list(mlik =  -10000 + rnorm(1,0,1),waic =10000+ rnorm(1,0,1) , dic =  10000+ rnorm(1,0,1),summary.fixed =list(mean = 1)))


  fmla.proc=as.character(formula)[2:3]
  fobserved = fmla.proc[1]
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam =stri_split_fixed(str = fmla.proc[2],pattern = "+",omit_empty = F)[[1]]
  
  p = length(fparam)
  if(p>k.max)
  {
    return(list(mlik = -10000+ rnorm(1,0,1),waic = 10000+ rnorm(1,0,1) , dic =  10000+ rnorm(1,0,1),summary.fixed =list(mean = 0)))
  }
  
  out<-NULL
  capture.output({tryCatch(capture.output({
    out <-inla(family = family,silent = 2L,data = data,formula = formula,control.compute = list(dic = TRUE, waic = TRUE, mlik = TRUE))
  }))})

  # use dic and aic as bic and aic correspondinly
  
  coef<-out$summary.fixed$mode
  #coef[1]<-coef[1]+out$summary.hyperpar$mode[1]
  
  sj=(stri_count_fixed(str = fparam, pattern = "&"))
  sj=sj+(stri_count_fixed(str = fparam, pattern = "|"))
  sj=sj+1
  Jprior = prod(factorial(sj)/((m^sj)*2^(2*sj-2)))
  
  if(is.null(out))
    return(list(mlik = -10000+log(r)*(sj),waic =  10000 , dic = 10000, summary.fixed =list(mean = NULL)))
  
  if(length(out$waic[1]$waic)==0)
    return(list(mlik = -10000+log(r)*(sj),waic =  10000 , dic = 10000, summary.fixed =list(mean = NULL)))
  
  mlik = (out$mlik[1]+log(Jprior) + p*log(r)+n) + bias
  if(mlik==-Inf)
    mlik = -10000 + rnorm(1,0,1)
  return(list(mlik = mlik,waic =  out$waic[1]$waic , dic = out$dic[1]$dic, summary.fixed =list(mean = coef)))
}

estimate.logic.inla(data = Xp,formula = Y ~ age+ 1 + X1 + X2 + X3 + X4 + X5)

Xp$pos = 1:length(Xp$X1)
Xp$pos1 = 1:length(Xp$X1)
Xp$pos2 = 1:length(Xp$X1)
Xp$pos3 = 1:length(Xp$X1)

#estimate.logic.inla(data = Xp,formula = Y ~ age+ 1 + X1 + X2 + X3 + X4 + X5 + f(Xp$pos,model="ar1") + f(data.example$pos1,model="rw1")+f(data.example$pos2,model="iid"))

data.example = Xp
res.mixed = pinferunemjmcmc(n.cores = 15, report.level =  0.2, num.mod.best = 1000,simplify = T,predict = F,test.data = as.data.frame(test),link.function = g,runemjmcmc.params = list(formula = formula1,latnames = c("I(age)","f(data.example$pos,model=\"ar1\")","f(data.example$pos1,model=\"rw1\")","f(data.example$pos2,model=\"iid\")","f(data.example$pos3,model=\"ou\")"),data = data.example,estimator = estimate.logic.inla,estimator.args = list(data = data.example, n = 1000, r = 1, m =stri_count_fixed(as.character(formula1)[3],"+"),k.max =20),recalc_margin = 249, save.beta = T,interact = T,outgraphs=F,relations=c("to25","expi","logi","to35","troot","sigmoid"),relations.prob =c(0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=1,mutation_rate = 250,last.mutation=10000, max.tree.size = 5, Nvars.max =20,p.allow.replace=0.9,p.allow.tree=0.01,p.nor=0.01,p.and = 0.9),n.models = 5000,unique =T,max.cpu = 4,max.cpu.glob = 4,create.table = F,presearch = F, create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
  max.N.glob=as.integer(10),
  min.N.glob=as.integer(5),
  max.N=as.integer(3),
  min.N=as.integer(1),
  printable = F)))

res.mixed$threads.stats[[1]]$fparam
res.mixed$threads.stats[[1]]$p.post

res.mixed$allposteriors


#now use qtl package to generate data in a more realistic scenario


library(qtl)
set.seed(040590)  #you can change this of course

# simulate 5 autosomes, each with 10 equally spaced markers 
# with different chromosomal lengths 
map = sim.map(c(100,90,80,60,40), M/5, include.x=FALSE, eq.spacing=TRUE)
plotMap(map)

# The longer the chromosomal length the less correlated are markers
# (at least theoretically)

# Now simulate data from a backcross design using the sim.cross function
n.ind = 1000   #sample size

simbc = sim.cross(map, type="bc", n.ind=n.ind,
                    model=rbind(c(1,45,1), c(5,20,1), c(5,50,1)))

# The relevant genotype data is in the structure geno 
str(simbc$geno)   #it is a bit complicated

# Get an X matrix for each chromosome
X.list = vector("list", 5)
for (chr in 1:5){
  X.list[[chr]] = pull.geno(simbc, chr)
}

#Check the correlations between markers within the same chromosome
lapply(X.list, cor)


#Create one large X matrix which you can the use to make your own 
# simulations of Y data with a logic regression model
X = cbind(X.list[[1]],X.list[[2]],X.list[[3]],X.list[[4]],X.list[[5]])-1
#permute elements of X

#plot the correlation's structure 
melted_cormat = reshape2::melt(cor(as.matrix(X)))
ggplot2::ggplot(data = melted_cormat,
                ggplot2::aes(x=Var1, y=Var2, fill=value)) +
  ggplot2::geom_tile() +
  ggplot2::labs(fill = "Corr") +
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.title.y =  ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 legend.title = ggplot2::element_text( size = 26),
                 legend.text = ggplot2::element_text( size = 17)
  )
X = X[,sample.int(n = M,size = M,replace = F)]


X = data.frame(X)
Y=rnorm(n = 1000,mean = 1+0.7*(X[,1]*X[,4]) + 0.8896846*(X[,8]*X[,11])+1.434573*(X[,5]*X[,9]),sd = 1)
X$Y=Y

print(paste0("L1 = ", paste(names(X)[c(1,4)],collapse = ", ",sep = "")))
print(paste0("L2 = ", paste(names(X)[c(8,11)],collapse = ", ",sep = "")))
print(paste0("L3 = ", paste(names(X)[c(5,9)],collapse = ", ",sep = "")))

#specify the initial formula
formula1 = as.formula(paste(colnames(X)[M+1],"~ 1 +",paste0(colnames(X)[-c(M+1)],collapse = "+")))
data.example = as.data.frame(X)

#run the inference with robust g prior
res4G = LogicRegr(formula = formula1,data = data.example,family = "Gaussian",prior = "G",report.level = 0.5,d = 15,cmax = 2,kmax = 15,p.and = 0.9,p.not = 0.01,p.surv = 0.2,ncores = 32)
print(res4G$feat.stat)
#run the inference with Jeffrey's prior
res4J = LogicRegr(formula = formula1,data = data.example,family = "Gaussian",prior = "J",report.level = 0.5,d = 15,cmax = 2,kmax = 15,p.and = 0.9,p.not = 0.01,p.surv = 0.2,ncores = 32)
print(res4J$feat.stat)

