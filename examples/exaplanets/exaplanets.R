data <- read.table(file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/exaplanets/oec.csv",sep = ",",header = T,fill=TRUE)
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package2.r")


nas<-unique(which(is.na(data[,-c(1,8,9,10,11,12,13,14,15,16,17,18,19,24)]), arr.ind=TRUE)[,1])
length(nas)

hist(log((data[,3])),breaks = 1000)



train<-data[-nas,-c(1,8,9,10,11,12,13,14,15,16,17,18,19,24)]
test  <- data[nas, -c(1,8,9,10,11,12,13,14,15,16,17,18,19,24)]
data.example <- as.data.frame(train[1,])


estimate.gamma.cpen <- function(formula, data,r = 1.0/223.0,logn=log(223.0),relat=c("cosi","sigmoid","tanh","atan","sini"))
{
  fparam<-NULL
  fmla.proc<-as.character(formula)[2:3]
  fobserved <- fmla.proc[1]
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]<-stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam <-stri_split_fixed(str = fmla.proc[2],pattern = "+I",omit_empty = F)[[1]]
  sj<-(stri_count_fixed(str = fparam, pattern = "*"))
  sj<-sj+(stri_count_fixed(str = fparam, pattern = "+"))
  for(rel in relat)
  sj<-sj+(stri_count_fixed(str = fparam, pattern = rel))
  sj<-sj+1
  tryCatch(capture.output({
  out <- glm(formula = formula,data = data, family = gaussian)
  # 1 for aic, 2 bic prior, else g.prior

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

pairs(train)
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
results<-array(0,dim=c(26,10))
train$PlanetaryDensJpt<-train[,2]/(train[,3]^3)
#names(train)[11]="PlanetaryDensJpt"
resp=2

#train[,2]<-log(train[,2])
data.example <- as.data.frame(train)

for(i in 1:length(resp))
{
  print(paste0("Model ",i," for ",names(test)[resp][i]))
  formula1 = as.formula(paste(colnames(train)[resp[i]],"~ 1 +",paste0(colnames(train)[-resp[i]],collapse = "+")))
  print(formula1)
  
  #print(estimate.gamma.cpen (formula=formula1, data=train))
  

  res = runemjmcmc(formula = formula1,data = data.example,estimator =estimate.gamma.cpen,estimator.args =  list(data = data.example),recalc_margin = 249, save.beta = F,interact = T,outgraphs=T,relations = c("","cosi","sigmoid","tanh","atan","sini"),relations.prob =c(0.9,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=2,mutation_rate = 250,last.mutation=5000, max.tree.size = 5, Nvars.max =15,p.allow.replace=0.9,p.allow.tree=0.2,p.nor=0.9,p.and = 0.9),n.models = 10000,unique =T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
    max.N.glob=as.integer(10),
    min.N.glob=as.integer(5),
    max.N=as.integer(3),
    min.N=as.integer(1),
    printable = F))

  
  ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
  ppp$p.post
  length(ppp$m.post)
  
  
  results[1,2*i-1]<-paste0("Parameter for ",colnames(train)[resp[i]])
  results[1,2*i]<-"Post.prob"
  idi<-order(ppp$p.post,decreasing =T)
  
  results[2:(length(mySearch$fparam[idi])+1),2*i-1]<-mySearch$fparam[idi]
  results[2:(length(mySearch$fparam[idi])+1),2*i]<-ppp$p.post[idi]
  
  
}


par(mar = c(15,0,4,0) + 4.1)
for(i in 1:length(resp))
{
  barplot(height=as.numeric(results[2:(length(mySearch$fparam[idi])+1),2*i]),density = 46,border="black",main = paste0("Marginal Inclusion (",names(test)[resp][i],"~ 1 + X'b)"),ylab="Probability",names.arg = results[2:21,2*i-1],las=3)
}

print(estimate.gamma.cpen (formula=PlanetaryMassJpt ~ 1 + I((I(sini(RadiusJpt)*(((I((RadiusJpt)*(((I((PlanetaryDensJpt)*((RadiusJpt))))))))))))), data=data.example))

print(estimate.gamma.cpen (formula=PlanetaryMassJpt ~ 1 + I((I((sini(I((RadiusJpt)*(((I(-(RadiusJpt)*((PlanetaryDensJpt)))))))))*((RadiusJpt))))), data=train))

print(estimate.gamma.cpen (formula=PlanetaryMassJpt ~ 1 + I(RadiusJpt*RadiusJpt*PlanetaryDensJpt) + I(RadiusJpt*PlanetaryDensJpt)+I(PlanetaryDensJpt), data=train))

print(estimate.gamma.cpen (formula=PlanetaryMassJpt ~ 1 + I((-(RadiusJpt))*(-RadiusJpt)*((-RadiusJpt))*(-PlanetaryDensJpt))+RadiusJpt, data=train))
print(estimate.gamma.cpen (formula=PlanetaryMassJpt ~ 1 + I((-(RadiusJpt))*(-RadiusJpt)*((-RadiusJpt))*(-PlanetaryDensJpt))+RadiusJpt + PeriodDays+SemiMajorAxisAU, data=train))


