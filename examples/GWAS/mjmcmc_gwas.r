#library(dplyr)
options("exppression" = 5000)
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")


pheno<-read.csv("/home/michaelh/SIMULATION_paper/data_S2_nocausal_5402/pimass/data.recode.pheno_1.txt",header = F)
geno<-t(read.csv("/home/michaelh/SIMULATION_paper/data_S2_nocausal_5402/pimass/data.recode.mean.geno.txt",header = F,stringsAsFactors = F))
names<-geno[1,]
geno<-as.data.frame(geno[-c(1,2,3),])
names(geno)<-names
geno$Y<-pheno$V1
geno <- as.data.frame(mclapply(geno, as.numeric))

cors<-cor(geno$Y,geno[,1:24602])
sum(abs(cors) > 0.05)


cov.names<-names[which(cors>0.03)]
sum<-summary(lm(as.formula(paste0("Y~1+",paste(cov.names,collapse = "+"))),data = geno))
cov.names<-names(sum$coefficients[-1,4])

formula1 <- as.formula(paste0("Y~1+",paste(cov.names,collapse = "+")))

estimate.lm.MBIC2 <- function(formula, data, n = 5402, m = 24602, c = 4,u=50)
{
  size<-stri_count_fixed(str = as.character(formula)[3],pattern = "+")
 
  if(size>u)
  {
    return(list(mlik = -20000,waic = 50000 , dic =  50000,summary.fixed =list(mean = array(0,size+1))))
  }
    
  out <- lm(formula = formula,data = data)
  logmarglik <- (2*logLik(out) - out$rank*log(m*m*n/c) + 2*log(factorial(out$rank)))/2
  # use dic and aic as bic and aic correspondinly
  return(list(mlik = logmarglik,waic = AIC(out) , dic =  BIC(out),summary.fixed =list(mean = coef(out))))
  
}


vect<-list(formula = formula1,outgraphs=F,data = geno,estimator = estimate.lm.MBIC2,presearch=T, locstop =T,estimator.args =  list(data = geno),recalc_margin = 249,gen.prob = c(1,0,0,0,0), save.beta = F,interact = T,relations=c("cos","sigmoid","tanh","atan","sin","erf"),relations.prob =c(0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=4,mutation_rate = 250,last.mutation = 15000, max.tree.size = 4, Nvars.max =40,p.allow.replace=0.7,p.allow.tree=0.2,p.nor=0,p.and = 0.9),n.models = 20000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,print.freq = 100,advanced.param = list(
  max.N.glob=as.integer(10),
  min.N.glob=as.integer(5),
  max.N=as.integer(3),
  min.N=as.integer(1),
  printable = F))

res<-do.call(runemjmcmc,args = vect)

detected<-sort(mySearch$fparam[which(mySearch$p.add>0.5)])
detected<-stri_replace(str = detected,fixed = "I(",replacement = "")
detected<-stri_replace(str = detected,fixed = ")",replacement = "")

detect.true.unique<-unique(dataNeigbourhoodS2$causSNPid[which(dataNeigbourhoodS2$SNPid %in% detected)])
detect.true<-which(detected %in% dataNeigbourhoodS2$SNPid)

detlen<-length(detect.true.unique)
totlen<-length(detected)-length(detect.true)+length(detect.true.unique)

detlen/20
(totlen-detlen)/totlen



