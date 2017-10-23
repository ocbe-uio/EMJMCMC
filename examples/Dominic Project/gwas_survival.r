
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")
library("survival")

load("/storage/edelmand/IMPACT2/Daten/betas_clinic_for_analysis.RData")

clinical<-cbind(clinic2[,-1], t(betas_ges))


if(FALSE){
  #function to be used in screening
  coxscreen<-function(x) tryCatch(return(summary(coxph(Surv(V1,V2) ~1+V3,data = as.data.frame(cbind(clinical$fudays_35fu,clinical$death_event_35fu
  ,clinical[[x]]))))$logtest[3]),error = function(err) {print(err)
    return(1)})
  
  screen<-mclapply(FUN = coxscreen, X = names(clinical)[-c(14:18)],mc.cores = 63)
}else
{
 screen<-read.csv(file = "pvalues.csv",header = F,row.names = F,col.names = F) 
}



ids<-NULL
for(i in 1:length(screen))
  if(screen[[i]]<=0.000001)
    ids<-c(ids,i)

ids=read.csv("indeces.csv")

clinical<-clinical[,c(ids,14:18)]

options("expressions"=500000)

fla= paste0("Surv(fudays_35fu,crc_death_35fu)~1 + ",paste(names[-c(16680:16684)],collapse = "+"))

formula=as.formula(paste0("Surv(fudays_35fu,crc_death_35fu)~1 + ",paste(names[-c(16680:16684)],collapse = "+")))


#estimator function (for speed inla is not adressed here currently)

#estimator function (for speed inla is not adressed here currently)
estimate.lm.BIC <- function(formula, data, n = 1125, m = 1000, c = 16,u=150)
{
  size<-stri_count_fixed(str = as.character(formula)[3],pattern = "+")
  
  if(size>u||size ==0)
  {
    return(list(mlik = (-50000 + rnorm(1,0,1)),waic = 50000 + rnorm(1,0,1), dic =  50000+ rnorm(1,0,1),summary.fixed =list(mean = array(0,size+1))))
  }
  out<-NULL
  
  capture.output(tryCatch({out <- (coxph(formula, data=data))}))
  
  if(is.null(out))
  {
    return(list(mlik = (-50000 + rnorm(1,0,1)),waic = 50000 + rnorm(1,0,1), dic =  50000+ rnorm(1,0,1),summary.fixed =list(mean = array(0,size+1))))
  }
  
  
  logmarglik <- (2*out$loglik[1]- size*log(n))/2# for mbic would be out$loglik[1] -  size*log(m*m*n/c) ) + 2*log(factorial(size)))/2#
  return(list(mlik = logmarglik,waic = -logmarglik , dic = -logmarglik, summary.fixed = c(0,out$coefficients)))
}


data.example = clinical
rm(clinical)
gc()


vect<-list(formula = formula, outgraphs=F,data = data.example,max.time = 150, estimator = estimate.lm.BIC,presearch=F, locstop =F,estimator.args =  list(data = data.example),recalc_margin = 100,gen.prob = c(1,0,0,0,0), save.beta = F,interact = F,relations=c("cos"),relations.prob =c(0.1),interact.param=list(allow_offsprings=3,mutation_rate = 101,last.mutation = 150000, max.tree.size = 4, Nvars.max =100,p.allow.replace=0.7,p.allow.tree=0.0005,p.nor=0,p.and = 0.9),n.models = 200000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,print.freq = 100,advanced.param = list(
  max.N.glob=as.integer(100),
  min.N.glob=as.integer(10),
  max.N=as.integer(3),
  min.N=as.integer(1),
  printable = F))


res<-do.call(runemjmcmc,args =vect)


res$p.post[which(res$p.post>0.05)]
mySearch$fparam[which(res$p.post>0.05)]

