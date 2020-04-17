#load the libraries
library(RCurl)
library(latex2exp)
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")


#define the working directory
workdir=""

#get and prepare the data
M=5
size=1

data.all = read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/Epigenetic%20Data/epigen.txt"),sep = ",",header = T)[,2:30]
data.all$strfact=0
data.all$strfact[which(data.all$strand=="+")]=1
names(data.all)[c(21,23,24)] = c("Ma","Mg","Md")
data.all$maexpr=data.all$Ma*data.all$express_noisy
data.all$mgexpr=data.all$Mg*data.all$express_noisy
data.all$mdexpr=data.all$Md*data.all$express_noisy
data.all$ppos = data.all$pos 

workdir = ""

#indentification data
data.indentify=data.all[which(data.all$total_bases<=2),]


data.train=data.all[which(data.all$total_bases>2),]


#inference data
plot(data.train$total_bases, col = 4, pch = 16,ylim = c(0,30),main = "Epigenetic data sample",xlab = "nucleobase id",ylab = "reads")
points(data.train$methylated_bases,col=2, pch = 17)

#alldata
plot(data.all$total_bases, col = 4, pch = 16,ylim = c(0,30),main = "Epigenetic data sample",xlab = "nucleobase id",ylab = "reads")
points(data.all$methylated_bases,col=2, pch = 17)


#specify covariates
covariates= colnames(data.all)[c(8:10,12:17,21,23,26,29,30,31,32,33)]
#specify observations
observations = colnames(data.all)[5]

formula1 = as.formula(paste(observations ,"~ 1 +",paste0(covariates,collapse = "+")))


#specify the function for computing the marginal likelihoods with inla
estimate.inla.ar1 = function(formula, data, family = "binomial",Ntrials = data.train$total_bases,control.compute = list(waic = T,dic = T))
{
  
  out=NULL
  #print(formula)
  capture.output({withRestarts(tryCatch(capture.output({out = inla(formula = formula,Ntrials =Ntrials,control.compute = control.compute,family = family,data = data) })), abort = function(){onerr=TRUE;out=NULL})})
 
  #on error return extremely small values of the marginal lihelihoods
  if(is.null(out))
  {
    return(list(mlik = -10000,waic =10000, dic =10000, summary.fixed = 0))
  }
  #otherwise return the values for a given model
  coef=out$summary.fixed$mode
  print(out$mlik[1])
  return(list(mlik = out$mlik[1],waic = out$waic[1] , dic = out$dic[1], summary.fixed =list(mean =coef)))
  
}




toplot = data.all[covariates]
names(toplot)= c(("DIST"),("CHG"),("CGH"),("DT1"),("DT2"),("DT3"),("DT4"),("DT5"),("DT6:20"),("Ma"),("Mg"),("CODE"),("EXPR"),("STRD"),("EXPRa"),("EXPRg"),("EXPRd"))
# prepare the correlation matrix in the melted format
melted_cormat = reshape2::melt(cor(toplot))
# plot the heat-map of the correlations
pdf(file=paste(workdir,"","ex5corr.pdf",sep = ""),width  =14,
    height    = 10)
ggplot2::ggplot(data = melted_cormat, 
                ggplot2::aes(x=Var1, y=Var2, fill=value)) + 
  ggplot2::geom_tile() +  
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.title.y =  ggplot2::element_blank())
                # axis.text.y = ggplot2::element_text(c(("DIST"),("CHG"),("CGH"),("DT1"),("DT2"),("DT3"),("DT4"),("DT5"),("DT6:20"),("Ma"),("Mg"),("CODE"),("EXPR"),("STRD"),("EXPRa"),("EXPRg"),("EXPRd"))),
                 #axis.text.x = ggplot2::element_blank())
dev.off()


#run the mjmcmc inference
data.example = as.data.frame(data.train)
res = runemjmcmc(formula = formula1,data = data.train,recalc_margin = 2^10,latent = "+ f(data.example$ppos,model=\"iid\") + f(data.example$pos,model=\"rw1\")", estimator =estimate.inla.ar1,estimator.args =  list(family = "binomial",data = data.train, Ntrials = data.train$total_bases,control.compute = list(waic = T,dic = T)),save.beta = T,interact = F,relations = c("","sin","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=2,mutation_rate = 100, max.tree.size = 200000, Nvars.max = 95,p.allow.replace=0.9,p.allow.tree=0.5,p.nor=0.3,p.and = 0.7),n.models = 10000,unique = T,max.cpu = 20,max.cpu.glob = 10,create.table = T,create.hash = F,pseudo.paral = F,burn.in = 100,print.freq = 100)


pdf(file=paste(workdir,"","barrm.pdf",sep = ""),width  = 12,
    height    = 5)
barplot(res$p.post,names.arg = c(latex2exp("X_{DIST}"),latex2exp("X_{CHG}"),latex2exp("X_{CGH}"),latex2exp("X_{DT1}"),latex2exp("X_{DT2}"),latex2exp("X_{DT3}"),latex2exp("X_{DT4}"),latex2exp("X_{DT5}"),latex2exp("X_{DT6:20}"),latex2exp("X_{M_{a}}"),latex2exp("X_{M_{g}}"),latex2exp("X_{CODE}"),latex2exp("X_{EXPR}"),latex2exp("X_{STRD}"),latex2exp("X_{EXPR,a}"),latex2exp("X_{EXPR,g}"),latex2exp("X_{EXPR,d}")),ylab = "Marginal inclusion probabilities",las=1,density = 46,border="black")
dev.off()


#analyze the results from the search
ppp=mySearch$post_proceed_results(statistics1 = statistics1)
truth = ppp$p.post # make sure it is equal to Truth column from the article
truth.m = ppp$m.post
truth.prob = ppp$s.mass
ordering = sort(ppp$p.post,index.return=T)
print("pi truth")
sprintf("%.10f",truth[ordering$ix])
sprintf(covariates[ordering$ix])



# n/b visualization iisues on Windows! To be fixed!
template = "INLA"
mySearch$visualize_results(statistics1, template, 200, crit=list(mlik = T, waic = T, dic = T),draw_dist =F)
View(statistics1[,])

#look at the best model
best=covariates[which(!is.na(statistics1[which(statistics1[,1]==max(statistics1[,1],na.rm = T))[1],17:(16+length(covariates))]))]
#logit link function
g=function(x)
{
  return((x = 1/(1+exp(-x))))
}
# condiotion on posterior mode in  the model space
data.train3=data.all
idtest=which(data.all$total_bases<=2)
idtrain=which(data.all$total_bases>2)
data.train3$methylated_bases[idtest]=NA
data.train3$ppos = data.train3$pos
#example of the underlying model within INLA
formula2 = as.formula(paste0("methylated_bases ~ 1+",paste(best,collapse = "+"), " + f(data.train3$ppos,model=\"iid\")", " + f(data.train3$pos,model=\"rw1\")"))
args=list(family = "binomial",data = data.train3, Ntrials = data.train3$total_bases, control.compute = list(waic = T,dic = T))
fm5=do.call(inla, c(args,formula = formula2))
summary(fm5)
fm5$summary.hyperpar$mean[1]
coef=fm5$summary.fixed$mode
coef[1]=coef[1]+fm5$summary.hyperpar$mean[1]

#plot the latent state across the genome
png(file=paste(workdir,"","latent.png",sep = ""),width     = 10,
    height    = 5,
    units     = "in",
    res       = 500)
plot(fm5$summary.random$`data.train3|S|pos`$mean+fm5$summary.random$`data.train3|S|ppos`$mean,  ylab = "RW(1) + IG", xlab = "",xaxt = "n")
axis(1, at=seq(1,length(data.all$pos), by = 50), labels=data.all$pos[seq(1,length(data.all$pos), by = 50)],las=2)
dev.off()

#now check glm results (no random latent GPs) corresoding to the covariates from the best found glmm model 
formula3 = as.formula(paste0("cbind(methylated_bases,unmethylated_bases)~ 1 +",paste(best,collapse = "+")))
fit1 = glm( formula = formula3,family=binomial,data=data.train)
summary(fit1)
#here all of the found covariates are also significant

#now plot the overall results
pdf(file=paste(workdir,"classify.pdf",sep = ""),width     = 10,
    height    = 5)
plot(data.all$total_bases, col = 4, pch = 16,ylim = c(-13,30),xlab = "",ylab = "Data and methylation probabilities",xaxt = "n",yaxt = "n")
axis(1, at=seq(1,length(data.all$pos), by = 50), labels=data.all$pos[seq(1,length(data.all$pos), by = 50)],las=2)
axis(2, at=c(seq(0,max(data.all$total_bases)+5,5)), labels = c(seq(0,max(data.all$total_bases)+5,5)))
axis(2, at=c(-13,-3), labels = c(0,1), las = 2)
points(data.all$methylated_bases,col=2, pch = 15)
lines(rep(2.6,length(data.all$total_bases)),col=3,lwd=4,lty = 1)
lines(rep(- 2.8,length(data.all$total_bases)))
lines(10*data.all$methylated_bases/data.all$total_bases-13,col = 5,lwd=2)
lines(10*g(fm5$summary.random[[2]]$mean)-13,lwd=4,col = "sienna4")
dev.off()


#Appedix

#compute the values of marginal likelihood for the alternative latent models (one might to restart R prior to that)
mliks.full = NULL
mliks.best = NULL
mliks.null = NULL
for(proc in c("ou","rw1","ar1","iid"))
{
  
  mliks.full[[proc]] = estimate.inla.ar1(formula = methylated_bases ~ 1 + I(base_dist) + I(CHG) + I(CG) + I(DT1) + 
                                           I(DT2) + I(DT3) + I(DT4) + I(DT5) + I(DT6_20) + I(Ma) + I(Mg) + 
                                           I(coding) + I(express_noisy) + I(strfact) + I(maexpr) + I(mgexpr) + 
                                           I(mdexpr)+ f(data.train$ppos,model="iid")+ f(data.train$pos,model=proc),data = data.train)$mlik
  
  mliks.null[[proc]] = estimate.inla.ar1(formula = methylated_bases ~ 1 + f(data.train$ppos,model="iid")+ f(data.train$pos,model=proc),data = data.train)$mlik 
  mliks.best[[proc]] = estimate.inla.ar1(formula = methylated_bases ~ 1 + I(CHG) + I(CG) + I(coding) + f(data.train$ppos,model="iid")+ f(data.train$pos,model=proc),data = data.train)$mlik 
}

#compute the values of marginal likelihood for the alternative latent models
for(proc in c(2,3,4))
{
  
  mliks.full[[paste0("ar",proc)]] = estimate.inla.ar1(formula = methylated_bases ~ 1 + I(base_dist) + I(CHG) + I(CG) + I(DT1) + 
                                                        I(DT2) + I(DT3) + I(DT4) + I(DT5) + I(DT6_20) + I(Ma) + I(Mg) + 
                                                        I(coding) + I(express_noisy) + I(strfact) + I(maexpr) + I(mgexpr) + 
                                                        I(mdexpr)+ f(data.train$ppos,model="iid")+ f(data.train$pos,model="ar",order = proc),data = data.train)$mlik
  
  gc()
  mliks.null[[paste0("ar",proc)]] = estimate.inla.ar1(formula = methylated_bases ~ 1 + f(data.train$ppos,model="iid")+ f(data.train$pos,model="ar",order = proc),data = data.train)$mlik 
  
  gc()
  mliks.best[[paste0("ar",proc)]] = estimate.inla.ar1(formula = methylated_bases ~ 1 + I(CHG) + I(CG) + I(coding) + f(data.train$ppos,model="iid")+ f(data.train$pos,model="ar",order = proc),data = data.train)$mlik 
  
  gc()
}
