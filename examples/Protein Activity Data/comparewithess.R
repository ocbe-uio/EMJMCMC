

formula.prot.1 = as.formula(paste0(fobserved.example,"~",paste(fparam.example[c(11,38,46,68)],collapse = "+")))

formula.prot.2 = as.formula(paste0(fobserved.example,"~",paste(fparam.example[c(5,6,10,14,16,20,27,38,46,54,69,72)],collapse = "+")))

formula.prot.3 = as.formula(paste0(fobserved.example,"~",paste(fparam.example[c(3,6,7,8,10,11,16,26,40,58,64,69,76,86)],collapse = "+")))

formula.prot.4 = as.formula(paste0(fobserved.example,"~",paste(fparam.example[c(6,10,11,38,46,69,86)],collapse = "+")))


estimate.bas.lm(formula = formula.prot.1,data = data.example,prior = 3, n = 96, g = 96)$mlik + 149.3995468795

estimate.bas.lm(formula = formula.prot.2,data = data.example,prior = 3, n = 96, g = 96)$mlik + 150.5537667327

estimate.bas.lm(formula = formula.prot.3,data = data.example,prior = 3, n = 96, g = 96)$mlik + 151.3277283846

estimate.bas.lm(formula = formula.prot.4,data = data.example,prior = 3, n = 96, g = 96)$mlik + 136.9203682823



marg.inclusions.ess <- read.table("Z:/ESS/Prot/Output/prot_Example_110000_sweeps_output_marg_prob_incl.txt",header = T)[,2]

marg.inclusions.mjm <-  read.table("Z:/ESS/Prot/Output/pp3.rs.csv",header = F)[,1]

plot(marg.inclusions.ess,marg.inclusions.mjm)


models.ess<- read.table("Z:/ESS/Prot/Output/prot_Example_110000_sweeps_output_models_history.txt",header = T,sep = ",")[,1]


models.ess<-stri_split_fixed(str = models.ess,pattern = "\t")  

mliks.ess<-array(0,length(models.ess))
mliks.mjm<-array(0,length(models.ess))

for(i in 1:length(models.ess))
{
  mliks.ess[i]<-as.numeric(models.ess[i][[1]][3])
  covs<-stri_split_fixed(str = models.ess[i][[1]][length(models.ess[i][[1]])],pattern = " ")[[1]]
  covs<-covs[-length(covs)]
  formula.prot = as.formula(paste0(fobserved.example,"~",paste(fparam.example[as.numeric(covs)],collapse = "+")))
  
  mliks.mjm[i]<-estimate.bas.lm(formula = formula.prot,data = data.example,prior = 3, n = 96, g = 96)$mlik
}
  

plot(mliks.ess,mliks.mjm)

mean(mliks.ess-mliks.mjm)
var(mliks.ess-mliks.mjm)


max(mliks.mjm)