

source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")

x11(display = "localhost:11.0")

setwd("/mn/sarpanitu/ansatte-u2/aliaksah/PycharmProjects/EMJMCMC2016/examples/fromT")

data_all = read.table("Data_50mostFrequentCovariates", check.names = FALSE)
data_all = data.frame(data_all)
data_all = data_all[,c(2:dim(data_all)[2], 1)]
names(data_all)[dim(data_all)[2]] = "Y"

# Preparing data
#data_all = data_all[sample(nrow(d_singles)),] # shuffle data
n = nrow(data_all)
split = floor(0.2 * n)

#cors = cor(data_all)[,251]

#which(abs(cors[1:250])>0.05)
trid= sample.int(n,0.2*n, replace = F)
d_train = (data_all[trid,c(51,1:50)])
test = (data_all[-trid,c(51,1:50)])

cors = cor(d_train)[,1]

d_test = data_all[(split + 1):n,]

row.names(d_train) = 1:dim(d_train)[1]

# Adapted from the tutorial
#M = 250
# specify the initial formula
formula1 = as.formula(paste(colnames(d_train)[1],"~ 1 + ",
                            paste0(colnames(d_train)[-c(1)],collapse = "+")))


data.example = d_train
# Bayesian logic regression with the Jeffrey's prior
res4J = LogicRegr(formula = formula1, data = data.example ,
                  family = "Bernoulli", prior = "G", report.level = 0.1,
                  d = 20, cmax = 3, kmax = 15, p.and = 0.9, p.not = 0.1, p.surv = 0.2,
                  ncores = 2)


print(res4J$feat.stat)

res4J$allposteriors

g=function(x)
{
  return((x = 1/(1+exp(-x))))
}
# specify the parameters of the custom estimator function
estimator.args = list(data = data.example, n = dim(data.example)[1],
                      m =stri_count_fixed(as.character(formula1)[3],"+"),k.max = 15)
# specify the parameters of gmjmcmc algorithm
gmjmcmc.params = list(allow_offsprings=1,mutation_rate = 250,
                      last.mutation=10000, max.tree.size = 3, Nvars.max =15,
                      p.allow.replace=0.9,p.allow.tree=0.01,p.nor=0.1,p.and = 0.9)
# specify some advenced parameters of mjmcmc
mjmcmc.params = list(max.N.glob=10, min.N.glob=5, max.N=3, min.N=1,
                     printable = F)
# run the inference of BLR with a non-binary covariate and predicions
res.alt = pinferunemjmcmc(n.cores = 4, report.level =  0.2,
                          num.mod.best = 500,simplify = T,predict = T,test.data = test,
                          link.function = g,
                          runemjmcmc.params = list(formula = formula1,
                                                   data = data.example,estimator = estimate.logic.bern.tCCH,
                                                   estimator.args =estimator.args,
                                                   recalc_margin = 249, save.beta = T,interact = T,outgraphs=F,
                                                   interact.param = gmjmcmc.params,
                                                   n.models = 20000,unique = F,max.cpu = 4,max.cpu.glob = 4,
                                                   create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,
                                                   print.freq = 1000,
                                                   advanced.param = mjmcmc.params))



print(base::rbind(c("expressions","probabilities"),res.alt$feat.stat))


print(1-mean(abs((as.integer(res.alt$predictions>0.875)-test$Y))))

res.alt$feat.stat