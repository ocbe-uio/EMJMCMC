source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")


setwd("/mn/sarpanitu/ansatte-u2/aliaksah/PycharmProjects/EMJMCMC2016/examples/fromT")

data_all = read.table("Data", check.names = FALSE)
data_all = data.frame(data_all)
data_all = data_all[,c(2:dim(data_all)[2], 1)]
names(data_all)[dim(data_all)[2]] = "Y"

# Preparing data
data_all = data_all[sample(nrow(d_singles)),] # shuffle data
n = nrow(data_all)
split = floor(0.1 * n)
d_train = data_all[1:split,]
d_test = data_all[(split + 1):n,]

# Adapted from the tutorial
M = 250
# specify the initial formula
formula1 = as.formula(paste(colnames(d_train)[M+1],"~ 1 + ",
                            paste0(colnames(d_train)[-c(M+1)],collapse = "+")))

# Bayesian logic regression with the Jeffrey's prior
res4J = LogicRegr(formula = formula1, data = d_train,
                  family = "Bernoulli", prior = "J", report.level = 0.5,
                  d = 50, cmax = 2, kmax = 20, p.and = 0.9, p.not = 0.1, p.surv = 0.2,
                  ncores = 20)