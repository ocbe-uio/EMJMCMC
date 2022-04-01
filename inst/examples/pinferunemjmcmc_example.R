
# inference

X <- read.csv("inst/extdata/exa1.csv")
data.example <- as.data.frame(X)

# specify the initial formula
formula1 <- as.formula(
  paste(colnames(X)[5], "~ 1 +", paste0(colnames(X)[-5], collapse = "+"))
)

# define the number or cpus
M <- 4
# define the size of the simulated samples
NM <- 1000
# define \k_{max} + 1 from the paper
compmax <- 16
# define treshold for preinclusion of the tree into the analysis
th <- (10)^(-5)
# define a final treshold on the posterior marginal probability for reporting a
# tree
thf <- 0.05
# specify tuning parameters of the algorithm for exploring DBRM of interest
# notice that allow_offsprings=3 corresponds to the GMJMCMC runs and
# allow_offsprings=4 -to the RGMJMCMC runs
\dontrun{
  res1 <- pinferunemjmcmc(
    n.cores = M, report.level = 0.5, num.mod.best = NM, simplify = TRUE,
    runemjmcmc.params = list(
      formula = formula1, data = data.example, estimator = estimate.gamma.cpen_2,
      estimator.args = list(data = data.example), recalc_margin = 249,
      save.beta = FALSE, interact = TRUE, outgraphs = FALSE,
      relations = c("to23", "expi", "logi", "to35", "sini", "troot", "sigmoid"),
      relations.prob = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
      interact.param = list(allow_offsprings = 3, mutation_rate = 250,
      last.mutation = 10000, max.tree.size = 5, Nvars.max = 15,
      p.allow.replace = 0.9, p.allow.tree = 0.01, p.nor = 0.9, p.and = 0.9),
      n.models = 10000, unique = TRUE, max.cpu = 4, max.cpu.glob = 4,
      create.table = FALSE, create.hash = TRUE, pseudo.paral = TRUE,
      burn.in = 100, print.freq = 1000,
      advanced.param = list(
        max.N.glob = as.integer(10),
        min.N.glob = as.integer(5),
        max.N = as.integer(3),
        min.N = as.integer(1),
        printable = FALSE
      )
    )
  )
  print(res1$feat.stat)
}

# prediction

compmax <- 21

# read in the train and test data sets
test <- read.csv("inst/extdata/breast_cancer_test.csv", header = TRUE, sep = ",")[, -1]
train <- read.csv("inst/extdata/breast_cancer_train.csv", header = TRUE, sep = ",")[, -1]

# transform the train data set to a data.example data.frame that EMJMCMC class
# will internally use
data.example <- as.data.frame(train)

# specify the link function that will be used in the prediction phase
g <- function(x) {
  return((x <- 1 / (1 + exp(-x))))
}

formula1 <- as.formula(
  paste(
    colnames(data.example)[31], "~ 1 +",
    paste0(colnames(data.example)[-31], collapse = "+")
  )
)

\dontrun{
  res <- pinferunemjmcmc(
    n.cores = 30, report.level = 0.5, num.mod.best = NM, simplify = TRUE,
    predict = TRUE, test.data = as.data.frame(test), link.function = g,
    runemjmcmc.params = list(
      formula = formula1, data = data.example, gen.prob = c(1, 1, 1, 1, 0),
      estimator = estimate.bas.glm.cpen,
      estimator.args = list(
        data = data.example, prior = BAS::aic.prior(), family = binomial(),
        yid = 31, logn = log(143), r = exp(-0.5)
      ), recalc_margin = 95, save.beta = TRUE, interact = TRUE,
      relations = c("gauss", "tanh", "atan", "sin"),
      relations.prob = c(0.1, 0.1, 0.1, 0.1),
      interact.param = list(
        allow_offsprings = 4, mutation_rate = 100, last.mutation = 1000,
        max.tree.size = 6, Nvars.max = 20, p.allow.replace = 0.5,
        p.allow.tree = 0.4, p.nor = 0.3, p.and = 0.9
      ), n.models = 7000, unique = TRUE, max.cpu = 4, max.cpu.glob = 4,
      create.table = FALSE, create.hash = TRUE, pseudo.paral = TRUE,
      burn.in = 100, print.freq = 1000,
      advanced.param = list(
        max.N.glob = as.integer(10), min.N.glob = as.integer(5),
        max.N = as.integer(3), min.N = as.integer(1), printable = FALSE
      )
    )
  )

  print(auc(prob = res$predictions, y = test$X))
  for (i in 1:M) {
    print(auc(prob = res$threads.stats[[i]]$preds, y = test$X))
    print(res$threads.stats[[i]]$post.populi)
  }

  for (jjjj in 1:10)
  {
    resw <- as.integer(res$predictions >= 0.1 * jjjj)
    prec <- (1 - sum(abs(resw - test$X), na.rm = T) / length(resw))
    print(prec)
    # FNR
    ps <- which(test$X == 1)
    fnr <- sum(abs(resw[ps] - test$X[ps])) / (sum(abs(resw[ps] - test$X[ps])) + length(ps))

    # FPR
    ns <- which(test$X == 0)
    fpr <- sum(abs(resw[ns] - test$X[ns])) / (sum(abs(resw[ns] - test$X[ns])) + length(ns))
  }
}
