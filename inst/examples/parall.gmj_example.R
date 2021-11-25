j <- 1
M <- 4
X4 <- as.data.frame(array(data = rbinom(n = 50 * 1000, size = 1, prob = runif(n = 50 * 1000, 0, 1)), dim = c(1000, 50)))
Y4 <- rnorm(n = 1000, mean = 1 + 7 * (X4$V4 * X4$V17 * X4$V30 * X4$V10) + 7 * (X4$V50 * X4$V19 * X4$V13 * X4$V11) + 9 * (X4$V37 * X4$V20 * X4$V12) + 7 * (X4$V1 * X4$V27 * X4$V3) + 3.5 * (X4$V9 * X4$V2) + 6.6 * (X4$V21 * X4$V18) + 1.5 * X4$V7 + 1.5 * X4$V8, sd = 1)
X4$Y4 <- Y4

formula1 <- as.formula(paste(colnames(X4)[51], "~ 1 +", paste0(colnames(X4)[-c(51)], collapse = "+")))
data.example <- as.data.frame(X4)

vect <- list(formula = formula1, outgraphs = FALSE, data = X4, estimator = EMJMCMC:::estimate.logic.lm, estimator.args = list(data = data.example, n = 100, m = 50), recalc_margin = 249, save.beta = FALSE, interact = TRUE, relations = c("", "lgx2", "cos", "sigmoid", "tanh", "atan", "erf"), relations.prob = c(0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), interact.param = list(allow_offsprings = 1, mutation_rate = 250, last.mutation = 15000, max.tree.size = 4, Nvars.max = 40, p.allow.replace = 0.7, p.allow.tree = 0.2, p.nor = 0, p.and = 0.9), n.models = 20000, unique = TRUE, max.cpu = 4, max.cpu.glob = 4, create.table = FALSE, create.hash = TRUE, pseudo.paral = TRUE, burn.in = 50, print.freq = 1000, advanced.param = list(
  max.N.glob = as.integer(10),
  min.N.glob = as.integer(5),
  max.N = as.integer(3),
  min.N = as.integer(1),
  printable = FALSE
))

params <- list(vect)[rep(1, M)]

for (i in 1:M)
{
  params[[i]]$cpu <- i
  params[[i]]$NM <- 1000
  params[[i]]$simlen <- 21
}
gc()
print(paste0("begin simulation ", j))
results <- EMJMCMC:::parall.gmj(X = params, M = 1) # increase M if you are not on Windows!
