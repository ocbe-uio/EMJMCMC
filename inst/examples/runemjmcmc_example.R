X4 <- as.data.frame(
  array(
    data = rbinom(n = 50 * 1000, size = 1, prob = runif(n = 50 * 1000, 0, 1)),
    dim = c(1000, 50)
  )
)
Y4 <- rnorm(
  n = 1000,
  mean = 1 +
    7 * (X4$V4 * X4$V17 * X4$V30 * X4$V10) +
    7 * (((X4$V50 * X4$V19 * X4$V13 * X4$V11) > 0)) +
    9 * (X4$V37 * X4$V20 * X4$V12) +
    7 * (X4$V1 * X4$V27 * X4$V3) +
    3.5 * (X4$V9 * X4$V2) +
    6.6 * (X4$V21 * X4$V18) +
    1.5 * X4$V7 +
    1.5 * X4$V8,
  sd = 1
)
X4$Y4 <- Y4
data.example <- as.data.frame(X4)

# specify the initial formula
formula1 <- as.formula(
  paste(colnames(X4)[51], "~ 1 +", paste0(colnames(X4)[-c(51)], collapse = "+"))
)


# specify tuning parameters of the algorithm for exploring DBRM of interest
# notice that allow_offsprings=3 corresponds to the GMJMCMC runs and
# allow_offsprings=4 -to the RGMJMCMC runs
\donttest{
  res <- runemjmcmc(
    formula = formula1, outgraphs = FALSE, data = X4,
    estimator = estimate.gamma.cpen, estimator.args = list(data = data.example),
    recalc_margin = 249, save.beta = FALSE, interact = TRUE,
    relations = c("cos", "sigmoid", "tanh", "atan", "sin", "erf"),
    relations.prob = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
    interact.param = list(
      allow_offsprings = 4, mutation_rate = 250, last.mutation = 15000,
      max.tree.size = 4, Nvars.max = 40, p.allow.replace = 0.7,
      p.allow.tree = 0.2, p.nor = 0, p.and = 0.9
    ), n.models = 20000, unique = TRUE, max.cpu = 4, max.cpu.glob = 4,
    create.table = FALSE, create.hash = TRUE, pseudo.paral = TRUE, burn.in = 50,
    print.freq = 1000,
    advanced.param = list(
      max.N.glob = as.integer(10),
      min.N.glob = as.integer(5),
      max.N = as.integer(3),
      min.N = as.integer(1),
      printable = FALSE
    )
  )
}
