# #***********************IMPORTANT******************************************************
# # if a multithreaded backend openBLAS for matrix multiplications
# # is installed on your machine, please force it to use 1 thread explicitly
# # unless ncores in LogrRegr is reasonably small, in the latter case
# # you might want to experiment with the combinations of blas_set_num_threads and ncores
# library(RhpcBLASctl)
# blas_set_num_threads(1)
# omp_set_num_threads(1)
# #***********************IMPORTANT******************************************************

n_threads <- 1L
set.seed(040590)

# construct a correlation matrix for M = 50 variables
M <- 50
m <- clusterGeneration::rcorrmatrix(M, alphad = 2.5)

# simulate 1000 binary variables with this correlation matrix
sample_size <- 1000L
invisible(
  capture.output(
    X <- suppressWarnings(
      bindata::rmvbin(sample_size, margprob = rep(0.5, M), bincorr = m)
    )
  )
)

# prepare the correlation matrix in the melted format
melted_cormat <- reshape2::melt(cor(X))

# simulate Gaussian responses from a model with two-way interactions
# which is considered in S.4 of the paper
df <- data.frame(X)
df$Y <- rnorm(
  n = sample_size,
  mean = 1 + 1.43 * (df$X5 * df$X9) + 0.89 * (df$X8 * df$X11) + 0.7 * (df$X1 * df$X4),
  sd = 1
)

# specify the initial formula
formula1 <- as.formula(
  paste(
    colnames(df)[M + 1],
    "~ 1 + ",
    paste0(colnames(df)[-c(M + 1)], collapse = "+")
  )
)

# Bayesian logic regression with the robust-g-prior
res4G <- LogicRegr(
  formula = formula1, data = df, family = "Gaussian", prior = "G",
  report.level = 0.5, d = 15, cmax = 2, kmax = 15, p.and = 0.9, p.not = 0.1,
  p.surv = 0.2, ncores = n_threads, print.freq = 0L
)

# Bayesian logic regression with the Jeffreys prior
res4J <- LogicRegr(
  formula = formula1, data = df, family = "Gaussian", prior = "J",
  report.level = 0.5, d = 15, cmax = 2, kmax = 15, p.and = 0.9, p.not = 0.1,
  p.surv = 0.2, ncores = n_threads, print.freq = 0L
)

# NULLs are expected because predict = FALSE on LogicRegr

# # print the results for the robust g-prior
test_that("Results for the G-prior are sensible", {
  for (i in seq_len(nrow(res4G$allposteriors))) {
    expect_gte(res4G$allposteriors[i, 2], 0)
    expect_lte(res4G$allposteriors[i, 2], 1)
  }
  expect_gte(res4G$threads.stats[[1]]$post.populi, 0)
  expect_gt(res4G$threads.stats[[1]]$cterm, 1000)
  expect_equal(res4G$threads.stats[[1]]$preds, NULL)
  expect_length(res4G, 4L)
})

# #print the results for the Jeffreys prior
test_that("Results for the Jeffrey's prior are sensible", {
  for (i in seq_len(nrow(res4J$allposteriors))) {
    expect_gte(res4J$allposteriors[i, 2], 0)
    expect_lte(res4J$allposteriors[i, 2], 1)
  }
  expect_gte(res4J$threads.stats[[1]]$post.populi, 0)
  expect_gt(res4J$threads.stats[[1]]$cterm, -1000)
  expect_equal(res4J$threads.stats[[1]]$preds, NULL)
  expect_length(res4J, 4L)
})

# simulate Gaussian responses from a model with two-way interactions
# and an age effect which is an extension of S.4 of the paper
Xp <- data.frame(X)
Xp$age <- rpois(sample_size,lambda = 34)
Xp$Y <- rnorm(
  n = sample_size,
  mean = 1 + 0.7 * (Xp$X1 * Xp$X4) + 0.89 * (Xp$X8 * Xp$X11) + 1.43 * (Xp$X5 * Xp$X9) + 2 * Xp$age,
  sd = 1
)
teid  <- sample.int(size = 100, n = sample_size, replace = FALSE)
test  <- Xp[teid, ]
train <- Xp[-teid, ]

# specify the initial formula
formula1 = as.formula(
  paste("Y ~ 1 +", paste0(colnames(test)[-c(51, 52)], collapse = "+"))
)

# specify the link function
g = function(x) x

# specify the parameters of the custom estimator function
estimator.args <- list(
  data = train,
  n = dim(train)[1],
  m = stringi::stri_count_fixed(as.character(formula1)[3],"+"),
  k.max = 15
)

# specify the parameters of gmjmcmc algorithm
gmjmcmc.params <- list(
  allow_offsprings = 1, mutation_rate = 250, last.mutation = 10000,
  max.tree.size = 5, Nvars.max = 15, p.allow.replace = 0.9, p.allow.tree = 0.01,
  p.nor = 0.01, p.and = 0.9
)

# specify some advenced parameters of mjmcmc
mjmcmc.params <- list(
  max.N.glob = 10, min.N.glob = 5, max.N = 3, min.N = 1, printable = FALSE
)

if (interactive()) {
  # run the inference of BLR with a non-binary covariate and predicions
  set.seed(4)
  res.alt <- suppressMessages(
    pinferunemjmcmc(
      n.cores = n_threads, report.level = 0.2,
      num.mod.best = 10, simplify = FALSE,
      predict = FALSE, test.data = test, link.function = g,
      runemjmcmc.params = list(
        formula = formula1, latnames = c("I(age)"), data = train,
        estimator = estimate.logic.lm.jef, estimator.args = estimator.args,
        recalc_margin = 249, save.beta = TRUE, interact = TRUE, outgraphs = FALSE,
        interact.param = gmjmcmc.params, n.models = 10000, unique = FALSE,
        max.cpu = n_threads, max.cpu.glob = n_threads, create.table = FALSE, create.hash = TRUE,
        pseudo.paral = FALSE, burn.in = 1000, print.freq = 0,
        advanced.param = mjmcmc.params
      )
    )
  )

  test_that("Output with non-binary convariance is correct", {
    expect_equal(
      res.alt$feat.stat[, 1],
      c(
        "I(age)", "I(X5)", "I(X8)", "I(X1)", "I(((X9))&((X11)))", "I(X9)",
        "I(X4)"
      )
    )
    expect_equal(
      res.alt$feat.stat[, 2],
      c(
        "1", "0.999999994807055", "0.999956729898551", "0.999577926669365",
        "0.98582433936886", "0.98377325490383", "0.934544180132695"
      )
    )
    expect_equal(
      sqrt(mean((res.alt$predictions - test$Y) ^ 2)), 1.184527, tolerance = 1e-6
    )
    tmean <- 1 + 2 * test$age + 0.7 * (test$X1 * test$X4) + 0.89 * (test$X8 * test$X11) + 1.43 * (test$X5 * test$X9)
    expect_equal(sqrt(mean((tmean -test$Y)^2)), 1.06036, tolerance = 1e-6)
    expect_equal(mean(abs((tmean -test$Y))), .8067262, tolerance = 1e-6)
  })
}
