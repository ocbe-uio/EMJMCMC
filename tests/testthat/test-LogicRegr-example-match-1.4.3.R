set.seed(265508)
n_cores <- 1L
n_row <- 100L
n_col <- 11L
n_tot <- n_row * n_col
X1 <- as.data.frame(
  array(
    data = rbinom(n = n_tot, size = 1, prob = runif(n = n_tot)),
    dim = c(n_row, n_col)
  )
)
Y1 <- rnorm(
  n = n_row,
  mean = 1 +
    0.7 * (X1$V1 * X1$V4) +
    0.8896846 * (X1$V8 * X1$V11) +
    1.434573 * (X1$V5 * X1$V9),
  sd = 1
)
X1$Y1 <- Y1

# specify the initial formula
formula1 <- as.formula(
  paste(
    colnames(X1)[n_col + 1L], "~ 1 +",
    paste0(colnames(X1)[-c(n_col + 1L)], collapse = "+")
  )
)
data.example <- as.data.frame(X1)

# run the inference with robust g prior
res4G <- EMJMCMC::LogicRegr(
  formula = formula1, data = data.example, family = "Gaussian", prior = "G",
  report.level = 0.5, d = 15, cmax = 2, kmax = 15, p.and = 0.9, p.not = 0.01,
  p.surv = 0.2, ncores = n_cores, print.freq = 0L
)

# run the inference with Jeffrey's prior
res4J <- EMJMCMC::LogicRegr(
  formula = formula1, data = data.example, family = "Gaussian", prior = "J",
  report.level = 0.5, d = 15, cmax = 2, kmax = 15, p.and = 0.9, p.not = 0.01,
  p.surv = 0.2, ncores = n_cores, print.freq = 0L
)

test_that("LogicRegr output matches version 1.4.3", {
  obs_4G <- as.numeric(res4G$feat.stat[, 2])
  obs_4J <- as.numeric(res4J$feat.stat[, 2])
  expect_equal(ncol(res4G$feat.stat), 2L)
  expect_equal(ncol(res4G$feat.stat), 2L)
  expect_true(all(obs_4G >= 0) && all(obs_4G <= 1))
  expect_true(all(obs_4J >= 0) && all(obs_4J <= 1))
})
