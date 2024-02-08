# simulate Gaussian responses
threads <- 1L
n_dims <- c(11, 100)
set.seed(040590)
X1 <- as.data.frame(
  array(
    data = rbinom(n = prod(n_dims), size = 1, prob = runif(n = prod(n_dims), 0, 1)),
    dim = rev(n_dims)
  )
)
Y1 <- rnorm(
  n = n_dims[2],
  mean = 1 + 0.7 * (X1$V1 * X1$V4) + 0.8896846 * (X1$V8 * X1$V11) + 1.434573 * (X1$V5 * X1$V9),
  sd = 1
)
X1$Y1 <- Y1

# specify the initial formula
formula1 <- as.formula(
  paste(colnames(X1)[ncol(X1)], "~ 1 +", paste0(colnames(X1)[-c(ncol(X1))], collapse = "+"))
)
data.example <- as.data.frame(X1)

# run the inference with robust g prior
res4G <- suppressMessages(
  LogicRegr(
    formula = formula1, data = data.example, family = "Gaussian", prior = "G",
    report.level = 0.5, d = 15, cmax = 2, kmax = 15, p.and = 0.9, p.not = 0.01,
    p.surv = 0.2, ncores = threads, print.freq = 0L
  )
)
# run the inference with Jeffrey's prior
res4J <- suppressMessages(
  LogicRegr(
    formula = formula1, data = data.example, family = "Gaussian", prior = "J",
    report.level = 0.5, d = 15, cmax = 2, kmax = 15, p.and = 0.9, p.not = 0.01,
    p.surv = 0.2, ncores = threads, print.freq = 0L
  )
)

# change to Bernoulli responses
X1 <- as.data.frame(
  array(data = rbinom(n = prod(n_dims), size = 1, prob = 0.3), dim = rev(n_dims))
)
Y1 <- -0.7 + 1 * ((1 - X1$V1) * (X1$V4)) + 1 * (X1$V8 * X1$V11) + 1 * (X1$V5 * X1$V9)
X1$Y1 <- round(1.0 / (1.0 + exp(-Y1)))

# specify the initial formula
formula1 <- as.formula(
  paste(colnames(X1)[ncol(X1)], "~ 1 +", paste0(colnames(X1)[-c(ncol(X1))], collapse = "+"))
)
data.example <- as.data.frame(X1)

# run the inference with robust g prior
res1G <- suppressWarnings(
  suppressMessages(
    LogicRegr(
      formula = formula1, data = data.example, family = "Bernoulli", prior = "G",
      report.level = 0.5, d = 15, cmax = 2, kmax = 15, p.and = 0.9, p.not = 0.2,
      p.surv = 0.2, ncores = threads, print.freq = 0L
    )
  )
)

# run the inference with Jeffrey's prior
res1J <- suppressWarnings(
  suppressMessages(
    LogicRegr(
      formula = formula1, data = data.example, family = "Bernoulli", prior = "J",
      report.level = 0.5, d = 15, cmax = 2, kmax = 15, p.and = 0.9, p.not = 0.2,
      p.surv = 0.2, ncores = threads, print.freq = 0L
    )
  )
)
test_that("outputs are correct", {
  expect_equal(
    res4G$feat.stat,
    matrix(c("I(((V5))&((V9)))", "0.998464672344289"), nrow = 1)
  )
  expect_equal(
    res4J$feat.stat,
    matrix(c( "I(((V5))&((V9)))", "0.995269065749787"), nrow = 1)
  )
  expect_equal(
    res1G$feat.stat,
    matrix(
      c(
        "I((1-(V1))&((V4)))", "0.999999996715506",
        "I(V8)", "0.999999302544405",
        "I(V11)", "0.999997483629185",
        "I(((((V5))))&((V9)))", "0.958531256229637"
      ),
      nrow = 4, ncol = 2, byrow = TRUE
    )
  )
  expect_equal(
    res1J$feat.stat,
    matrix(
      c(
        "I(((V11))&((V8)))", "0.999994827992723",
        "I(V1)", "0.97262338705771",
        "I(V5)", "0.949549067791939",
        "I(((V9)))", "0.941316258317002"
      ),
      nrow = 4, ncol = 2, byrow = TRUE
    )
  )
})
