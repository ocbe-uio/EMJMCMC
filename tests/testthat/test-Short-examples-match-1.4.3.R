# Meta-objects =========================================== #

set.seed(7070865)
tol <- 1e-3

# Data =================================================== #

j <- 1
M <- 4
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
formula1 <- as.formula(
  paste(colnames(X4)[51], "~ 1 +", paste0(colnames(X4)[-c(51)], collapse = "+"))
)
formula2 <- as.formula(
  paste(
    colnames(data.example)[1], "~ 1 +", paste0(colnames(data.example)[-1],
    collapse = "+")
  )
)
X1 <- as.data.frame(
  array(data = rbinom(n = 50 * 1000, size = 1, prob = 0.3), dim = c(1000, 50))
)
Y1 <- -0.7 + 1 * ((1 - X1$V1) * (X1$V4)) + 1 * (X1$V8 * X1$V11) + 1 * (X1$V5 * X1$V9)
X1$Y1 <- round(1.0 / (1.0 + exp(-Y1)))
formula3 <- as.formula(
  paste(colnames(X1)[51], "~ 1 +", paste0(colnames(X1)[-c(51)], collapse = "+"))
)

# Individual function tests ============================== #

test_that("estimate.bas.blm output matches version 1.4.3", {
  data.example$Y4 <- as.integer(data.example$Y > mean(data.example$Y))
  out <- estimate.bas.glm(
    formula = formula1,
    data = data.example,
    prior = BAS::aic.prior(),
    logn = 47,
    family = binomial()
  )
  expect_equal(out$mlik, -51)
  expect_equal(out$waic, -102)
  expect_equal(out$dic, -2397)
  expect_equal(mean(out$summary.fixed$mean), 0.5209033, tolerance = tol)
})

test_that("estimate.bas.lm output matches version 1.4.3", {
  out <- estimate.bas.lm(
    formula = formula1,
    data = data.example,
    prior = 2,
    n = 47
  )
  expect_equal(out$mlik, -103.0293, tolerance = tol)
  expect_equal(out$waic, 5735.152, tolerance = tol)
  expect_equal(out$dic, 5990.356, tolerance = tol)
  expect_equal(mean(out$summary.fixed$mean), 0.4194573, tolerance = tol)
})

test_that("estimate.bigm output matches version 1.4.3", {
  out <- estimate.bigm(
    formula = formula2, data = data.example, n = 47, prior = "BIC", maxit = 20,
    chunksize = 1000000, family = gaussian()
  )
  expect_equal(out$mlik, -2629.5, tolerance = tol)
  expect_equal(out$waic, 334.4996, tolerance = tol)
  expect_equal(out$dic, 2629.5, tolerance = tol)
  expect_equal(mean(out$summary.fixed$mean), 0.005676248, tolerance = tol)
})

test_that("estimate.glm output matches version 1.4.3", {
  out <- estimate.glm(
    formula = formula2, data = data.example, prior = 2, family = gaussian()
  )
  expect_equal(out$mlik, -1738.214, tolerance = tol)
  expect_equal(out$waic, 1483.01, tolerance = tol)
  expect_equal(out$dic, 1738.214, tolerance = tol)
  expect_equal(mean(out$summary.fixed$mean), 0.005676248, tolerance = tol)
})

test_that("estimate.logic.glm output matches version 1.4.3", {
  out <- estimate.logic.glm(
    formula = formula3, data = X1, family = binomial(), n = 1000, m = 50
  )
  expect_equal(out$mlik, 384.3592, tolerance = tol)
  expect_equal(out$waic, -581.9598, tolerance = tol)
  expect_equal(out$dic, -832.2553, tolerance = tol)
  expect_equal(mean(out$summary.fixed$mean), 0.1606305, tolerance = tol)
})

test_that("estimate.logic.lm output matches version 1.4.3", {
  out <- estimate.logic.lm(formula = formula1, data = X4, n = 1000, m = 50)
  expect_equal(out$mlik, -2694.691, tolerance = tol)
  expect_equal(out$waic, 4735.152, tolerance = tol)
  expect_equal(out$dic, 4990.356, tolerance = tol)
  expect_equal(mean(out$summary.fixed$mean), 0.4194573, tolerance = tol)
})

test_that("estimate.speedglm output matches version 1.4.3", {
  out <- estimate.speedglm(
    formula = formula1, data = data.example, prior = "BIC", logn = log(47),
    family = gaussian()
  )
  expect_equal(out$mlik, -196.3575, tolerance = tol)
  expect_equal(out$waic, -102, tolerance = tol)
  expect_equal(out$dic, -196.3575, tolerance = tol)
  expect_equal(mean(out$summary.fixed$mean), 0.01960784, tolerance = tol)
})

test_that("simplify.formula output matches version 1.4.3", {
  out <- simplify.formula(fmla = formula1, names = colnames(X1))
  expect_equal(
    out$fparam,
    c(
      "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11",
      "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21",
      "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30", "V31",
      "V32", "V33", "V34", "V35", "V36", "V37", "V38", "V39", "V40", "V41",
      "V42", "V43", "V44", "V45", "V46", "V47", "V48", "V49", "V50"
    )
  )
  expect_equal(out$fobserved, "Y4")
})
