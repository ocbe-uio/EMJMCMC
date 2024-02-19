set.seed(80334)
n_cores <- min(parallel::detectCores() - 1L, 20L)
M <- 1L
n_row <- 100L
n_col <- 50L
n_tot <- n_row * n_col
X4 <- as.data.frame(
  array(
    data = rbinom(n = n_tot, size = 1, prob = runif(n = n_tot)),
    dim = c(n_row, n_col)
  )
)
Y4 <- rnorm(
  n = n_row,
  mean = 1 +
    7 * (X4$V4 * X4$V17 * X4$V30 * X4$V10) +
    7 * (X4$V50 * X4$V19 * X4$V13 * X4$V11) +
    9 * (X4$V37 * X4$V20 * X4$V12) +
    7 * (X4$V1 * X4$V27 * X4$V3) +
    3.5 * (X4$V9 * X4$V2) +
    6.6 * (X4$V21 * X4$V18) +
    1.5 * X4$V7 +
    1.5 * X4$V8,
  sd = 1
)
X4$Y4 <- Y4

formula1 <- as.formula(
  paste(
    colnames(X4)[n_col + 1L], "~ 1 +",
    paste0(colnames(X4)[-c(n_col + 1L)], collapse = "+")
  )
)
data.example <- as.data.frame(X4)

vect <- list(
  formula = formula1, outgraphs = FALSE, data = X4,
  estimator = estimate.logic.lm,
  estimator.args = list(data = data.example, n = 100, m = n_col),
  recalc_margin = 1000L, save.beta = FALSE, interact = TRUE,
  relations = c("", "lgx2", "cos", "sigmoid", "tanh", "atan", "erf"),
  relations.prob = c(0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
  interact.param = list(
    allow_offsprings = 1, mutation_rate = 250, last.mutation = 1000,
    max.tree.size = 4, Nvars.max = 40, p.allow.replace = 0.7,
    p.allow.tree = 0.2, p.nor = 0, p.and = 0.9
  ), n.models = 20000, unique = TRUE, max.cpu = n_cores, max.cpu.glob = n_cores,
  create.table = FALSE, create.hash = TRUE, pseudo.paral = TRUE,
  burn.in = 50, print.freq = 1000,
  advanced.param = list(
    max.N.glob = 10L, min.N.glob = 5L, max.N = 3L, min.N = 1L, printable = FALSE
  )
)

params <- list(vect)[rep(1, M)]

for (i in seq_len(M)) {
  params[[i]]$cpu <- i
  params[[i]]$NM <- n_row
  params[[i]]$simlen <- 21
}

results <- parall.gmj(X = params, M = n_cores)

test_that("parall.gmj output matches version 1.4.3", {
  expect_length(results, M)
  for (i in seq_len(M)) {
    expect_length(results[[i]], 4L)
    expect_named(results[[i]], c("post.populi", "p.post", "cterm", "fparam"))
    expect_length(results[[i]][[2]], 40L)
    expect_true(all(results[[i]][[2]] <= 1.0))
    expect_length(results[[i]][[3]], 1L)
    expect_length(results[[i]][[4]], 40L)
  }
})
