X <- read.csv(system.file("extdata", "exa1.csv", package="EMJMCMC"))
data.example <- as.data.frame(X)

# specify the initial formula
formula1 <- as.formula(
  paste(colnames(X)[5], "~ 1 +", paste0(colnames(X)[-5], collapse = "+"))
)

M <- 2 # define the number or cpus
NM <- 100 # define the size of the simulated samples
compmax <- 16 # define \k_{max} + 1 from the paper
th <- (10)^(-5) # define treshold for preinclusion of the tree into the analysis
thf <- 0.05 # final treshold on the posterior marginal prob for reporting a tree
# specify tuning parameters of the algorithm for exploring DBRM of interest
# notice that allow_offsprings=3 corresponds to the GMJMCMC runs and
# allow_offsprings=4 -to the RGMJMCMC runs
set.seed(9239838)
res1 <- suppressMessages(
  pinferunemjmcmc(
    n.cores = M, report.level = 0.5, num.mod.best = NM, simplify = TRUE,
    runemjmcmc.params = list(
      formula = formula1, data = data.example, estimator = estimate.gamma.cpen,
      estimator.args = list(data = data.example), recalc_margin = 249,
      save.beta = FALSE, interact = TRUE, outgraphs = FALSE,
      relations = c("to23", "expi", "logi2", "to35", "sini", "troot", "sigmoid"),
      relations.prob = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
      interact.param = list(
        allow_offsprings = 3, mutation_rate = 250,
        last.mutation = 10000, max.tree.size = 5, Nvars.max = 15,
        p.allow.replace = 0.9, p.allow.tree = 0.01, p.nor = 0.9, p.and = 0.9
      ),
      n.models = 1000, unique = TRUE, max.cpu = 4, max.cpu.glob = 4,
      create.table = FALSE, create.hash = TRUE, pseudo.paral = TRUE,
      burn.in = 10, print.freq = 1000,
      advanced.param = list(
        max.N.glob = as.integer(10),
        min.N.glob = as.integer(5),
        max.N = as.integer(3),
        min.N = as.integer(1),
        printable = FALSE
      )
    )
  )
)

test_that("pinferunemjmcmc output matches version 1.4.3", {
  expect_named(
    res1,
    c("feat.stat", "predictions", "allposteriors", "threads.stats")
  )
  expect_length(res1, 4)
  expect_equal(ncol(res1$feat.stat), 2L)
  expect_equal(mean(res1$allposteriors$posterior), 0.3, tolerance = 1e-1)
  expect_true(all(res1$threads.stats[[1]]$p.post >= 0))
  expect_true(all(res1$threads.stats[[1]]$p.post <= 1))
})
