X <- read.csv(system.file("extdata", "exa1.csv", package="EMJMCMC"))
data.example <- as.data.frame(X)

# specify the initial formula
formula1 <- as.formula(
  paste(colnames(X)[5], "~ 1 +", paste0(colnames(X)[-5], collapse = "+"))
)

# a set of nonlinearities that will be used in the DBRM model as well as the
# "relations" argument of pinferunemjmcmc
sini <- function(x) sin(x / 180 * pi)
expi <- function(x) exp(-abs(x))
logi <- function(x) log(abs(x) + 1)
troot <- function(x) abs(x)^(1 / 3)
to23 <- function(x) abs(x)^(2.3)
to35 <- function(x) abs(x)^(3.5)

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
      relations = c("to23", "expi", "logi", "to35", "sini", "troot", "sigmoid"),
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
  expect_length(res1$feat.stat, 10)
  expect_equal(mean(res1$allposteriors$posterior), 0.309510, tolerance = 1e-4)
  expect_equal(mean(res1$threads.stats[[1]]$p.post), 0.391782, tolerance = 1e-4)
})
