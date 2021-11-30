# TODO #8: fill empty tests

# inference

X <- read.csv("inst/extdata/exa1.csv")
data.example <- as.data.frame(X)

# specify the initial formula
formula1 <- as.formula(
  paste(colnames(X)[5], "~ 1 +", paste0(colnames(X)[-5], collapse = "+"))
)

# a set of nonlinearities that will be used in the DBRM model
sini <- function(x) sin(x / 180 * pi)
expi <- function(x) exp(-abs(x))
logi <- function(x) log(abs(x) + 1)
troot <- function(x) abs(x)^(1 / 3)
to23 <- function(x) abs(x)^(2.3)
to35 <- function(x) abs(x)^(3.5)

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
res1 <- pinferunemjmcmc(
  n.cores = M, report.level = 0.5, num.mod.best = NM, simplify = TRUE,
  runemjmcmc.params = list(
    formula = formula1, data = data.example, estimator = estimate.gamma.cpen,
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

test_that("pinferunemjmcmc output matches version 1.4.3", {

})
