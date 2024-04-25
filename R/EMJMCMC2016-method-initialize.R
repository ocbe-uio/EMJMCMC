# class constructor
EMJMCMC2016$methods(initialize = function(estimator.function = estimate.gamma.cpen, estimator.args.list = list(
                        family = "gaussian", data = data.example, control.fixed = list(prec = list(default = 0.00001), prec.intercept = 0.00001, mean = list(default = 0), mean.intercept = 0),
                        control.family = list(hyper = list(prec = list(prior = "loggamma", param = c(0.00001, 0.00001), initial = 0))),
                        control.compute = list(dic = TRUE, waic = TRUE, mlik = TRUE)
                      ), search.args.list = NULL, latent.formula = "") {
  estimator <<- estimator.function
  estimator.args <<- estimator.args.list
  latent.formula <<- latent.formula
  temp.file <- gsub("\\W", "", tempfile(tmpdir = ""))
  g.results <<- big.matrix(
    nrow = 4, ncol = 2,
    backingpath = tempdir(),
    backingfile = paste0(temp.file, ".bak"),
    descriptorfile = paste0(temp.file, ".desc"),
  )
  g.results[1, 1] <<- -Inf
  g.results[1, 2] <<- 1
  g.results[2, 1] <<- Inf
  g.results[2, 2] <<- 1
  g.results[3, 1] <<- Inf
  g.results[3, 2] <<- 1
  g.results[4, 1] <<- 0
  g.results[4, 2] <<- 0

  if (is.null(search.args.list)) {
    max.cpu <<- as.integer(Nvars * 0.05 + 1)
    objective <<- as.integer(1)
    if (Sys.info()["sysname"] == "Windows") {
      parallelize <<- lapply
      parallelize.global <<- lapply
      parallelize.hyper <<- lapply
    } else {
      parallelize <<- parallel::mclapply
      parallelize.global <<- parallel::mclapply
      parallelize.hyper <<- parallel::mclapply
    }
    Nvars <<- as.integer(length(fparam.example))
    min.N <<- as.integer(Nvars / 6)
    min.N.glob <<- as.integer(Nvars / 3)
    max.N.glob <<- as.integer(Nvars / 2)
    max.N <<- as.integer(Nvars / 5)
    switch.type.glob <<- as.integer(2)
    min.N.randomize <<- as.integer(1)
    max.N.randomize <<- as.integer(1)
    deep.method <<- as.integer(1)
    type.randomize <<- as.integer(3)
    pool.cor.prob <<- FALSE
    prand <<- 0.01
    max.cpu.glob <<- as.integer(Nvars * 0.05 + 1)
    max.cpu.hyper <<- as.integer(2)
    sup.large.n <<- as.integer(1000)
    save.beta <<- FALSE
    filtered <<- vector(mode = "character", length = 0)
    printable.opt <<- FALSE
    keep.origin <<- FALSE
    thin_rate <<- as.integer(-1)
    p.allow.tree <<- 0.6
    p.epsilon <<- 0.0001
    latnames <<- ""
    p.add.default <<- 1
    p.allow.replace <<- 0.3
    sigmas <<- c("", "sin", "cos", "sigmoid", "tanh", "atan", "erf")
    sigmas.prob <<- c(0.4, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
    del.sigma <<- 0.5
    pool.cross <<- 0.9
    gen.prob <<- c(1, 1, 1, 1, 1)
    p.nor <<- 0.3
    p.and <<- 0.7
    max.tree.size <<- as.integer(15)
    Nvars.max <<- as.integer(Nvars)
    Nvars.init <<- as.integer(Nvars)
    allow_offsprings <<- as.integer(0)
    mutation_rate <<- as.integer(100)
    locstop.nd <<- FALSE
    double.hashing <<- (Nvars > 20)
    hash.length <<- as.integer(25)
    filtered <<- vector(mode = "character", length = 0)
    aa <<- 0.9
    cc <<- 0.0
    M.nd <<- as.integer(min(Nvars, 100))
    M.mcmc <<- as.integer(5)
    SA.param <<- list(t.min = 0.0001, t.init = 10, dt = 3, M = as.integer(min(Nvars / 5 + 1, 20)))
    fobserved <<- fobserved.example
    switch.type <<- as.integer(2)
    n.size <<- as.integer(10)
    LocImprove <<- as.array(c(50, 50, 50, 50, 150))
    isobsbinary <<- as.array(0:(length(fparam.example) - 1))
    fparam <<- fparam.example
    fparam.pool <<- fparam.example
    p.add <<- array(data = 0.5, dim = Nvars)
    if (exists("statistics")) {
      recalc.margin <<- 2^Nvars
    } else if (exists("statistics1")) {
      recalc.margin <<- 2^Nvars
    } else {
      recalc.margin <<- 2^Nvars
    }
    last.mutation <<- as.integer(2^(Nvars / 2) * 0.01)
    seed <<- as.integer(stats::runif(n = 1, min = 1, max = 10000))
    p.prior <<- stats::runif(n = Nvars, min = 0.5, max = 0.5)
  } else {
    max.cpu <<- as.integer(search.args.list$max.cpu)
    objective <<- as.integer(search.args.list$objective)
    parallelize <<- search.args.list$parallelize
    latnames <<- search.args.list$latnames
    parallelize.global <<- search.args.list$parallelize.global
    parallelize.hyper <<- search.args.list$parallelize.hyper
    p.prior <<- search.args.list$p.prior
    min.N <<- as.integer(search.args.list$min.N)
    printable.opt <<- search.args.list$printable.opt
    min.N.glob <<- as.integer(search.args.list$min.N.glob)
    max.N.glob <<- as.integer(search.args.list$max.N.glob)
    switch.type.glob <<- as.integer(search.args.list$switch.type.glob)
    min.N.randomize <<- as.integer(search.args.list$min.N.randomize)
    max.N.randomize <<- as.integer(search.args.list$max.N.randomize)
    type.randomize <<- as.integer(search.args.list$type.randomize)
    max.cpu.glob <<- as.integer(search.args.list$max.cpu.glob)
    locstop.nd <<- search.args.list$locstop.nd
    max.cpu.hyper <<- as.integer(search.args.list$max.cpu.hyper)
    save.beta <<- search.args.list$save.beta
    aa <<- search.args.list$lambda.a
    prand <<- search.args.list$prand
    p.add.default <<- search.args.list$ p.add.default
    sup.large.n <<- search.args.list$sup.large.n
    thin_rate <<- search.args.list$thin_rate
    keep.origin <<- search.args.list$keep.origin
    cc <<- search.args.list$lambda.c
    pool.cor.prob <<- search.args.list$pool.cor.prob
    M.nd <<- as.integer(search.args.list$stepsGreedy)
    M.mcmc <<- as.integer(search.args.list$stepsLocMCMC)
    SA.param <<- search.args.list$SA.params
    fobserved <<- search.args.list$fobserved
    switch.type <<- as.integer(search.args.list$fswitch.type)
    n.size <<- as.integer(search.args.list$n.size)
    LocImprove <<- as.array(search.args.list$prior.optimizer.freq)
    max.N <<- as.integer(search.args.list$max.N)
    fparam <<- search.args.list$fparam
    fparam.pool <<- search.args.list$fparam
    isobsbinary <<- as.array(0:(length(fparam) - 1))
    p.add <<- as.array(search.args.list$p.add)
    recalc.margin <<- search.args.list$recalc.margin
    Nvars <<- as.integer(length(fparam))
    seed <<- search.args.list$seed
    max.tree.size <<- as.integer(search.args.list$max.tree.size)
    Nvars.max <<- as.integer(search.args.list$Nvars.max)
    Nvars.init <<- as.integer(search.args.list$Nvars)
    allow_offsprings <<- as.integer(search.args.list$allow_offsprings)
    mutation_rate <<- as.integer(search.args.list$mutation_rate)
    p.allow.tree <<- search.args.list$p.allow.tree
    p.epsilon <<- search.args.list$p.epsilon
    p.allow.replace <<- search.args.list$p.allow.replace
    last.mutation <<- as.integer(search.args.list$last.mutation)
    p.nor <<- search.args.list$p.nor
    p.and <<- search.args.list$p.and
    deep.method <<- as.integer(search.args.list$deep.method)
    sigmas <<- search.args.list$sigmas
    sigmas.prob <<- search.args.list$sigmas.prob
    del.sigma <<- search.args.list$del.sigma
    pool.cross <<- search.args.list$pool.cross
    gen.prob <<- search.args.list$gen.prob
    double.hashing <<- search.args.list$double.hashing
    hash.length <<- as.integer(search.args.list$hash.length)
  }
})
