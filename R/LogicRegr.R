#' @title A wrapper for running the Bayesian logic regression based inference
#' in a easy to use way
#' @description A wrapper for running the Bayesian logic regression based
#' inference in a easy to use way
#' @param formula a formula object for the model to be addressed
#' @param data a data frame object containing variables and observations
#' corresponding to the formula used
#' @param family a string taking values of either "Gaussian" or "Bernoulli"
#' correspodning to the linear or logistic Bayesian logic regression contexts
#' @param prior character values "J" or "G" corresponing either to Jeffey's
#' or robust g prior
#' @param report.level a numeric value in (0,1) specifying the treshold for
#' detections based on the marginal inclusion probabilities
#' @param d population size for the GMJMCMC algorithm
#' @param cmax the maximal allowed depth of logical expressions to be considered
#' @param kmax the maximal number of logical expressions per model
#' @param p.and probability of AND parameter of GMJMCMC algorithm
#' @param p.not probability of applying logical NOT in GMJMCMC algorithm
#' @param p.surv minimal survival probabilities for the features to be allowed to enter the next population
#' @param ncores the maximal number of cores (and GMJMCMC threads) to be
#' addressed in the analysis
#' @param n.mods the number of the best models in the thread to calculate
#' marginal inclusion probabilities
#' @param advanced should only be adrresed by experienced users to tune advanced
#' parameters of GMJMCMC, advanced corresponds to the vector of tuning
#' parameters of runemjmcmc function
#' @param print.freq printing frequency of the intermediate results
#' @return a list of
#' \describe{
#'  \item{feat.stat}{detected logical expressions and their marginal inclusion
#'    probabilities}
#'  \item{predictions}{NULL currently, since LogrRegr function is not designed
#'    for predictions at the moment, which is still possible in its expert
#'    mother function pinferunemjmcmc}
#'  \item{allposteriors}{all visited by GMJMCMC logical expressions and their
#'    marginal inclusion probabilities}
#'  \item{threads.stats}{a vector of detailed outputs of individual ncores
#'    threads of GMJMCMC run}
#' }
#' @seealso runemjmcmc pinferunemjmcmc
#' @keywords methods models
#' @example /inst/examples/LogicRegr_example.R
#' @export
LogicRegr = function(
  formula, data, family = "Gaussian",prior = "J",report.level = 0.5, d = 20,
  cmax = 5, kmax = 20, p.and = 0.9, p.not = 0.05, p.surv = 0.1, ncores = -1,
  n.mods = 1000, print.freq = 1000L,
  advanced = list(
    presearch = TRUE,locstop = FALSE,
    estimator = estimate.logic.bern.tCCH,
    estimator.args =  list(data = data.example,n = 1000, m = 50,r=1),
    recalc_margin = 250, save.beta = FALSE, interact = TRUE,
    relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),
    relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),
    interact.param = list(
      allow_offsprings=1, mutation_rate = 300, last.mutation = 5000,
      max.tree.size = 1, Nvars.max = 100, p.allow.replace=0.9, p.allow.tree=0.2,
      p.nor=0.2, p.and = 1
    ),
    n.models = 10000, unique = TRUE, max.cpu = ncores, max.cpu.glob = ncores,
    create.table = FALSE, create.hash = TRUE, pseudo.paral = TRUE, burn.in = 50,
    outgraphs = FALSE, print.freq = print.freq,
    advanced.param = list(
      max.N.glob=as.integer(10),
      min.N.glob=as.integer(5),
      max.N=as.integer(3),
      min.N=as.integer(1),
      printable = FALSE
    )
  )
) {
  data.example = data
  advanced$formula = formula
  advanced$data = data
  advanced$interact.param$Nvars.max = d
  advanced$interact.param$max.tree.size = cmax - 1
  advanced$interact.param$p.and = p.and
  advanced$interact.param$p.nor = p.not
  advanced$interact.param$p.allow.tree = p.surv
  if(!prior %in% c("J","G"))
  {
    warning("Wrong prior supplied. J (for Jeffrey's) and G (for robust g) priors are allowd only. Setting J as default.")
    prior = "J"
  }
  if(!family %in% c("Gaussian","Bernoulli"))
  {
    warning("Wrong familty supplied. Gaussian and Bernoulli families are allowd only. Setting Gaussian as default.")
    family = "Gaussian"
  }
  if(family == "Gaussian")
  {
    if(prior == "J")
    {
      advanced$estimator = estimate.logic.lm.jef
      advanced$estimator.args = list(data = data, n = dim(data)[1], m =stringi::stri_count_fixed(as.character(formula)[3],"+"),k.max = kmax)
    }else{
      advanced$estimator = estimate.logic.lm.tCCH
      advanced$estimator.args = list(data = data, n = dim(data)[1], m =stringi::stri_count_fixed(as.character(formula)[3],"+"),k.max = kmax)
    }
  }else{

    if(prior == "J")
    {
      advanced$estimator = estimate.logic.bern
      advanced$estimator.args =  list(data = data, n = dim(data)[1], m =stringi::stri_count_fixed(as.character(formula)[3],"+"),k.max = kmax)
    }else{
      advanced$estimator = estimate.logic.bern.tCCH
      advanced$estimator.args =  list(data = data, n = dim(data)[1], m =stringi::stri_count_fixed(as.character(formula)[3],"+"),k.max = kmax)
    }

  }

  if(ncores<1)
    ncores = parallel::detectCores()

  return(
    pinferunemjmcmc(
      n.cores = ncores, mcgmj = mcgmjpar, report.level = report.level, simplify = TRUE,
      num.mod.best = n.mods, predict = FALSE, runemjmcmc.params = advanced
    )
  )
}
