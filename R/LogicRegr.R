LogicRegr = function(formula, data, family = "Gaussian",prior = "J",report.level = 0.5, d = 20, cmax = 5, kmax = 20, p.and = 0.9, p.not = 0.05, p.surv = 0.1,ncores = -1, n.mods = 1000 ,advanced = list(presearch = T,locstop = F ,estimator = estimate.logic.bern.tCCH,estimator.args =  list(data = data.example,n = 1000, m = 50,r=1),recalc_margin = 250, save.beta = F,interact = T,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=1,mutation_rate = 300,last.mutation = 5000, max.tree.size = 1, Nvars.max = 100,p.allow.replace=0.9,p.allow.tree=0.2,p.nor=0.2,p.and = 1),n.models = 10000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,outgraphs=F,print.freq = 1000,advanced.param = list(
  max.N.glob=as.integer(10),
  min.N.glob=as.integer(5),
  max.N=as.integer(3),
  min.N=as.integer(1),
  printable = F)))
{
  data.example = data
  advanced$formula = formula
  advanced$data = data
  advanced$interact.param$Nvars.max = d
  advanced$interact.param$ max.tree.size= cmax - 1
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

  return(pinferunemjmcmc(n.cores = ncores,report.level =  report.level, simplify = T, num.mod.best = n.mods, predict = F, runemjmcmc.params = advanced))

}
