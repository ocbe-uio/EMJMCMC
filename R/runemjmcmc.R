#' @title Mode jumping MJMCMC or Genetically Modified Mode jumping MCMC or Reversible Genetically Modified Mode jumping MCMC
#' for variable selection, Bayesian model averaging and feature engineering
#' @description A function that creates an EMJMCMC2016 object with
#' specified values of some parameters and default values of other parameters.
#' @details The algorithm is an extended Metropolis-Hastings algorithm
#' (or its Genetically modified version) mixing single site changes with
#' occasionally large jumps. The models are described through the gamma vector,
#' a binary vector indicating which variables that are included in the model.
#' @param formula a typical formula for specifying a model with all potential covariates included
#' @param data a data frame containing both covariates and response
#' @param secondary a character vector of names other covariates excluded from those defined in formula (relevant for GMJMCMC only)
#' @param latnames a character vector of names other covariates excluded from populations of GMJMCMC, for example for continuous covariates to be combined with BLR (relevant for GMJMCMC only) or the names of latent Gaussian variables to be selected in BGNLMM
#' @param save.beta a boolean parameter defining if beta coefficients for the models should be stored (must be set to TRUE if one is interested in predictions)
#' @param deep.method an integer in \{1, 2, 3, 4\} defining the method of estimating the alpha parameters of BGNLM, details to be found in https://www.jair.org/index.php/jair/article/view/13047
#' @param estimator a function returning a list with marginal likelihood, waic, dic and coefficients of the addressed model. The list should be of a format: list(mlik = mlik,waic = waic , dic = dic,summary.fixed =list(mean = coefficients))
#' @param estimator.args a list of arguments of estimator functions to be used (formula parameter has to be omitted, see the example)
#' @param n.models maximal number of models to be estimated during the search
#' @param unique defines whether n.models allows repetitions of the same models (unique=FALSE) or not (unique=TRUE)
#' @param locstop.nd Defines whether local greedy optimizers stop at the first local optima found (locstop.nd=TRUE) or not (locstop.nd=FALSE)
#' @param latent a latent random field to be addressed (to be specifically used when estimator = INLA, currently unsupported)
#' @param create.table a Boolean variable defining if a big.memory based hash table (only available for MJMCMC with no feature engineering, allows data sharing between CPUs) or the original R hash data structure (available for all algorithm, does not allow data sharing between CPUs) is used for storing of the results
#' @param hash.length a parameter defining hash size for the big.memory based hash table as 2^hash.length (only relevant when create.table = TRUE)
#' @param pseudo.paral defines if lapply or mclapply is used for local vectorized computations within the chain (can only be TRUE if create.table= TRUE)
#' @param max.cpu maximal number of CPUs in MJMCMC when within chain parallelization is allowed pseudo.paral = FALSE
#' @param max.cpu.glob maximal number of CPUs in global moves in MJMCMC when within chain parallelization is allowed pseudo.paral = FALSE
#' @param presearch a boolean parameter defining if greedy forward and backward regression steps are used for initialization of initial approximations of marginal inclusion probabilities
#' @param locstop a boolean parameter defining if the presearch is stopped at the first local extremum visited
#' @param interact a boolean parameter defining if feature engineering is allowed in the search
#' @param relations a vector of allowed modification functions (only relevant when feature engineering is enabled by means of interact = TRUE)
#' @param relations.prob probability distribution of addressing modifications defined in relations parameter (both vectors must be of the same length)
#' @param gen.prob a vector of probabilities for different operators in GMJMCMC or RGMJMCMC in the deep regression context (hence only relevant if \code{interact.param$allow_offsprings} is either 3 or 4)
#' @param pool.cross a parameter defining the probability of addressing covariates from the current pool of covariates in GMJMCMC (covariates from the set of filtered covariates can be addressed with probability 1-pool.cross) (only relevant when interact = TRUE)
#' @param p.add a default marginal inclusion probability parameter to be changed during the search to the true value
#' @param p.add.default a parameter defining sparsity after filtrations in GMJMCMC as initial marginal inclusion probabilities vector for parameters in the current pool
#' @param p.epsilon a parameter to define minimal deviations from 0 and 1 probabilities when allowing adaptive MCMC based on marginal inclusion probabilities
#' @param del.sigma a parameter describing probability of deleting each of the function from the selected feature in the reduction operator(only relevant for the deep regression models context)
#' @param pool.cor.prob a boolean parameter indicating if inclusion of the filtered covariates during mutations are based on probabilities proportional to the absolute values of correlations of these parameters and the observations (should not be addressed for multivariate observations, e.g. survival studies with Cox regression)
#' @param interact.param a list of parameters for GMJMCMC, where allow_offsprings is 1 for logic regression context, 2 for the old version of GMJMCMC for deep regressions, 3 for the new version of GMJMCMC for deep regressions and 4 for the RGMJMCMC for the deep regressions; mutation_rate defines how often changes of the search space are allowed in terms of the number of MJMCMC iterations per search space; last.mutation defines the iteration after which changes of search space are no longer allowed; max.tree.size is a parameter defining maximal depth of features; Nvars.max is a parameter defining maximal number of covariates in the search space after the first filtration; p.allow.replace is a parameter defining the upper bound on the probability allowing the replacement of corresponding features with marginal inclusion probabilities below it; p.allow.tree is a lower bound for the probability of not being filtered out after initializing steps of MJMCMC in GMJMCMC; p.nor is a parameter for not operator in the logic regression context (allow_offsprings==1); p.and = is the probability of & crossover in the logic regression context (allow_offsprings==1)
#' @param prand probability of changes of components in randomization kernels of RGMJMCMC
#' @param keep.origin a boolean parameter defining if the initially unfiltered covariates can leave the search space afterwards (TRUE) or not (FALSE)
#' @param sup.large.n omitted currently
#' @param recalc_margin a parameter defining how often marginal inclusion probabilities would be recalculated
#' @param create.hash a parameter defining if by default the results are stored in a hash table
#' @param interact.order omitted currently
#' @param burn.in number of burn-in steps for (R)(G)MJMCMC
#' @param eps omitted, not to be changed
#' @param max.time maximal time for the run of (R)(G)MJMCMC algorithm in minutes
#' @param max.it maximal number of (R)(G)MJMCMC iterations
#' @param print.freq printing frequency of the intermediate results
#' @param outgraphs a boolean variable defining if the graphics on the marginal inclusion probabilities should be drawn (must not be used inside mclapply wrapper of runemjmcmc since otherwise errors can occur)
#' @param advanced.param omitted currently
#' @param distrib_of_neighbourhoods a matrix defining probability distribution on 7 types of neighbourhoods within 4 possible local search strategies as well as within global moves
#' @param distrib_of_proposals probability distribution up to a constant of proportionality for addressing different local search strategies after large jumps or no large jumps (5th component)
#' @param quiet defaults to \code{FALSE}. If \code{TRUE}, prints intermediate
#' messages
#' @details See Hubin & Storvik (2016),Hubin, Storvik & Frommlet (2017),
#' Hubin & Storvik (2017) details. The local optimization is performed through
#' stepwise search within a neighborhood in the current gamma vector, allowing
#' one component to be changed at a time.
#' @return a list containing
#' \describe{
#'  \item{p.post}{a vector of posterior probabilities of the final vector of active covariates (features)}
#'  \item{m.post}{a vector of posterior probabilities of the models from the search space induced by the final vector of active covariates (features)}
#'  \item{s.mass}{sum of marginal likelihoods times the priors from the explored part of the search space induced by the final vector of active covariates (features)}
#' }
#' @references Hubin & Storvik (2016),Hubin, Storvik & Frommlet (2017), Hubin & Storvik (2017)
#' @author Aliaksandr Hubin
#' @seealso global objects statistics1 (if create.table== TRUE) or hashStat (if create.table== FALSE) contain all marginal likelihoods and two other model selection criteria as well as all of the beta coefficients for the models (if save.beta== TRUE)
#' @example /inst/examples/runemjmcmc_example.R
#' @keywords methods models
#' @export
runemjmcmc<-function(
  formula, data, secondary = vector(mode="character", length=0), latnames="",
  estimator,estimator.args = "list",n.models,p.add.default = 1,p.add = 0.5,
  unique = FALSE,save.beta= FALSE, locstop.nd = FALSE, latent="",max.cpu=4,max.cpu.glob=2,
  create.table= TRUE, hash.length = 20, presearch= TRUE, locstop = FALSE ,pseudo.paral = FALSE,
  interact = FALSE,deep.method =1,
  relations = c("","sin","cos","sigmoid","tanh","atan","erf"),
  relations.prob =c(0.4,0.1,0.1,0.1,0.1,0.1,0.1),gen.prob = c(1,10,5,1,0),
  pool.cross = 0.9,p.epsilon = 0.0001, del.sigma = 0.5,pool.cor.prob = FALSE,
  interact.param=list(allow_offsprings=2,mutation_rate = 100,last.mutation=2000,
  max.tree.size = 10000, Nvars.max = 100, p.allow.replace = 0.7,
  p.allow.tree=0.1,p.nor=0.3,p.and = 0.7), prand = 0.01,keep.origin = TRUE,
  sup.large.n = 5000, recalc_margin = 2^10, create.hash= FALSE,interact.order=1,
  burn.in=1, eps = 10^6, max.time = 120,max.it = 25000, print.freq = 100,
  outgraphs= FALSE,advanced.param=NULL,
  distrib_of_neighbourhoods=t(array(data = c(7.6651604,16.773326,14.541629,12.839445,2.964227,13.048343,7.165434, 0.9936905,15.942490,11.040131,3.200394,15.349051,5.466632,14.676458, 1.5184551,9.285762,6.125034,3.627547,13.343413,2.923767,15.318774, 14.5295380,1.521960,11.804457,5.070282,6.934380,10.578945,12.455602, 6.0826035,2.453729,14.340435,14.863495,1.028312,12.685017,13.806295),dim = c(7,5))),
  distrib_of_proposals = c(76.91870,71.25264,87.68184,60.55921,15812.39852),
  quiet = TRUE)
{

#first create the object
global_env <- as.environment(1L)
assign("data.example",data, envir=global_env)
variables <- simplify.formula(formula,names(data.example))
assign("fparam.example",variables$fparam, envir=global_env)
assign("fobserved.example",variables$fobserved, envir=global_env)

#for(i in 1:length(fparam.example))
#{
#  fparam.example[i]<<-paste("I(",variables$fparam[i],")",sep = "")
#}
fparam.tmp<- as.vector(sapply(FUN = paste,"I(",variables$fparam,")",sep="")[,1])
if(latnames[1]!="")
  assign("fparam.example", c(fparam.tmp,latnames), envir=global_env)
else
  assign("fparam.example", fparam.tmp, envir=global_env)
assign("mySearch",methods::new(structure("EMJMCMC2016", package = "EMJMCMC")), envir=global_env)
if(length(secondary)>0)
  mySearch$filtered <- sapply(FUN = paste,"I(",secondary,")",sep="")
mySearch$estimator <- estimator
mySearch$latnames <- latnames
mySearch$estimator.args <<- estimator.args
mySearch$latent.formula <- latent
mySearch$save.beta <- save.beta
mySearch$prand <- prand
mySearch$p.add <- array(p.add,length(fparam.example))
mySearch$p.add.default <- p.add.default
mySearch$recalc.margin <<- as.integer(recalc_margin)
mySearch$max.cpu <<- as.integer(max.cpu)
mySearch$locstop.nd <- locstop.nd
mySearch$pool.cor.prob <- pool.cor.prob
mySearch$sup.large.n <- as.integer(sup.large.n)
mySearch$max.cpu.glob <<- as.integer(max.cpu.glob)
mySearch$deep.method <- as.integer(deep.method)
if(interact)
{
  mySearch$allow_offsprings <<- as.integer(interact.param$allow_offsprings)
  mySearch$mutation_rate <<- as.integer(interact.param$mutation_rate)
  mySearch$Nvars.max <<- as.integer(interact.param$Nvars.max)
  mySearch$max.tree.size <- as.integer(interact.param$max.tree.size)
  mySearch$p.allow.replace <- interact.param$p.allow.replace
  mySearch$p.allow.tree <<-  interact.param$p.allow.tree
  mySearch$p.epsilon <- p.epsilon
  mySearch$keep.origin <- keep.origin
  mySearch$sigmas<<-relations
  mySearch$sigmas.prob <<- relations.prob
  mySearch$del.sigma <- del.sigma
  mySearch$pool.cross <- pool.cross
  mySearch$gen.prob<<-gen.prob
  mySearch$p.nor <- interact.param$p.nor
  mySearch$p.and <- interact.param$p.and
  mySearch$last.mutation <- as.integer(interact.param$last.mutation)
}

if(!is.null(advanced.param))
{
  mySearch$max.N.glob<<-as.integer(advanced.param$max.N.glob)
  mySearch$min.N.glob<<-as.integer(advanced.param$min.N.glob)
  mySearch$max.N<<-as.integer(advanced.param$max.N)
  mySearch$min.N <- as.integer(advanced.param$min.N)
  mySearch$printable.opt <- advanced.param$printable
}

if(exists("hashStat"))
{
  hash::clear(hashStat)
  remove(hashStat,envir=global_env)
}
if(exists("statistics1"))
{
  remove(statistics,envir=global_env )
  remove(statistics1,envir=global_env)
}
if(exists("hash.keys1"))
{
  remove(hash.keys,envir=global_env)
  remove(hash.keys1,envir=global_env)
}

if(create.table)
{
  if(pseudo.paral) mySearch$parallelize <- lapply
  #carry the search (training out)
  assign("statistics1",bigmemory::big.matrix(nrow = 2 ^min((length(fparam.example)),hash.length)+1, ncol =  16+length(fparam.example)*save.beta,init = NA, type = "double"), envir=global_env)
  assign("statistics",bigmemory::describe(statistics1), envir=global_env)
  mySearch$g.results[4,1] <- 0
  mySearch$g.results[4,2] <- 0
  mySearch$p.add <- array(data = 0.5,dim = length(fparam.example))
  if((length(fparam.example))>20)
  {
    mySearch$hash.length <- as.integer(hash.length)
    mySearch$double.hashing <- TRUE
    assign(hash.keys1, bigmemory::big.matrix(nrow = 2 ^(hash.length)+1, ncol = length(fparam.example),init = 0, type = "char"), global_env)
    assign(hash.keys, bigmemory::describe(hash.keys1), global_env)
  }

}else if(create.hash)
{

  assign("hashStat", hash::hash(), envir=global_env)
  mySearch$parallelize <<- lapply
  mySearch$hash.length <- as.integer(20)
  mySearch$double.hashing<<- FALSE
}
# now as the object is created run the algorithm
initsol=stats::rbinom(n = length(fparam.example),size = 1,prob = 0.5)
if(unique)
  resm <- mySearch$modejumping_mcmc(list(varcur=initsol,locstop=locstop,presearch=presearch,statid=5, distrib_of_proposals =distrib_of_proposals,distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = eps, trit = n.models*100, trest = n.models, burnin = burn.in, max.time = max.time, maxit = max.it, print.freq = print.freq))
else
  resm <- mySearch$modejumping_mcmc(list(varcur=initsol,locstop=locstop,presearch=presearch,statid=5, distrib_of_proposals =distrib_of_proposals,distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = eps, trit =  n.models, trest = n.models*100, burnin = burn.in, max.time = max.time, maxit = max.it, print.freq = print.freq))
ppp<-1
if (!quiet) message("MJMCMC is completed")
if(create.table)
{
  ppp<-mySearch$post_proceed_results(statistics1 = statistics1)
  truth = ppp$p.post # make sure it is equal to Truth column from the article
  truth.m = ppp$m.post
  truth.prob = ppp$s.mass
  ordering = sort(ppp$p.post,index.return= TRUE)
}
else if(create.hash)
{
  ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
  truth = ppp$p.post # make sure it is equal to Truth column from the article
  truth.m = ppp$m.post
  truth.prob = ppp$s.mass
  ordering = sort(ppp$p.post,index.return= TRUE)
}
if (!quiet) {
  message("Post Proceed Results")
  cat(
    "pi truth",
    sprintf("%.10f",truth[ordering$ix]),
    sprintf(fparam.example[ordering$ix]),
    sep = "\n"
  )
}

if(outgraphs)
{
  withr::local_par(mar = c(10,4,4,2) + 4.1)
  graphics::barplot(resm$bayes.results$p.post,density = 46,border="black",main = "Marginal Inclusion (RM)",ylab="Probability",names.arg = mySearch$fparam,las=2)
  graphics::barplot(resm$p.post,density = 46,border="black",main = "Marginal Inclusion (MC)",ylab="Probability",names.arg = mySearch$fparam,las=2)
}

return(ppp)
}
