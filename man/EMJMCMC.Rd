\name{runemjmcmc}
\alias{runemjmcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mode jumping MJMCMC or Genetically Modified Mode jumping MCMC for variable selection, Bayesian model averaging and feature engineering}
\description{The algorithm is an extended Metropolis-Hastings algorithm (or its Genetically modified version) mixing single site changes with occationally large jumps. The models are described through the gamma vector, a binary vector indicating which variables that are included in the model.}
\usage{
runemjmcmc<-function(formula, data, secondary = vector(mode="character", length=0),
                     estimator,estimator.args = "list",n.models, unique = F,save.beta=F, locstop.nd = F, latent="",max.cpu=4,max.cpu.glob=2,create.table=T, hash.length = 20, presearch=T, locstop =F ,pseudo.paral = F,interact = F,relations = c("","sin","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.1,0.1,0.1,0.1,0.1,0.1),gen.prob = c(1,10,5,1,1),pool.cross = 0.9,p.epsilon = 0.0001, del.sigma = 0.5,pool.cor.prob = F, interact.param=list(allow_offsprings=2,mutation_rate = 100,last.mutation=2000, max.tree.size = 10000, Nvars.max = 100, p.allow.replace = 0.7,p.allow.tree=0.1,p.nor=0.3,p.and = 0.7), prand = 0.01,keep.origin = T, sup.large.n = 5000, recalc_margin = 2^10, create.hash=F,interact.order=1,burn.in=100, eps = 10^6, max.time = 120,max.it = 25000, print.freq = 100,outgraphs=F,advanced.param=NULL, distrib_of_neighbourhoods,  distrib_of_proposals = c(76.91870,71.25264,87.68184,60.55921,15812.39852))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A typical formula for specifying a model with all potential covariates included}
  \item{data}{A data frame containing both covariates and response}
  \item{secondary}{A charactor vector of names other covariates excluded from those defined in formula (relevant for GMJMCMC only)}
  \item{estimator}{A function returning a list with marginal likelihood, waic, dic and coefficients of the addressed model. The list should be of a format: list(mlik = mlik,waic = waic , dic = dic,summary.fixed =list(mean = coefficients))}
  \item{estimator.args}{A list of arguments of estimator functions to be used (formula parameter has to be omitted, see the example)}
  \item{n.models}{Maximal number of models to be estimated during the search}
  \item{unique}{Defines whether n.models allows repetitions of the same models (unique=FALSE) or not (unique=TRUE)}
  \item{locstop.nd}{Defines whether local greedy optimizers stop at the first local optima found (locstop.nd=TRUE) or not (locstop.nd=FALSE)}
  \item{latent}{A latent random field to be adressed (to be specifucally used when estimator = INLA)}
  \item{create.table}{A boolean variable defining if a big.memory based hash table (only available for MJMCMC with no feature engineering, allows data sharing between CPUs) or the original R hash data structure (available for all algorithm, does not allow data sharing between CPUs) is used for storing of the results}
  \item{hash.length}{A paramter defining hash size for the big.memory based hash table as 2^hash.length (only relevant when create.table = T)}
  \item{pseudo.paral}{Defines if lapply or mclapply is used for local vectorized computations within the chain (can only be TRUE if create.table=T)}
  \item{max.cpu}{Maximal number of cpus in MJMCMC when within chain parallelization is allowed pseudo.paral = F}
  \item{max.cpu.glob}{Maximal number of cpus in global moves in MJMCMC when within chain parallelization is allowed pseudo.paral = F}
  \item{presearch}{A boolean parameter defining if greedy forward and backward regression steps are used for initialization of initial approximations of marginal inclusion probabilities}
  \item{locstop}{A boolean parameter defining if the presearch is stopped at the first local exremum visited}
  \item{interact}{A boolean parameter defining if feature engineering is allowed in the search}
  \item{relations}{A vector of allowed modification functions (only relevant when feature engineering is enabled by means of interact = T)}
  \item{relations.prob}{Probability distribution of addressing modifications defined in relations parameter (both vectors must be of the same length)}
  \item{gen.prob}{A vector of probabilities for different operators in GMJMCMC or RGMJMCMC in the deep regression context (hence only relevant if interact.param$allow_offsprings is either 3 or 4)}
  \item{pool.cross}{A parameter defining the probability of adressing covariates from the current pool of covariates in GMJMCMC (covariates from the set of filtered covariates can be addressed with probability 1-pool.cross) (only relevant when interact = TRUE)}
  \item{p.epsilon}{A parameter to define minimal deviations from 0 and 1 probabilities when allowing adaptive MCMC based on marginal inclusion probabilities}
  \item{del.sigma}{A parameter decribing probability of deleting each of the function from the selected faeture in the reduction operator(only relevant for the deep regression models context)}
  \item{pool.cor.prob}{A boolean parameter indicating if inclusion of the filtered covariates during mutations are based on probabilities proportional to the absoulute values of correlations of these parameters and the observations (should not be addressed for multivariate observations, e.g. survival studies with Cox regression)}
  \item{interact.param}{a list of parameters for GMJMCMC, where allow_offsprings is 1 for logic regression context, 2 for the old version of GMJMCMC for deep regressions, 3 for the new version of GMJMCMC for deep regressions and 4 for the RGMJMCMC for the deep regressions; mutation_rate defines how often changes of the search space are allowed in terms of the number of MJMCMC iterations per search space; last.mutation defines the iteration after which changes of search space are no longer allowed; max.tree.size is a parameter defining maximal depth of features; Nvars.max is a parameter defining maximal number of covariates in the search space after the first filtration; p.allow.replace is a parameter defining the upepr bound on the probability allowing the replacement of correspondin features with marginal inclusion probabilities below it; p.allow.tree is a lower bound for the probability of not being filtered out after initializing steps of MJMCMC in GMJMCMC; p.nor is a parameter for not operator in the logic regression context (allow_offsprings==1); p.and = is the probability of & crossover in the logic regression context (allow_offsprings==1)}
  \item{prand}{Probability of changes of components in randomization kernels of RGMJMCMC}
  \item{keep.origin}{A boolean parameter defining if the initially unfiltered covariates can leave the search space afterwards (TRUE) or not (FALSE)}
  \item{sup.large.n}{Omitted currently}
  \item{recalc_margin}{A parameter defining how often marginal inclusion probabilities whould be recalculated}
  \item{create.hash}{A parameter defining if by default the results are stored in a hash table}
  \item{interact.order}{Omitted currently}
  \item{burn.in}{Number of burn-in steps for (R)(G)MJMCMC}
  \item{eps}{Currently omitted, not to be changed}
  \item{max.time}{Maximal time for the run of (R)(G)MJMCMC algorithm in minutes}
  \item{max.it}{Maximal number of (R)(G)MJMCMC iterations}
  \item{print.freq}{Printing frequency of the intermediate results}
  \item{outgraphs}{A boolean variable defining if the graphics on the marginal inclusion probabilities should be drawn (must not be used inside mclapply wrapper of runemjmcmc since otherwise errors can occur)}
  \item{advanced.param}{Currently omitted}
  \item{distrib_of_neighbourhoods}{A matrix defining probability distribution on 7 types of neighbourhoods within 4 possible local search strategies as well as within global moves}
  \item{distrib_of_proposals}{Probability distribution up to a constant of proportionality for addressing different local search strategies after large jumps or no large jumps (5th component)}
}
\details{See Hubin & Storvik (2016),Hubin, Storvik & Frommlet (2017), Hubin & Storvik (2017) details. The local optimization is performed through stepwise search within a neighborhood in the current gamma vector, allowing one component to be changed at a time.}
\value{A list containing
  \item{p.post}{A vector of posterior probabilities of the final vector of active covariates (features)}
  \item{m.post}{A vector of posterior probabilities of the models from the search space induced by the final vector of active covariates (features)}
  \item{s.mass}{Sum of marginal likelihoods times the priors from the explored part of the search space induced by the final vector of active covariates (features)}
}
\references{Hubin & Storvik (2016),Hubin, Storvik & Frommlet (2017), Hubin & Storvik (2017)}
\author{Aliaksandr Hubin}

\seealso{global objects statistics1 (if create.table==T) or hashStat (if create.table==F) contain all marginal likelihoods and two other model selection criteria as well as all of the beta coefficients for the models (if save.beta==T)}
\examples{
data(crimedata)


if(!("RCurl" %in% rownames(installed.packages()))) 
  install.packages("RCurl")
library(RCurl)
simx <- read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/Protein%20Activity%20Data/proteincen.txt"),sep = " ")
data.example <- as.data.frame(simx)
names(data.example)[89]="Y"

system.time({
  
  formula1 = as.formula(paste(colnames(data.example)[89],"~ 1 +",paste0(colnames(data.example)[-89],collapse = "+")))
  
  res = runemjmcmc(formula = formula1,data = data.example,recalc_margin = 2^10,estimator =estimate.bas.lm,estimator.args =  list(data = data.example,prior = 3, g = 96 ,n=96),save.beta = T,interact = T,relations = c("","sin","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=2,mutation_rate = 100, max.tree.size = 200000, Nvars.max = 95,p.allow.replace=0.9,p.allow.tree=0.5,p.nor=0.3,p.and = 0.7),n.models = 50000,unique = T,max.cpu = 10,max.cpu.glob = 10,create.table = F,create.hash = T,pseudo.paral = F,burn.in = 100,print.freq = 100,advanced.param = list(
    max.N.glob=as.integer(20),
    min.N.glob=as.integer(5),
    max.N=as.integer(3),                                                                                                    min.N=as.integer(1),
    printable = F))
})

}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line