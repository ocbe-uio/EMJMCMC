\name{runemjmcmc}
\alias{runemjmcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mode jumping MJMCMC or Genetically Modified Mode jumping MCMC for variable selection, Bayesian model averaging and feature engineering}
\description{The algorithm is an extended Metropolis-Hastings algorithm (or its Genetically modified version) mixing single site changes with occationally large jumps. The models are described through the gamma vector, a binary vector indicating which variables that are included in the model.}
\usage{
runemjmcmc(formula, data, secondary = vector(mode="character", length=0),
                     estimator,estimator.args = "list",n.models, unique = FALSE, save.beta= FALSE,  locstop.nd = FALSE,  latent="",max.cpu=4,max.cpu.glob=2,create.table= TRUE,  hash.length = 20, presearch= TRUE,  locstop =F ,pseudo.paral = FALSE, interact = FALSE, relations = c("","sin","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.1,0.1,0.1,0.1,0.1,0.1),gen.prob = c(1,10,5,1,1),pool.cross = 0.9,p.epsilon = 0.0001,p.add.default = 1,p.add = 0.5, del.sigma = 0.5,pool.cor.prob = FALSE,  interact.param=list(allow_offsprings=2,mutation_rate = 100,last.mutation=2000, max.tree.size = 10000, Nvars.max = 100, p.allow.replace = 0.7,p.allow.tree=0.1,p.nor=0.3,p.and = 0.7), prand = 0.01,keep.origin = TRUE,  sup.large.n = 5000, recalc_margin = 2^10, create.hash= FALSE, interact.order=1,burn.in=100, eps = 10^6, max.time = 120,max.it = 25000, print.freq = 100,outgraphs= FALSE, advanced.param=NULL, distrib_of_neighbourhoods,  distrib_of_proposals = c(76.91870,71.25264,87.68184,60.55921,15812.39852))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a typical formula for specifying a model with all potential covariates included}
  \item{data}{a data frame containing both covariates and response}
  \item{secondary}{a charactor vector of names other covariates excluded from those defined in formula (relevant for GMJMCMC only)}
  \item{estimator}{a function returning a list with marginal likelihood, waic, dic and coefficients of the addressed model. The list should be of a format: list(mlik = mlik,waic = waic , dic = dic,summary.fixed =list(mean = coefficients))}
  \item{estimator.args}{a list of arguments of estimator functions to be used (formula parameter has to be omitted, see the example)}
  \item{n.models}{maximal number of models to be estimated during the search}
  \item{unique}{defines whether n.models allows repetitions of the same models (unique=FALSE) or not (unique=TRUE)}
  \item{locstop.nd}{Defines whether local greedy optimizers stop at the first local optima found (locstop.nd=TRUE) or not (locstop.nd=FALSE)}
  \item{latent}{a latent random field to be adressed (to be specifucally used when estimator = INLA)}
  \item{create.table}{a boolean variable defining if a big.memory based hash table (only available for MJMCMC with no feature engineering, allows data sharing between CPUs) or the original R hash data structure (available for all algorithm, does not allow data sharing between CPUs) is used for storing of the results}
  \item{hash.length}{a paramter defining hash size for the big.memory based hash table as 2^hash.length (only relevant when create.table = T)}
  \item{pseudo.paral}{defines if lapply or mclapply is used for local vectorized computations within the chain (can only be TRUE if create.table=T)}
  \item{max.cpu}{maximal number of cpus in MJMCMC when within chain parallelization is allowed pseudo.paral = F}
  \item{max.cpu.glob}{maximal number of cpus in global moves in MJMCMC when within chain parallelization is allowed pseudo.paral = F}
  \item{presearch}{a boolean parameter defining if greedy forward and backward regression steps are used for initialization of initial approximations of marginal inclusion probabilities}
  \item{locstop}{a boolean parameter defining if the presearch is stopped at the first local exremum visited}
  \item{interact}{a boolean parameter defining if feature engineering is allowed in the search}
  \item{relations}{a vector of allowed modification functions (only relevant when feature engineering is enabled by means of interact = T)}
  \item{relations.prob}{probability distribution of addressing modifications defined in relations parameter (both vectors must be of the same length)}
  \item{gen.prob}{a vector of probabilities for different operators in GMJMCMC or RGMJMCMC in the deep regression context (hence only relevant if interact.param$allow_offsprings is either 3 or 4)}
  \item{pool.cross}{a parameter defining the probability of adressing covariates from the current pool of covariates in GMJMCMC (covariates from the set of filtered covariates can be addressed with probability 1-pool.cross) (only relevant when interact = TRUE)}
  \item{p.add}{a default marginal inclusion probability parameter to be changed during the search to the true value}
  \item{o.add.default}{a parameter defining sparcity after filtrations in gmjmcmc as initial marginal inclusion probabilities vector for paramters in the current pool}
  \item{p.epsilon}{a parameter to define minimal deviations from 0 and 1 probabilities when allowing adaptive MCMC based on marginal inclusion probabilities}
  \item{del.sigma}{a parameter decribing probability of deleting each of the function from the selected faeture in the reduction operator(only relevant for the deep regression models context)}
  \item{pool.cor.prob}{a boolean parameter indicating if inclusion of the filtered covariates during mutations are based on probabilities proportional to the absoulute values of correlations of these parameters and the observations (should not be addressed for multivariate observations, e.g. survival studies with Cox regression)}
  \item{interact.param}{a list of parameters for GMJMCMC, where allow_offsprings is 1 for logic regression context, 2 for the old version of GMJMCMC for deep regressions, 3 for the new version of GMJMCMC for deep regressions and 4 for the RGMJMCMC for the deep regressions; mutation_rate defines how often changes of the search space are allowed in terms of the number of MJMCMC iterations per search space; last.mutation defines the iteration after which changes of search space are no longer allowed; max.tree.size is a parameter defining maximal depth of features; Nvars.max is a parameter defining maximal number of covariates in the search space after the first filtration; p.allow.replace is a parameter defining the upepr bound on the probability allowing the replacement of correspondin features with marginal inclusion probabilities below it; p.allow.tree is a lower bound for the probability of not being filtered out after initializing steps of MJMCMC in GMJMCMC; p.nor is a parameter for not operator in the logic regression context (allow_offsprings==1); p.and = is the probability of & crossover in the logic regression context (allow_offsprings==1)}
  \item{prand}{probability of changes of components in randomization kernels of RGMJMCMC}
  \item{keep.origin}{a boolean parameter defining if the initially unfiltered covariates can leave the search space afterwards (TRUE) or not (FALSE)}
  \item{sup.large.n}{omitted currently}
  \item{recalc_margin}{a parameter defining how often marginal inclusion probabilities whould be recalculated}
  \item{create.hash}{a parameter defining if by default the results are stored in a hash table}
  \item{interact.order}{omitted currently}
  \item{burn.in}{number of burn-in steps for (R)(G)MJMCMC}
  \item{eps}{omitted, not to be changed}
  \item{max.time}{maximal time for the run of (R)(G)MJMCMC algorithm in minutes}
  \item{max.it}{maximal number of (R)(G)MJMCMC iterations}
  \item{print.freq}{printing frequency of the intermediate results}
  \item{outgraphs}{a boolean variable defining if the graphics on the marginal inclusion probabilities should be drawn (must not be used inside mclapply wrapper of runemjmcmc since otherwise errors can occur)}
  \item{advanced.param}{omitted currently}
  \item{distrib_of_neighbourhoods}{a matrix defining probability distribution on 7 types of neighbourhoods within 4 possible local search strategies as well as within global moves}
  \item{distrib_of_proposals}{probability distribution up to a constant of proportionality for addressing different local search strategies after large jumps or no large jumps (5th component)}
}
\details{See Hubin & Storvik (2016),Hubin, Storvik & Frommlet (2017), Hubin & Storvik (2017) details. The local optimization is performed through stepwise search within a neighborhood in the current gamma vector, allowing one component to be changed at a time.}
\value{a list containing
  \item{p.post}{a vector of posterior probabilities of the final vector of active covariates (features)}
  \item{m.post}{a vector of posterior probabilities of the models from the search space induced by the final vector of active covariates (features)}
  \item{s.mass}{sum of marginal likelihoods times the priors from the explored part of the search space induced by the final vector of active covariates (features)}
}
\references{Hubin & Storvik (2016),Hubin, Storvik & Frommlet (2017), Hubin & Storvik (2017)}
\author{Aliaksandr Hubin}

\seealso{global objects statistics1 (if create.table==T) or hashStat (if create.table==F) contain all marginal likelihoods and two other model selection criteria as well as all of the beta coefficients for the models (if save.beta==T)}
\examples{

    X4= as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))
    Y4=rnorm(n = 1000,mean = 1+7*(X4$V4*X4$V17*X4$V30*X4$V10)+7*(((X4$V50*X4$V19*X4$V13*X4$V11)>0)) + 9*(X4$V37*X4$V20*X4$V12)+ 7*(X4$V1*X4$V27*X4$V3)
             +3.5*(X4$V9*X4$V2) + 6.6*(X4$V21*X4$V18) + 1.5*X4$V7 + 1.5*X4$V8,sd = 1)
    X4$Y4=Y4
    data.example = as.data.frame(X4)

    #specify the initial formula
    formula1 = as.formula(paste(colnames(X4)[51],"~ 1 +",paste0(colnames(X4)[-c(51)],collapse = "+")))


    #specify tuning parameters of the algorithm for exploring DBRM of interest
    #notice that allow_offsprings=3 corresponds to the GMJMCMC runs and
    #allow_offsprings=4 -to the RGMJMCMC runs
     res = EMJMCMC:::runemjmcmc(formula = formula1,outgraphs= FALSE, data = X4,estimator = estimate.gamma.cpen,estimator.args =  list(data = data.example),recalc_margin = 249, save.beta = FALSE, interact = TRUE, relations=c("cos","sigmoid","tanh","atan","sin","erf"),relations.prob =c(0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=4,mutation_rate = 250,last.mutation = 15000, max.tree.size = 4, Nvars.max =40,p.allow.replace=0.7,p.allow.tree=0.2,p.nor=0,p.and = 0.9),n.models = 20000,unique = TRUE, max.cpu = 4,max.cpu.glob = 4,create.table = FALSE, create.hash = TRUE, pseudo.paral = TRUE, burn.in = 50,print.freq = 1000,advanced.param = list(
      max.N.glob=as.integer(10),
      min.N.glob=as.integer(5),
      max.N=as.integer(3),
      min.N=as.integer(1),
      printable = F))

}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line
