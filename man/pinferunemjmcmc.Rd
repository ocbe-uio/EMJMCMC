\name{pinferunemjmcmc}
\alias{pinferunemjmcmc}
\title{A wrapper for running the GLMM, BLR, or DBRM based inference and predictions in an expert but rather easy to use way}
\usage{pinferunemjmcmc(n.cores = 4, mcgmj = mcgmjpar, report.level =  0.5, simplify = F, num.mod.best = 1000, predict = F,  test.data = 1, link.function = function(z)z, runemjmcmc.params)}
\arguments{
\item{runemjmcmc.params}{a vector of parameters of runemjmcmc function, see the help of runemjmcmc for details}
\item{mcgmj}{an mclapply like function for perfroming for perfroming parallel computing, do not change the default unless you are using Windows}
\item{simplify}{a logical value specifying in simplification of the features is to be done after the search is completed}
\item{predict}{a logical value specifying if predictions should be done by the run of pinferunemjmcmc}
\item{test.data}{covariates data.frame to be used for predictions}
\itme{link.function}{the link functions to be used to make predictions}
\item{report.level}{a numeric value in (0,1) specifying the treshold for detections based on the marginal inclusion probabilities}
\item{n.cores}{the maximal number of cores (and (R)(G)MJMCMC threads) to be addressed in the analysis}
\item{num.mod.best}{the number of the best models in the thread to calculate marginal inclusion probabilities}
}
\value{a list of
\item{feat.stat}{detected features or logical expressions and their marginal inclusion probabilities}
\item{predictions}{predicted values if they are required, NULL otherwise}
\item{allposteriors}{all visited by (R)(G)MJMCMC features and logical expressions and their marginal inclusion probabilities}
\item{threads.stats}{a vector of detailed outputs of individual n.cores threads of (R)(G)MJMCMC run}
}
\seealso{EMJMCMC::runemjmcmc, EMJMCMC::LogrRegr, EMJMCMC::DeepRegr, EMJMCMC::LinRegr}
\examples{

#inference

#prediction
pinferunemjmcmc(n.cores =30, report.level =  0.8 , num.mod.best = NM,simplify = T, predict = T,test.data = as.data.frame(test),link.function = g, runemjmcmc.params =list(formula = formula1,data = data.example,gen.prob = c(1,1,1,1,0),estimator =estimate.bas.glm.cpen,estimator.args =  list(data = data.example,prior = aic.prior(),family = binomial(),yid=31, logn = log(143),r=exp(-0.5)),recalc_margin = 95, save.beta = T,interact = T,relations = c("gauss","tanh","atan","sin"),relations.prob =c(0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=4,mutation_rate = 100,last.mutation=1000, max.tree.size = 6, Nvars.max = 20,p.allow.replace=0.5,p.allow.tree=0.4,p.nor=0.3,p.and = 0.9),n.models = 7000,unique =T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(max.N.glob=as.integer(10), min.N.glob=as.integer(5), max.N=as.integer(3), min.N=as.integer(1), printable = F)))

}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line
