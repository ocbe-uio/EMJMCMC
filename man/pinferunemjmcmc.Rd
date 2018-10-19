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
\item{feat.stat}{detected logical expressions and their marginal inclusion probabilities}
\item{predictions}{NULL currently, since LogrRegr function is not designed for predictions at the moment, which is still possible in its expert mother function pinferunemjmcmc}
\item{allposteriors}{all visited by GMJMCMC logical expressions and their marginal inclusion probabilities}
\item{threads.stats}{a vector of detailed outputs of individual ncores threads of GMJMCMC run}
}
\seealso{EMJMCMC::runemjmcmc, EMJMCMC::LogrRegr, EMJMCMC::DBRM}
\examples{

#inference

#prediction

}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line
