\name{parallelize}
\alias{parallelize}
\title{An example of user defined parallelization (cluster based) function for within an MJMCMC chain calculations (mclapply or lapply are used by default depending on specification and OS).}
\usage{parallelize(X,FUN)}
\arguments{
\item{X}{a vector (atomic or list) or an expressions vector. Other objects (including classed objects) will be coerced by as.list}
\item{FUN}{the function to be applied to each element of X or v, or in parallel to X}
}
\value{
\item{parallelize(X,FUN)}{a list of the same length as X and named by X}
}
\details{Only allowed when working with big.memory based hash table within MJMCMC (see runemjmcmc for more details)}
\seealso{parLapply, clusterMap, mclapply, lapply, etc.}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line