\name{estimate.elnet }
\alias{estimate.elnet }
\title{A test function to work with elastic networks in future, be omitted so far}
\usage{estimate.elnet(formula,response, data, family,alpha)}
\arguments{
\item{formula}{a formula object for the model to be addressed}
\item{data}{a data frame object containing variables and observations corresponding to the formula used}
\item{response}{response in a formula}
\item{family}{distribution of the response family object}
\item{alpha}{regularization parameter in [0,1]}
}
\value{a list of
\item{mlik}{marginal likelihood of the model}
\item{waic}{AIC model selection criterion}
\item{dic}{BIC model selection criterion}
\item{summary.fixed$mean}{a vector of posterior modes of the parameters}
}
\seealso{glmnet::glmnet}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line