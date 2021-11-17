\name{factorial}
\alias{factorial}
\title{truncated factorial to avoid stack overflow for huge values}
\usage{factorial(x)}
\arguments{
\item{x}{a non-negative integer number}
}
\value{
\item{factorial(x)}{truncated facctorial as min(x!,171!)}
}
\examples{EMJMCMC:::factorial(10)}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line