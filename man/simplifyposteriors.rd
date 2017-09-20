\name{simplifyposteriors}
\alias{simplifyposteriors}
\title{A function that ads up posteriors for the same expression written in different character form in different parallel runs of the algorithm (mainly for Logic Regression and Deep Regression contexts)}
\usage{simplifyposteriors(X,posteriors,th=0.0001,thf=0.2, resp)}
\arguments{
\item{X}{a data.frame containing the data on the covariates}
\item{posteriors}{a data.frame with expressions in the first column and their posteriors in the second column from all of the runs}
\item{th}{initial filtering before summarization treshold}
\item{thf}{treshold for final filtering after summarization}
\item{resp}{the response to be addressed}
}
\value{
\item{res}{a data.frame with the summirized across runs expressions and their posteriors}
}
\seealso{runemjmcmc}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line