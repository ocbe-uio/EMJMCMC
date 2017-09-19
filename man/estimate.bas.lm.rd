\name{estimate.bas.lm}
\alias{estimate.bas.lm}
\title{Obtaining Bayesian estimators of interest from a LM model}
\usage{estimate.bas.lm(formula, data, prior, n, g = 0)}
\arguments{
\item{formula}{a formula object for the model to be addressed}
\item{data}{a data frame object containing variables and observations corresponding to the formula used}
\item{prior}{integers 1, 2 or 3 are allowed corresponding to AIC, BIC or Gelman's g-prior}
\item{n}{sample size}
}
\value{a list of
\item{mlik}{marginal likelihood of the model}
\item{waic}{AIC model selection criterion}
\item{dic}{BIC model selection criterion}
\item{summary.fixed$mean}{a vector of posterior modes of the parameters}
}
\seealso{BAS::bayesglm.fit}
\examples{


library(RCurl)

simx <- read.table(text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/US%20Data/simcen-x1.txt"),sep = ",")
simy <- read.table(text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/US%20Data/simcen-y1.txt"))
data.example <- cbind(simy,simx)
names(data.example)[1]="Y"
formula1 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-1],collapse = "+")))
estimate.bas.lm(formula = formula1, data = data.example,prior = 2, n=47)

}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line