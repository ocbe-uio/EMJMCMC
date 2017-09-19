\name{estimate.glm}
\alias{estimate.glm}
\title{Obtaining Bayesian estimators of interest from a GLM model}
\usage{estimate.glm(formula, data, family, prior, n=1, g = 0)}
\arguments{
\item{formula}{a formula object for the model to be addressed}
\item{data}{a data frame object containing variables and observations corresponding to the formula used}
\item{family}{distribution family for the responces}
\item{prior}{integers 1,2 or 3 corresonding to AIC, BIC or Gelman's g-prior}
\item{n}{sample size}
\item{g}{g parameter of Gelman's g prior}
}
\value{a list of
\item{mlik}{marginal likelihood of the model}
\item{waic}{AIC model selection criterion}
\item{dic}{BIC model selection criterion}
\item{summary.fixed$mean}{a vector of posterior modes of the parameters}
}
\seealso{glm}
\examples{

library(RCurl)

simx <- read.table(text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/US%20Data/simcen-x1.txt"),sep = ",")
simy <- read.table(text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/US%20Data/simcen-y1.txt"))
data.example <- cbind(simy,simx)
names(data.example)[1]="Y"
formula1 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-1],collapse = "+")))
estimate.glm(formula = formula1, data = data.example,prior = 2, family = gaussian())
}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line