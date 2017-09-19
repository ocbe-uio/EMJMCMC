\name{estimate.bigm}
\alias{estimate.bigm}
\title{Obtaining Bayesian estimators of interest from a GLM model}
\usage{estimate.bigm(formula, data, family, prior, maxit = 2,chunksize = 1000000)}
\arguments{
\item{formula}{a formula object for the model to be addressed}
\item{data}{a data frame object containing variables and observations corresponding to the formula used}
\item{family}{distribution family foe the responces}
\item{prior}{either "AIC" or "BIC"}
\item{maxit}{maximum number of Fisher scoring iterations}
\item{chunksize}{size of chunks for processng the data frame}
}
\value{a list of
\item{mlik}{marginal likelihood of the model}
\item{waic}{AIC model selection criterion}
\item{dic}{BIC model selection criterion}
\item{summary.fixed$mean}{a vector of posterior modes of the parameters}
\item{n}{sample size}
}
\seealso{biglm::bigglm}
\examples{

library(RCurl)

simx <- read.table(text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/US%20Data/simcen-x1.txt"),sep = ",")
simy <- read.table(text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/US%20Data/simcen-y1.txt"))
data.example <- cbind(simy,simx)
names(data.example)[1]="Y"
formula1 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-1],collapse = "+")))
estimate.bigm(formula = formula1, data = data.example,n=47,prior = "BIC", maxit = 20,chunksize = 1000000, family = gaussian())
}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line