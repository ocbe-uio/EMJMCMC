\name{estimate.inla}
\alias{estimate.inla}
\title{Obtaining Bayesian estimators of interest from a GLM model}
\usage{estimate.inla(formula, args)}
\arguments{
\item{formula}{a formula object for the model to be addressed}
\item{args}{inla arguments}
}
\value{a list of
\item{mlik}{marginal likelihood of the model}
\item{waic}{WAIC model selection criterion}
\item{dic}{DIC model selection criterion}
\item{summary.fixed$mean}{a vector of posterior modes of the parameters}
}
\seealso{INLA::inla}
\examples{

library(RCurl)

simx <- read.table(text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/US%20Data/simcen-x1.txt"),sep = ",")
simy <- read.table(text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/US%20Data/simcen-y1.txt"))
data.example <- cbind(simy,simx)
names(data.example)[1]="Y"
formula1 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-1],collapse = "+")))
estimate.inla(formula = formula1, args = list(data = data.example,control.compute=list(dic = T, waic=T)))
}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line