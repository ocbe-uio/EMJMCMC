\name{estimate.bas.glm}
\alias{estimate.bas.glm}
\title{Obtaining Bayesian estimators of interest from a GLM model}
\usage{estimate.bas.glm(formula, data, family, prior, logn)}
\arguments{
\item{formula}{a formula object for the model to be addressed}
\item{data}{a data frame object containing variables and observations corresponding to the formula used}
\item{family}{either poisson() or binomial(), that are currently adopted within this function}
\item{prior}{aic.prior(),bic.prior() or ic.prior() are allowed}
\item{logn}{log sample size}
}
\value{A list of
\item{mlik}{marginal likelihood of the model}
\item{waic}{AIC model selection criterion}
\item{dic}{BIC model selection criterion}
\item{summary.fixed$mean}{a vector of posterior modes of the parameters}
}
\seealso{BAS::bayesglm.fit}
\examples{

if(!("RCurl" %in% rownames(installed.packages()))) 
  install.packages("RCurl")
library(RCurl)

simx <- read.table(text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/US%20Data/simcen-x1.txt"),sep = ",")
simy <- read.table(text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/US%20Data/simcen-y1.txt"))
data.example <- cbind(simy,simx)
names(data.example)[1]="Y"
data.example$Y=as.integer(data.example$Y>mean(data.example$Y))
formula1 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-1],collapse = "+")))
estimate.bas.glm(formula = formula1, data = data.example,prior = aic.prior(), logn=47, family = binomial())
}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line