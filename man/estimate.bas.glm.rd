\name{estimate.bas.glm}
\alias{estimate.bas.glm}
\title{Obtaining Bayesian estimators of interest from a GLM model}
\usage{estimate.bas.glm(formula, data, family, prior, logn)}
\arguments{
\item{formula}{a formula object for the model to be addressed}
\item{data}{a data frame object containing variables and observations corresponding to the formula used}
\item{family}{either poisson() or binomial(), that are currently adopted within this function}
\item{prior}{BAS::aic.prior(),bic.prior() or ic.prior() are allowed}
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

X4= as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))
Y4=rnorm(n = 1000,mean = 1+7*(X4$V4*X4$V17*X4$V30*X4$V10)+7*(((X4$V50*X4$V19*X4$V13*X4$V11)>0)) + 9*(X4$V37*X4$V20*X4$V12)+ 7*(X4$V1*X4$V27*X4$V3)
             +3.5*(X4$V9*X4$V2) + 6.6*(X4$V21*X4$V18) + 1.5*X4$V7 + 1.5*X4$V8,sd = 1)
X4$Y4=Y4
data.example = as.data.frame(X4)
data.example$Y4=as.integer(data.example$Y>mean(data.example$Y))
formula1 = as.formula(paste(colnames(X4)[51],"~ 1 +",paste0(colnames(X4)[-c(51)],collapse = "+")))

EMJMCMC:::estimate.bas.glm(formula = formula1, data = data.example,prior = BAS::aic.prior(), logn=47, family = binomial())
}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line
