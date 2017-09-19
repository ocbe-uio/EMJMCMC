\name{estimate.logic.glm}
\alias{estimate.logic.glm}
\title{Obtaining Bayesian estimators of interest from a GLM model in a logic regression context}
\usage{estimate.logic.glm(formula, data, family, n, m, r = 1)}
\arguments{
\item{formula}{a formula object for the model to be addressed}
\item{data}{a data frame object containing variables and observations corresponding to the formula used}
\item{family}{either poisson() or binomial(), that are currently adopted within this function}
\item{n}{sample size}
\item{m}{total number of input binary leaves}
\item{r}{omitted}
}
\value{a list of
\item{mlik}{marginal likelihood of the model}
\item{waic}{AIC model selection criterion}
\item{dic}{BIC model selection criterion}
\item{summary.fixed$mean}{a vector of posterior modes of the parameters}
}
\seealso{BAS::bayesglm.fit, estimate.logic.lm}
\examples{

X1<- as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = 0.3),dim = c(1000,50)))
Y1=-0.7+1*((1-X1$V1)*(X1$V4)) + 1*(X1$V8*X1$V11)+1*(X1$V5*X1$V9)
X1$Y1<-round(1.0/(1.0+exp(-Y1)))

formula1 = as.formula(paste(colnames(X1)[51],"~ 1 +",paste0(colnames(X1)[-c(51)],collapse = "+")))

estimate.logic.glm(formula = formula1, data = X1,family = binomial(),n = 1000, m = 50)

}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line