\name{simplify.formula}
\alias{simplify.formula}
\title{A function parsing the formula into the vectors of charactor arrays of responses and covariates}
\usage{simplify.formula(fmla,names)}
\arguments{
\item{fmla}{an R formula object}
\item{names}{all column names from the data.frame to be used with the formula}
}
\value{a list of
\item{fobserved}{a vector of charactor arrays corresponding to the observations}
\item{fparam}{a vector of charactor arrays corresponding to the covariates}
}
\seealso{formula, data.frame}
\examples{
X1<- as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = 0.3),dim = c(1000,50)))
Y1=-0.7+1*((1-X1$V1)*(X1$V4)) + 1*(X1$V8*X1$V11)+1*(X1$V5*X1$V9)
X1$Y1<-round(1.0/(1.0+exp(-Y1)))

formula1 = as.formula(paste(colnames(X1)[51],"~ 1 +",paste0(colnames(X1)[-c(51)],collapse = "+")))
names = colnames(X1)
simplify.formula(fmla=formula1,names=names)

}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line