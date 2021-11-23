estimate.bas.glm <- function(formula, data, family, prior, logn)
{

#only poisson and binomial families are currently adopted
X <- stats::model.matrix(object = formula,data = data)
out <- BAS::bayesglm.fit(x = X, y = data[,1], family=family,coefprior=prior)
# use dic and aic as bic and aic correspondinly
return(list(mlik = out$logmarglik,waic = -(out$deviance + 2*out$rank) , dic =  -(out$deviance + logn*out$rank),summary.fixed =list(mean = stats::coefficients(out))))

}
