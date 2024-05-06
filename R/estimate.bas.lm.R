#' @title Obtaining Bayesian estimators of interest from a LM model
#' @param formula a formula object for the model to be addressed
#' @param data a data frame object containing variables and observations corresponding to the formula used
#' @param prior integers 1, 2 or 3 are allowed corresponding to AIC, BIC or Zellner's g-prior
#' @param n sample size
#' @param g g
#' @return a list of
#' \describe{
#' \item{mlik}{marginal likelihood of the model}
#' \item{waic}{AIC model selection criterion}
#' \item{dic}{BIC model selection criterion}
#' \item{summary.fixed$mean}{a vector of posterior modes of the parameters}
#' }
#' @seealso BAS::bayesglm.fit
#' @example /inst/examples/estimate.bas.lm_example.R
#' @keywords methods models
#' @export
estimate.bas.lm <- function(formula, data, prior, n, g = 0)
{

out <- stats::lm(formula = formula,data = data)
# 1 for aic, 2 bic prior, else g.prior

p <- out$rank
if(prior == 1)
{
  ss<-sum(out$residuals^2)
  logmarglik <- -0.5*(log(ss)+2*p)
}
else if(prior ==2)
{
  ss<-sum(out$residuals^2)
  logmarglik <- -0.5*(log(ss)+log(n)*p)
}
else
{
  Rsquare <- summary(out)$r.squared
  #logmarglik =  .5*(log(1.0 + g) * (n - p -1)  - log(1.0 + g * (1.0 - Rsquare)) * (n - 1))*(p!=1)
  logmarglik =  .5*(log(1.0 + g) * (n - p)  - log(1.0 + g * (1.0 - Rsquare)) * (n - 1))*(p!=1)
}

# use dic and aic as bic and aic correspondinly
return(list(mlik = logmarglik,waic = stats::AIC(out) , dic =  stats::BIC(out),summary.fixed =list(mean = stats::coef(out))))

}
