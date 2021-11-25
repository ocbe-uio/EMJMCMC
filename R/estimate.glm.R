#' @title Obtaining Bayesian estimators of interest from a GLM model
#' @param formula a formula object for the model to be addressed
#' @param data a data frame object containing variables and observations
#' corresponding to the formula used
#' @param family distribution family for the responces
#' @param prior integers 1,2 or 3 corresonding to AIC, BIC or Gelman's g-prior
#' @param n sample size
#' @param g g parameter of Gelman's g prior
#' @return a list of
#' \describe{
#'  \item{mlik}{marginal likelihood of the model}
#'  \item{waic}{AIC model selection criterion}
#'  \item{dic}{BIC model selection criterion}
#'  \item{summary.fixed$mean}{a vector of posterior modes of the parameters}
#' }
#' @seealso glm
#' @example inst/examples/estimate.glm_example.R
#' @keywords methods models
#' @export
estimate.glm <- function(formula, data, family, prior, n=1, g = 0)
{

out <- stats::lm(formula = formula, family = family, data = data)
# 1 for aic, 2 bic prior, else g.prior

p <- out$rank
if(prior == 1)
{

  logmarglik <- -stats::AIC(out)
}
else if(prior ==2)
{
  logmarglik <- -stats::BIC(out)
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
