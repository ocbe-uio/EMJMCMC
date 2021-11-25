#' @title Obtaining Bayesian estimators of interest from a GLM model
#' @param formula a formula object for the model to be addressed
#' @param data a data frame object containing variables and observations
#' corresponding to the formula used
#' @param family either poisson() or binomial(), that are currently adopted
#' within this function
#' @param prior BAS::aic.prior(), bic.prior() or ic.prior() are allowed
#' @param logn log sample size
#' @return  A list of
#' \describe{
#'  \item{mlik}{marginal likelihood of the model}
#'  \item{waic}{AIC model selection criterion}
#'  \item{dic}{BIC model selection criterion}
#'  \item{summary.fixed$mean}{a vector of posterior modes of the parameters}
#' }
#' @seealso BAS::bayesglm.fit
#' @example /inst/examples/estimate.bas.glm_example.R
#' @keywords methods models
#' @export
estimate.bas.glm <- function(formula, data, family, prior, logn)
{

#only poisson and binomial families are currently adopted
X <- stats::model.matrix(object = formula,data = data)
out <- BAS::bayesglm.fit(x = X, y = data[,1], family=family,coefprior=prior)
# use dic and aic as bic and aic correspondinly
return(list(mlik = out$logmarglik,waic = -(out$deviance + 2*out$rank) , dic =  -(out$deviance + logn*out$rank),summary.fixed =list(mean = stats::coefficients(out))))

}
