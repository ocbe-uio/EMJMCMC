#' @title Obtaining Bayesian estimators of interest from a GLM model
#' @param formula a formula object for the model to be addressed
#' @param data a data frame object containing variables and observations
#' corresponding to the formula used
#' @param family distribution family foe the responses
#' @param prior either "AIC" or "BIC"
#' @param logn log sample size
#' @return
#' \describe{
#'  \item{mlik}{marginal likelihood of the model}
#'  \item{waic}{AIC model selection criterion}
#'  \item{dic}{BIC model selection criterion}
#'  \item{summary.fixed$mean}{a vector of posterior modes of the parameters}
#' }
#' @seealso speedglm::speedglm.wfit
#' @example inst/examples/estimate.logic.lm_example.R
#' @keywords methods models
#' @export
estimate.speedglm <- function(formula, data, family, prior, logn) # weird behaviour, bad control of singularity
{
# use dic and aic as bic and aic correspondingly
X <- stats::model.matrix(object = formula,data = data)
out <- speedglm::speedglm.wfit(y = data[,1], X = X, intercept=FALSE, family=family,eigendec = TRUE, method = "Cholesky")
if(prior == "AIC")
  return(list(mlik = -out$aic ,waic = -(out$deviance + 2*out$rank) , dic =  -(out$RSS),summary.fixed =list(mean = out$coefficients)))
if(prior=="BIC")
  return(list(mlik = -out$RSS-logn*out$rank ,waic = -(out$deviance + 2*out$rank) , dic =  -(out$RSS+logn*out$rank),summary.fixed =list(mean = out$coefficients)))
}
