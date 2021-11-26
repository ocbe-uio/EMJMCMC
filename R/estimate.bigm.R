#' @title Obtaining Bayesian estimators of interest from a GLM model
#' @param formula a formula object for the model to be addressed
#' @param data a data frame object containing variables and observations corresponding to the formula used
#' @param family distribution family foe the responces
#' @param prior either "AIC" or "BIC"
#' @param n sample size
#' @param maxit maximum number of Fisher scoring iterations
#' @param chunksize size of chunks for processng the data frame
#' @return a list of
#' \describe{
#'  \item{mlik}{marginal likelihood of the model}
#'  \item{waic}{AIC model selection criterion}
#'  \item{dic}{BIC model selection criterion}
#'  \item{summary.fixed$mean}{a vector of posterior modes of the parameters}
#'  \item{n}{sample size}
#' }
#' @seealso biglm::bigglm
#' @example /inst/examples/estimate.bigm_example.R
#' @keywords methods models
#' @export
estimate.bigm <- function(formula, data, family, prior, n, maxit = 2, chunksize = 1000000) # nice behaviour
{
  out <- biglm::bigglm(
    data = data, family = family, formula = formula, sandwich = F,
    maxit = maxit, chunksize = chunksize
  )
  if (prior == "AIC") {
    penalty <- 2
  }
  if (prior == "BIC") {
    penalty <- n
  }
  return(
    list(
      mlik = -stats::AIC(out, k = penalty),
      waic = stats::AIC(out, k = 2),
      dic = stats::AIC(out, k = n),
      summary.fixed = list(mean = stats::coef(out))
    )
  )
}
