#' @title Obtaining Bayesian estimators of interest from a GLM model
#' @param formula a formula object for the model to be addressed
#' @param args inla arguments
#' @return
#' \describe{
#'  \item{mlik}{marginal likelihood of the model}
#'  \item{waic}{AIC model selection criterion}
#'  \item{dic}{BIC model selection criterion}
#'  \item{summary.fixed$mean}{a vector of posterior modes of the parameters}
#' }
#' @seealso INLA::inla
#' @example inst/examples/estimate.inla_example.R
#' @keywords methods models
#' @export
estimate.inla <- function(formula, args)
{

out<-NULL
tryCatch(utils::capture.output({
out <- do.call(INLA::inla, c(args,formula = formula))}))
if(is.null(out))
  return(list(mlik = -10000,waic =  10000 , dic = 10000, summary.fixed =list(mean = NULL)))
# use dic and aic as bic and aic correspondinly
coef<-out$summary.fixed$mode
coef[1]<-coef[1]+out$summary.hyperpar$mode[1]
return(list(mlik = out$mlik[1],waic =  out$waic[1] , dic = out$dic[1], summary.fixed =list(mean = coef)))

}
