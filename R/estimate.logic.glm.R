#' @title Obtaining Bayesian estimators of interest from a GLM model in a
#' logic regression context
#' @param formula a formula object for the model to be addressed
#' @param data a data frame object containing variables and observations corresponding to the formula used
#' @param family either poisson() or binomial(), that are currently adopted within this function
#' @param n sample size
#' @param m total number of input binary leaves
#' @param r omitted
#' @return a list of
#' \describe{
#'  \item{mlik}{marginal likelihood of the model}
#'  \item{waic}{AIC model selection criterion}
#'  \item{dic}{BIC model selection criterion}
#'  \item{summary.fixed$mean}{a vector of posterior modes of the parameters}
#' }
#' @seealso BAS::bayesglm.fit estimate.logic.lm
#' @example inst/examples/estimate.logic.glm_example.R
#' @keywords methods models
#' @export
estimate.logic.glm <- function(formula, data, family, n, m, r = 1)
{
X <- stats::model.matrix(object = formula,data = data)
out <- BAS::bayesglm.fit(x = X, y = data[,51], family=family,coefprior=BAS::aic.prior())
p <- out$rank
fmla.proc<-as.character(formula)[2:3]
fobserved <- fmla.proc[1]
fmla.proc[2]<-stringi::stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
fmla.proc[2]<-stringi::stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
fparam <-stringi::stri_split_fixed(str = fmla.proc[2],pattern = "+",omit_empty = F)[[1]]
sj<-(stringi::stri_count_fixed(str = fparam, pattern = "&"))
sj<-sj+(stringi::stri_count_fixed(str = fparam, pattern = "|"))
sj<-sj+1
Jprior <- sum(log(truncfactorial(sj)/((m^sj)*2^(2*sj-2))))
mlik = (-(out$deviance + log(n)*(out$rank)) + 2*(Jprior))/2+n
if(mlik==-Inf)
  mlik = -10000
return(list(mlik = mlik,waic = -(out$deviance + 2*out$rank) , dic =  -(out$deviance + log(n)*out$rank),summary.fixed =list(mean = stats::coefficients(out))))
}
