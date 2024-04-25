#' @title Estimate marginal log posterior of a single BGNLM model 
#' @importFrom stringi stri_replace_all stri_split_fixed stri_count_fixed
#' @param formula formula
#' @param data dataset
#' @param r prior inclusion penalty parameter
#' @param logn logn
#' @param relat a set of nonlinear transformations in the class of BGNLMs of interest
#' @return  A list of
#' \describe{
#'  \item{mlik}{marginal likelihood of the model}
#'  \item{waic}{AIC model selection criterion}
#'  \item{dic}{BIC model selection criterion}
#'  \item{summary.fixed$mean}{a vector of posterior modes of the parameters}
#' }
#' @export
estimate.gamma.cpen <- function(formula, data, r = 1.0 / 1000.0, logn = log(1000.0), relat = c("cos", "sigmoid", "tanh", "atan", "sin", "erf")) {
  fparam <- NULL
  fmla.proc <- as.character(formula)[2:3]
  fobserved <- fmla.proc[1]
  fmla.proc[2] <- stringi::stri_replace_all(str = fmla.proc[2], fixed = " ", replacement = "")
  fmla.proc[2] <- stringi::stri_replace_all(str = fmla.proc[2], fixed = "\n", replacement = "")
  fparam <- stringi::stri_split_fixed(str = fmla.proc[2], pattern = "+I", omit_empty = FALSE)[[1]]
  sj <- (stringi::stri_count_fixed(str = fparam, pattern = "*"))
  sj <- sj + (stringi::stri_count_fixed(str = fparam, pattern = "+"))
  for (rel in relat) {
    sj <- sj + (stringi::stri_count_fixed(str = fparam, pattern = rel))
  }
  # sj<-sj+1
  tryCatch(utils::capture.output(
    {
      out <- stats::glm(formula = formula, data = data, family = stats::gaussian)
      # 1 for aic, 2 bic prior, else g.prior

      mlik <- (-(stats::BIC(out) - 2 * log(r) * sum(sj)) + 1000) / 2
      waic <- (out$deviance + 2 * out$rank) + 10000
      dic <- (out$deviance + logn * out$rank) + 10000
      summary.fixed <- list(mean = stats::coefficients(out))
    },
    error = function(err) {
      print(err)
      mlik <- -10000
      waic <- 10000
      dic <- 10000
      summary.fixed <- list(mean = array(0, dim = length(fparam)))
    }
  ))
  return(list(mlik = mlik, waic = waic, dic = dic, summary.fixed = summary.fixed))
}

#' @title Estimate marginal log posterior of a single BGNLM model with alternative defaults
#' @importFrom stringi stri_replace_all stri_split_fixed stri_count_fixed
#' @param formula formula
#' @param data dataset
#' @param r prior inclusion penalty parameter
#' @param logn logn
#' @param relat a set of nonlinear transformations in the class of BGNLMs of interest
#' @return  A list of
#' \describe{
#'  \item{mlik}{marginal likelihood of the model}
#'  \item{waic}{AIC model selection criterion}
#'  \item{dic}{BIC model selection criterion}
#'  \item{summary.fixed$mean}{a vector of posterior modes of the parameters}
#' }
#' @export
estimate.gamma.cpen_2 = function(formula, data,r = 1.0/223.0,logn=log(223.0),relat=c("to23","expi","logi","to35","sini","troot","sigmoid"))
{
  fparam=NULL
  fmla.proc=as.character(formula)[2:3]
  fobserved = fmla.proc[1]
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam =stri_split_fixed(str = fmla.proc[2],pattern = "+I",omit_empty = FALSE)[[1]]
  sj=(stri_count_fixed(str = fparam, pattern = "*"))
  sj=sj+(stri_count_fixed(str = fparam, pattern = "+"))
  for(rel in relat)
    sj=sj+(stri_count_fixed(str = fparam, pattern = rel))
  sj=sj+1
  tryCatch(capture.output({
    out = glm(formula = formula,data = data, family = gaussian)
    mlik = (-(out$deviance -2*log(r)*sum(sj)))/2
    waic = -(out$deviance + 2*out$rank)
    dic =  -(out$deviance + logn*out$rank)
    summary.fixed =list(mean = coefficients(out))

  }, error = function(err) {
    print(err)
    mlik = -10000
    waic = -10000
    dic =  -10000
    summary.fixed =list(mean = array(0,dim=length(fparam)))
  }))
  return(list(mlik = mlik,waic = waic , dic = dic,summary.fixed =summary.fixed))

}

contrelu <- function(x) log(1 + exp(x))
cosi <- function(x) cos(x / 180 * pi)
expi <- function(x) exp(-abs(x))
gauss <- function(x) exp(-x * x)
gfifth <- function(x) (abs(x)^(1 / 5))
gfquar <- function(x) as.integer(x < quantile(x, probs = 0.25))
glquar <- function(x) as.integer(x > quantile(x, probs = 0.75))
gmean <- function(x) as.integer(x > mean(x))
gmedi <- function(x) as.integer(x > median(x))
gone <- function(x) as.integer(x > 0)
grelu <- function(x) (x * (x > 0))
gthird <- function(x) (abs(x)^(1 / 3))
logi <- function(x) log(abs(x + 0.1))
logi2 <- function(x) log(abs(x) + 1)
sini <- function(x) sin(x / 180 * pi)
to23 <- function(x) abs(x)^(2.3)
to25 <- function(x) abs(x)^(2.5)
to35 <- function(x) abs(x)^(3.5)
troot <- function(x) abs(x)^(1 / 3)
