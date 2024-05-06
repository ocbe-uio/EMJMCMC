estimate.logic.lm.jef = function(formula= NULL, data, n, m, r = 1,k.max=21)
{
#define the function estimating parameters of a given Gaussian logic regression with Jeffrey's prior
  if(is.null(formula))
    return(list(mlik =  -10000,waic =10000 , dic =  10000,summary.fixed =list(mean = 1)))
  out = stats::lm(formula = formula,data = data)
  p = out$rank
  if(p>k.max)
  {
    return(list(mlik = -10000,waic = 10000 , dic =  10000,summary.fixed =list(mean = 0)))
  }
  fmla.proc=as.character(formula)[2:3]
  fobserved = fmla.proc[1]
  fmla.proc[2]=stringi::stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]=stringi::stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam =stringi::stri_split_fixed(str = fmla.proc[2],pattern = "+",omit_empty = FALSE)[[1]]
  sj=(stringi::stri_count_fixed(str = fparam, pattern = "&"))
  sj=sj+(stringi::stri_count_fixed(str = fparam, pattern = "|"))
  sj=sj+1
  Jprior = prod(truncfactorial(sj)/((m^sj)*2^(2*sj-2)))
  mlik = (-stats::BIC(out)+2*log(Jprior) + 2*p*log(r)+n)/2
  if(mlik==-Inf)
    mlik = -10000
  return(list(mlik = mlik,waic = stats::AIC(out)-n , dic =  stats::BIC(out)-n,summary.fixed =list(mean = stats::coef(out))))
}
