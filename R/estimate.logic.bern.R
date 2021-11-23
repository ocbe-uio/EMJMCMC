estimate.logic.bern = function(formula = NULL, data, family = stats::binomial(), n=1000, m=50, r = 1,k.max=21)
{
#define the function estimating parameters of a given Bernoulli logic regression with Jeffrey's prior
  if(length(formula)==0)
    return(list(mlik =  -10000 + stats::rnorm(1,0,1),waic =10000 , dic =  10000,summary.fixed =list(mean = 1)))

  if(is.na(formula))
    return(list(mlik =  -10000 + stats::rnorm(1,0,1),waic =10000 , dic =  10000,summary.fixed =list(mean = 1)))

  out = stats::lm(formula = formula,data = data, family=family)
  p = out$rank
  if(p>k.max)
  {
    return(list(mlik = -10000,waic = 10000 , dic =  10000,summary.fixed =list(mean = 0)))
  }
  fmla.proc=as.character(formula)[2:3]
  fobserved = fmla.proc[1]
  fmla.proc[2]=stringi::stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]=stringi::stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam =stringi::stri_split_fixed(str = fmla.proc[2],pattern = "+",omit_empty = F)[[1]]
  sj=(stringi::stri_count_fixed(str = fparam, pattern = "&"))
  sj=sj+(stringi::stri_count_fixed(str = fparam, pattern = "|"))
  sj=sj+1
  Jprior = sum(log(factorial(sj)/((m^sj)*2^(2*sj-2))))
  mlik = (-(out$deviance + log(n)*(out$rank)) + 2*(Jprior))/2+n
  if(mlik==-Inf)
    mlik = -10000
  return(list(mlik = mlik,waic = -(out$deviance + 2*out$rank) , dic =  -(out$deviance + log(n)*out$rank),summary.fixed =list(mean = stats::coefficients(out))))
}
