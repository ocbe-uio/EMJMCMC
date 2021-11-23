estimate.logic.lm.tCCH = function(formula = NULL, data, n=1000, m=50, r = 1, p.a = 1, p.b = 2, p.r = 1.5, p.s = 0, p.v=-1, p.k = 1,k.max=21)
{
#define the function estimating parameters of a given Gaussian logic regression with robust g prior
  if(is.null(formula))
    return(list(mlik =  -10000,waic =10000 , dic =  10000,summary.fixed =list(mean = 1)))
  fmla.proc=as.character(formula)[2:3]
  fobserved = fmla.proc[1]
  if(fmla.proc[2]=="-1")
    return(list(mlik =  -10000,waic =10000 , dic =  10000,summary.fixed =list(mean = 1)))
  out = stats::lm(formula = formula,data = data)
  p = out$rank
  if(p>k.max)
  {
    return(list(mlik = -10000,waic = 10000 , dic =  10000,summary.fixed =list(mean = 0)))
  }
  fmla.proc[2]=stringi::stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]=stringi::stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam =stringi::stri_split_fixed(str = fmla.proc[2],pattern = "+",omit_empty = F)[[1]]
  sj=(stringi::stri_count_fixed(str = fparam, pattern = "&"))
  sj=sj+(stringi::stri_count_fixed(str = fparam, pattern = "|"))
  sj=sj+1
  Jprior = prod(factorial(sj)/((m^sj)*2^(2*sj-2)))
  p.v = (n+1)/(p+1)
  R.2 = summary(out)$r.squared

  mlik = (-0.5*p*log(p.v) -0.5*(n-1)*log(1-(1-1/p.v)*R.2) + log(beta((p.a+p)/2,p.b/2)) - log(beta(p.a/2,p.b/2)) + log(BAS::phi1(p.b/2,(n-1)/2,(p.a+p.b+p)/2,p.s/2/p.v,R.2/(p.v-(p.v-1)*R.2))) - BAS::hypergeometric1F1(p.b/2,(p.a+p.b)/2,p.s/2/p.v,log = T)+log(Jprior) + p*log(r)+n)
  if(mlik==-Inf||is.na(mlik)||is.nan(mlik))
    mlik = -10000
  return(list(mlik = mlik,waic = stats::AIC(out)-n , dic =  stats::BIC(out)-n,summary.fixed =list(mean = stats::coef(out))))
}
