estimate.bigm <- function(formula, data, family, prior,n, maxit = 2,chunksize = 1000000) # nice behaviour
{

out <- stats::lm(data = data, family=family,formula = formula, sandwich = F,maxit = maxit, chunksize = chunksize)
if(prior == "AIC")
  return(list(mlik = -stats::AIC(out,k = 2) ,waic = stats::AIC(out,k = 2) , dic =  stats::AIC(out,k = n),summary.fixed =list(mean = stats::coef(out))))
if(prior=="BIC")
  return(list(mlik = -stats::AIC(out,k = n) ,waic = stats::AIC(out, k = 2) , dic =  stats::AIC(out,k = n),summary.fixed =list(mean = stats::coef(out))))
}
