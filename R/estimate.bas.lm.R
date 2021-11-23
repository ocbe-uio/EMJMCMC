estimate.bas.lm <- function(formula, data, prior, n, g = 0)
{

out <- stats::lm(formula = formula,data = data)
# 1 for aic, 2 bic prior, else g.prior

p <- out$rank
if(prior == 1)
{
  ss<-sum(out$residuals^2)
  logmarglik <- -0.5*(log(ss)+2*p)
}
else if(prior ==2)
{
  ss<-sum(out$residuals^2)
  logmarglik <- -0.5*(log(ss)+log(n)*p)
}
else
{
  Rsquare <- summary(out)$r.squared
  #logmarglik =  .5*(log(1.0 + g) * (n - p -1)  - log(1.0 + g * (1.0 - Rsquare)) * (n - 1))*(p!=1)
  logmarglik =  .5*(log(1.0 + g) * (n - p)  - log(1.0 + g * (1.0 - Rsquare)) * (n - 1))*(p!=1)
}

# use dic and aic as bic and aic correspondinly
return(list(mlik = logmarglik,waic = stats::AIC(out) , dic =  stats::BIC(out),summary.fixed =list(mean = stats::coef(out))))

}
