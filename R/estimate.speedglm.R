estimate.speedglm <- function(formula, data, family, prior, logn) # weird behaviour, bad control of singularity
{
# use dic and aic as bic and aic correspondinly
X <- stats::model.matrix(object = formula,data = data)
out <- speedglm::speedglm.wfit(y = data[,1], X = X, intercept=FALSE, family=family,eigendec = T, method = "Cholesky")
if(prior == "AIC")
  return(list(mlik = -out$aic ,waic = -(out$deviance + 2*out$rank) , dic =  -(out$RSS),summary.fixed =list(mean = out$coefficients)))
if(prior=="BIC")
  return(list(mlik = -out$RSS-logn*out$rank ,waic = -(out$deviance + 2*out$rank) , dic =  -(out$RSS+logn*out$rank),summary.fixed =list(mean = out$coefficients)))
}
