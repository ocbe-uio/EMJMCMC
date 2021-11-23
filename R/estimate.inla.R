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
