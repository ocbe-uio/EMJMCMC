estimate.inla.poisson <- function(formula, data,r = 1.0/200.0,logn=log(200.0), relat=c("*","+","cos","sigmoid","tanh","atan","sin","erf"))
{


out<-NULL
utils::capture.output({tryCatch(utils::capture.output({
  fmla.proc<-as.character(formula)[2:3]
  sj<-sum(stringi::stri_count_fixed(str = fmla.proc[2],pattern = relat))
  out <-INLA::inla(family = "poisson",silent = 2L,data = data,formula = formula,control.compute = list(dic = TRUE, waic = TRUE, mlik = TRUE))
    }))})
if(is.null(out))
  return(list(mlik = -10000+log(r)*(sj),waic =  10000 , dic = 10000, summary.fixed =list(mean = NULL)))

if(length(out$waic[1]$waic)==0)
  return(list(mlik = -10000+log(r)*(sj),waic =  10000 , dic = 10000, summary.fixed =list(mean = NULL)))
# use dic and aic as bic and aic correspondinly

coef<-out$summary.fixed$mode
#coef[1]<-coef[1]+out$summary.hyperpar$mode[1]
return(list(mlik = out$mlik[1]+log(r)*(sj),waic =  out$waic[1]$waic , dic = out$dic[1]$dic, summary.fixed =list(mean = coef)))

}
