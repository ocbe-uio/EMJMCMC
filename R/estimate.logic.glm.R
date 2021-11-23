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
Jprior <- sum(log(factorial(sj)/((m^sj)*2^(2*sj-2))))
mlik = (-(out$deviance + log(n)*(out$rank)) + 2*(Jprior))/2+n
if(mlik==-Inf)
  mlik = -10000
return(list(mlik = mlik,waic = -(out$deviance + 2*out$rank) , dic =  -(out$deviance + log(n)*out$rank),summary.fixed =list(mean = stats::coefficients(out))))
}
