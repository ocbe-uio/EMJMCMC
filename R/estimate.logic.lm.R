estimate.logic.lm <- function(formula, data, n, m, r = 1)
{
out <- stats::lm(formula = formula,data = data)
p <- out$rank
fmla.proc<-as.character(formula)[2:3]
fobserved <- fmla.proc[1]
fmla.proc[2]<-stringi::stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
fmla.proc[2]<-stringi::stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
fparam <-stringi::stri_split_fixed(str = fmla.proc[2],pattern = "+",omit_empty = F)[[1]]
sj<-(stringi::stri_count_fixed(str = fparam, pattern = "&"))
sj<-sj+(stringi::stri_count_fixed(str = fparam, pattern = "|"))
sj<-sj+1
Jprior <- prod(factorial(sj)/((m^sj)*2^(2*sj-2)))
#tn<-sum(stringi::stri_count_fixed(str = fmla.proc[2], pattern = "I("))
mlik = (-stats::BIC(out)+2*log(Jprior) + 2*p*log(r)+n)/2
if(mlik==-Inf)
  mlik = -10000
return(list(mlik = mlik,waic = stats::AIC(out)-n , dic =  stats::BIC(out)-n,summary.fixed =list(mean = stats::coef(out))))
}
