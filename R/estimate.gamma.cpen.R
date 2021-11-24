estimate.gamma.cpen <- function(formula, data,r = 1.0/1000.0,logn=log(1000.0),relat=c("cos","sigmoid","tanh","atan","sin","erf"))
{
fparam<-NULL
fmla.proc<-as.character(formula)[2:3]
fobserved <- fmla.proc[1]
fmla.proc[2]<-stringi::stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
fmla.proc[2]<-stringi::stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
fparam <-stringi::stri_split_fixed(str = fmla.proc[2],pattern = "+I",omit_empty = F)[[1]]
sj<-(stringi::stri_count_fixed(str = fparam, pattern = "*"))
sj<-sj+(stringi::stri_count_fixed(str = fparam, pattern = "+"))
for(rel in relat)
  sj<-sj+(stringi::stri_count_fixed(str = fparam, pattern = rel))
#sj<-sj+1
tryCatch(utils::capture.output({
  out <- stats::lm(formula = formula,data = data, family = stats::gaussian)
  # 1 for aic, 2 bic prior, else g.prior

  mlik = (-(stats::BIC(out) -2*log(r)*sum(sj))+1000)/2
  waic = (out$deviance + 2*out$rank)+10000
  dic =  (out$deviance + logn*out$rank)+10000
  summary.fixed =list(mean = stats::coefficients(out))

}, error = function(err) {
  print(err)
  mlik = -10000
  waic = 10000
  dic =  10000
  summary.fixed =list(mean = array(0,dim=length(fparam)))
}))
return(list(mlik = mlik,waic = waic , dic = dic,summary.fixed =summary.fixed))

}

#specify the estimator function returning p(Y|m)p(m), model selection criteria and the vector of the modes for the beta coefficients
# TODO: copied from pinferunemjmcmc example. Diff with above and remove
# estimate.gamma.cpen = function(formula, data,r = 1.0/223.0,logn=log(223.0),relat=c("to23","expi","logi","to35","sini","troot","sigmoid"))
# {
#   fparam=NULL
#   fmla.proc=as.character(formula)[2:3]
#   fobserved = fmla.proc[1]
#   fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
#   fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
#   fparam =stri_split_fixed(str = fmla.proc[2],pattern = "+I",omit_empty = F)[[1]]
#   sj=(stri_count_fixed(str = fparam, pattern = "*"))
#   sj=sj+(stri_count_fixed(str = fparam, pattern = "+"))
#   for(rel in relat)
#     sj=sj+(stri_count_fixed(str = fparam, pattern = rel))
#   sj=sj+1
#   tryCatch(capture.output({
#     out = glm(formula = formula,data = data, family = gaussian)
#     mlik = (-(out$deviance -2*log(r)*sum(sj)))/2
#     waic = -(out$deviance + 2*out$rank)
#     dic =  -(out$deviance + logn*out$rank)
#     summary.fixed =list(mean = coefficients(out))

#   }, error = function(err) {
#     print(err)
#     mlik = -10000
#     waic = -10000
#     dic =  -10000
#     summary.fixed =list(mean = array(0,dim=length(fparam)))
#   }))
#   return(list(mlik = mlik,waic = waic , dic = dic,summary.fixed =summary.fixed))

# }
