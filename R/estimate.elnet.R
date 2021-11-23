estimate.elnet <- function(formula,response, data, family,alpha)
{
#estimate elastic nets
X <- stats::model.matrix(object = formula,data = data)
if(dim(X)[2]<=2)
  return(list(mlik =  -10000,waic = 10000 , dic =  10000,summary.fixed =list(mean = array(0,dim=dim(X)[2]))))
out <-glmnet::glmnet(x=X[,-1], y = data[[response]], family=family, a=alpha)
dout<- as.numeric(-stats::deviance(out)[[out$dim[2]]])
return(list(mlik = dout,waic = -dout , dic =  -dout,summary.fixed =list(mean = stats::coef(out)[,out$dim[2]])))
}
