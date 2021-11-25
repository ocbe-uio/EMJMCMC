#' @title A test function to work with elastic networks in future, be omitted so far
#' @param formula a formula object for the model to be addressed
#' @param data a data frame object containing variables and observations corresponding to the formula used
#' @param response response in a formula
#' @param family distribution of the response family object
#' @param alpha regularization parameter in [0,1]
#' @return
#' \describe{
#'  \item{mlik}{marginal likelihood of the model}
#'  \item{waic}{AIC model selection criterion}
#'  \item{dic}{BIC model selection criterion}
#'  \item{summary.fixed$mean}{a vector of posterior modes of the parameters}
#' }
#' @keywords methods models
#' @seealso glmnet::glmnet
#' @export
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
