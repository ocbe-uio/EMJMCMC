#' @title An example of user defined parallelization (cluster based) function
#' for within an MJMCMC chain calculations (mclapply or lapply are used by
#' default depending on specification and OS).
#' @param X a vector (atomic or list) or an expressions vector. Other objects
#' (including classed objects) will be coerced by as.list
#' @param FUN the function to be applied to each element of X or v, or in
#' parallel to X
#' @return \code{parallelize(X,FUN)}, a list of the same length as X and named by X
#' @details Only allowed when working with big.memory based hash table within
#' MJMCMC (see runemjmcmc for more details)
#' @seealso parLapply clusterMap mclapply lapply
#' @keywords methods models
#' @importFrom bigmemory attach.resource
parallelize<-function(X,FUN)
{
max.cpu <- length(X)
cl <-parallel::makeCluster(max.cpu ,type = paral.type,outfile = "")#outfile = ""
parallel::clusterEvalQ(cl = cl,expr = c(library(bigmemory)))
if(exists("statistics"))
{
  parallel::clusterExport(cl=cl, "statistics")
  parallel::clusterEvalQ(cl,{statistics <- attach.resource(statistics);1})
}
res.par <- parallel::parLapply(cl = cl, X, FUN)
parallel::stopCluster(cl)
return(res.par)
}
