#' @importFrom bigmemory attach.resource
parallelize<-function(X,FUN)
{
max.cpu <- length(X)
cl <-parallel::makeCluster(max.cpu ,type = paral.type,outfile = "")#outfile = ""
parallel::clusterEvalQ(cl = cl,expr = c(library(INLA),library(bigmemory)))
if(exists("statistics"))
{
  parallel::clusterExport(cl=cl, "statistics")
  parallel::clusterEvalQ(cl,{statistics <- attach.resource(statistics);1})
}
res.par <- parallel::parLapply(cl = cl, X, FUN)
parallel::stopCluster(cl)
return(res.par)
}
