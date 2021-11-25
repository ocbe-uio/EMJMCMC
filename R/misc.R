# This file hosts small, internal functions that are used as auxiliary to other,
# user-level functions.

m<-function(a,b)a*b

factorial<-function(x) ifelse(x<=170,gamma(abs(x) + 1),gamma(171))

mcgmjpar = function(X,FUN,mc.cores) parallel::mclapply(X= X,FUN = FUN,mc.preschedule = T,mc.cores = mc.cores)

mcgmjpse = function(X,FUN,mc.cores) lapply(X,FUN)

sigmoid<- function(x) {
 return(1/(1+(exp(-x))))
}

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
