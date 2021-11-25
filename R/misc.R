# This file hosts small, internal functions that are used as auxiliary to other,
# user-level functions.

mcgmjpar = function(X,FUN,mc.cores) parallel::mclapply(X= X,FUN = FUN,mc.preschedule = T,mc.cores = mc.cores)

mcgmjpse = function(X,FUN,mc.cores) lapply(X,FUN)

sigmoid<- function(x) {
 return(1/(1+(exp(-x))))
}
