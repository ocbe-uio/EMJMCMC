# This file hosts small, internal functions that are used as auxiliary to other,
# user-level functions.

mcgmjpar = function(X,FUN,mc.cores) parallel::mclapply(X= X,FUN = FUN,mc.preschedule = T,mc.cores = mc.cores)

mcgmjpse = function(X,FUN,mc.cores) if(.Platform[[1]]=="unix") parallel::mclapply(X= X,FUN = FUN,mc.preschedule = T,mc.cores = mc.cores) else lapply(X,FUN)
