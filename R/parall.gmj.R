parall.gmj <- function(X, M = 16, preschedule = FALSE) {
  parallel::mclapply(
    X              = X,
    FUN            = do.call.emjmcmc,
    mc.preschedule = preschedule,
    mc.cores       = M,
    mc.cleanup     = TRUE
  )
}
