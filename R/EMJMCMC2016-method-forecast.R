EMJMCMC2016$methods(
  forecast = function(covariates, nvars, link.g) {
    ids <- which(!is.na(statistics1[, 15]))
    res <- 0
    for (i in ids)
    {
      res <- res + statistics1[i, 15] * link.g(sum(statistics1[i, 16:nvars] * covariates, na.rm = T))
    }
    return(list(forecast = res))
  }
)
