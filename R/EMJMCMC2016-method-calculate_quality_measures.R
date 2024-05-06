EMJMCMC2016$methods(
  calculate_quality_measures = function(vect, n, truth) {
    rmse.pi <- array(data = 0, Nvars)
    bias.pi <- array(data = 0, Nvars)
    for (i in 1:n)
    {
      bias.pi <- (bias.pi + (vect[, i] - truth))
      rmse.pi <- (rmse.pi + (vect[, i]^2 + truth^2 - 2 * vect[, i] * truth))
    }
    bias.pi <- bias.pi / n
    rmse.pi <- rmse.pi / n
    return(list(bias.pi = bias.pi, rmse.pi = rmse.pi))
  }
)
