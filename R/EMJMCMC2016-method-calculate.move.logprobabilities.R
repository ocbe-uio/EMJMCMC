EMJMCMC2016$methods(
  # calculate move probabilities
  calculate.move.logprobabilities = function(varold, varnew, switch.type, min.N, max.N) {
    if (switch.type == 1) # random size random N(x)
      {
        min.N <<- max.N

        log.mod.switch.prob <- ifelse(max.N < min.N, 0, log(1 / (max.N - min.N + 1))) # probability of having that many differences
        KK <- sum(abs(varold - varnew))

        log.mod.switch.prob <- 0 # always the same probabilities for moves within thenighbourhood for p=0.5
        log.mod.switchback.prob <- 0
      } else if (switch.type == 2) # fixed N(x) inverse operator
      {
        if (min.N != max.N) {
          min.N <<- max.N
        }
        log.mod.switch.prob <- ifelse(max.N < min.N, 0, log(1 / (max.N - min.N + 1)))
        KK <- max.N
        log.mod.switch.prob <- log.mod.switch.prob + KK * log(factorial(Nvars - KK + 1) / factorial(Nvars))
        log.mod.switchback.prob <- log.mod.switch.prob
      } else if (switch.type == 3) # random sized inverse N(x)
      {
        log.mod.switch.prob <- ifelse(max.N < min.N, 0, log(1 / (max.N - min.N + 1)))
        KK <- sum(abs(varold - varnew))
        log.mod.switch.prob <- log.mod.switch.prob + KK * log(factorial(Nvars - KK + 1) / factorial(Nvars))
        log.mod.switchback.prob <- log.mod.switch.prob
      } else if (switch.type == 4) # fixed N(x) for reverse from type 2 swaps
      {
        if (min.N != max.N) {
          min.N <<- max.N
        }
        log.mod.switch.prob <- ifelse(max.N < min.N, 0, log(1 / (max.N - min.N + 1)))
        KK <- max.N
        log.mod.switch.prob <- log.mod.switch.prob + KK * log(factorial(Nvars - KK + 1) / factorial(Nvars))
        log.mod.switchback.prob <- log.mod.switch.prob
      } else if (switch.type > 4) {
      log.mod.switch.prob <- log(x = 1)
      log.mod.switchback.prob <- log(x = 1)
    }
    return(list(log.switch.forw.prob = log.mod.switch.prob, log.switch.back.prob = log.mod.switchback.prob))
  }
)
