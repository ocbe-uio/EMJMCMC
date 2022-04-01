EMJMCMC2016$methods(
  # norm between probabilities vectors
  normprob = function(p1, p2) {
    nn <- abs(1 - sum((p1 + 0.1) / (p2 + 0.1)) / length(p1))
    if (is.na(nn)) {
      nn <- Inf
    }
    return(nn)
  }
)
