EMJMCMC2016$methods(
  # calculates posterior probabilities based on a current search
  post_proceed_results = function(statistics1) {
    xyz <- which(!is.na(statistics1[, 1]))
    g.results[4, 2] <<- length(xyz)
    xyz <- intersect(xyz, which(statistics1[, 1] != -10000))
    moddee <- which(statistics1[, 1] == max(statistics1[, 1], na.rm = TRUE))[1]
    zyx <- array(data = NA, dim = length(statistics1[, 1]))
    nconsum <- sum(exp(-statistics1[moddee, 1] + statistics1[xyz, 1]), na.rm = TRUE)

    if (nconsum > 0) {
      zyx[xyz] <- exp(statistics1[xyz, 1] - statistics1[moddee, 1]) / nconsum
    } else {
      nnnorm <- sum(statistics1[xyz, 4], na.rm = T)
      if (nnnorm == 0) {
        nnnorm <- 1
      }
      zyx[xyz] <- statistics1[xyz, 4] / nnnorm
    }
    statistics1[, 15] <- zyx

    lldd <- 2^(Nvars) + 1
    p.post <- array(data = 0, dim = Nvars)
    for (i in xyz)
    {
      vec <- dectobit(i - 1)
      varcur <- c(array(0, dim = (Nvars - length(vec))), vec)
      p.post <- (p.post + varcur * statistics1[i, 15])
    }

    if (!exists("p.post") || is.null(p.post) || sum(p.post, na.rm = T) == 0 || sum(p.post, na.rm = T) > Nvars) {
      p.post <- array(data = 0.5, dim = Nvars)
    }

    return(list(p.post = p.post, m.post = zyx, s.mass = sum(exp(statistics1[xyz, 1]), na.rm = TRUE)))
  }
)
