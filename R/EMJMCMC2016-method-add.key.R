EMJMCMC2016$methods(
  add.key = function(dec, bit, sum.one, levl) {
    lb <- length(bit)
    if (dec > 2^hash.length) {
      return(FALSE)
    }

    if (levl) {
      if (is.na(statistics1[dec, 1])) {
        statistics1[dec, 16] <- sum.one
        hash.keys1[dec, ] <- bit
        return(TRUE)
      }

      if (is.na(statistics1[dec, 16])) {
        statistics1[dec, 16] <- sum.one
        hash.keys1[dec, ] <- bit # c(array(0,dim = Nvars-lb),bit)
        return(TRUE)
      }

      if (statistics1[dec, 16] != sum.one) {
        return(FALSE)
      }
      i <- 1
      # print(lb)
      lb <- length(bit)
      while (i <= lb && (hash.keys1[dec, Nvars - i + 1]) == bit[lb - i + 1]) {
        i <- i + 1
      }
      return((i - 1) == lb)
    } else {
      if (is.na(statistics[dec, 1])) {
        statistics[dec, 16] <- sum.one
        hash.keys[dec, ] <- bit
        return(TRUE)
      }

      if (is.na(statistics[dec, 16])) {
        statistics[dec, 16] <- sum.one
        hash.keys[dec, ] <- bit # c(array(0,dim = Nvars-lb),bit)
        return(TRUE)
      }

      if (statistics[dec, 16] != sum.one) {
        return(FALSE)
      }
      i <- 1
      lb <- length(bit)
      while (i <= lb && (hash.keys[dec, Nvars - i + 1]) == bit[lb - i + 1]) {
        i <- i + 1
      }
      return((i - 1) == lb)
    }
  }
)
