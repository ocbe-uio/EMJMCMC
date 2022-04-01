EMJMCMC2016$methods(
# transform decimal numbers to binary
  dectobit = function(dec) # transform a natural number into a binary vector to correspond between vector of solutions and storage array
  {
    if (!double.hashing) {
      if (dec == 0) {
        return(0)
      }
      q <- dec
      bin <- NULL
      while (q != 0) {
        r <- q / 2
        q <- floor(r)
        bin <- c(as.integer(2 * (r - q)), bin)
      }
      return(bin)
    }

    return(dehash(dec + 1))
  }
)
EMJMCMC2016$methods(
  dectobit.alt = function(dec) # transform a natural number into a binary vector to correspond between vector of solutions and storage array
  {
    if (dec == 0) {
      return(0)
    }
    q <- dec
    bin <- NULL
    while (q != 0) {
      r <- q / 2
      q <- floor(r)
      bin <- c(as.integer(2 * (r - q)), bin)
    }
    return(bin)
  }
)
