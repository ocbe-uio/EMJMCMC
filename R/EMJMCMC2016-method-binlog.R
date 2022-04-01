EMJMCMC2016$methods(
  binlog = function(x) # bitwise logorithm (base 2) computations based on Al Kashi's algorithm
  {
    lx <- length(x)
    tol <- -lx
    y <- 0 # initialise output
    b <- 0.5 # initialise mantissa


    # arrange the input into a known range
    if (lx == 1) {
      if (x[1] == 0) {
        return(-Inf)
      }

      if (x[1] == 1) {
        return(0)
      }
    }
    # move one bit to the left end
    # elsewhise

    if (x[2] == 0 && lx == 2) {
      return(1)
    }

    powto <- 2^(-c(1:(lx - 1)))
    float.x <- sum(x[2:lx] * powto)
    x <- 1 + float.x
    y <- lx - 1


    f <- 0
    fb <- -1
    # move one bit to the right end
    # now x = 1.5
    # loop until desired tolerance met

    while (fb > tol) {
      x <- x * x

      # update the index
      if (x >= 2) {
        x <- x / 2
        y <- y + b
        f <- log(exp(f) + exp(fb))
      }
      # scale for the next bit
      b <- b / 2
      fb <- (fb - log(2))
    }
    return(list(y = y, z = lx - 1, f = f))
  }
)
