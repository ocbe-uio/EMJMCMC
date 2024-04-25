# transform binary numbers to decimal
EMJMCMC2016$methods(
  bittodec.alt = function(bit) # transform a binary vector into a natural number to correspond between vector of solutions and storage array
  {
    n <- length(bit)
    dec <- 0
    for (i in 1:n)
    {
      j <- n - i
      dec <- dec + ((2)^j) * bit[i]
    }
    return(dec)
  }
)
EMJMCMC2016$methods(
  bittodec = function(bit) # transform a binary vector into a natural number to correspond between vector of solutions and storage array
  {
    if (!double.hashing) {
      n <- length(bit)
      dec <- 0
      for (i in 1:n)
      {
        j <- n - i
        dec <- dec + ((2)^j) * bit[i]
      }

      return(dec)
    } else {
      if (exists("statistics1")) {
        hash.level <- 0
        dec <- hashing(bit) + 1
        # print(dec)
        sum.one <- sum(bit) * which.max(bit) + sum(bit[Nvars - 7:Nvars + 1]) * Nvars
        jjj <- 1
        while (!add.key(dec, bit, sum.one, TRUE)) {
          hash.level <- hash.level + 1
          bit1 <- dectobit.alt(2654435761 * (dec + 97 * sum.one + hash.level * 36599) + hash.level * 59 + hash.level)
          dec <- hashing(bit1) + 1
          # jjj<-jjj+1
          # print(dec)
        }
        # print(jjj)
        dec <- dec - 1
      } else if (exists("statistics")) {
        hash.level <- 0
        dec <- hashing(bit) + 1
        sum.one <- sum(bit) * which.max(bit) + sum(bit[Nvars - 7:Nvars + 1]) * Nvars
        while (!add.key(dec, bit, sum.one, FALSE)) {
          hash.level <- hash.level + 1
          bit1 <- dectobit.alt(2654435761 * (dec + 97 * sum.one + hash.level * 36599) + hash.level * 59 + hash.level)
          dec <- hashing(bit1) + 1
          # print(dec)
        }
        dec <- dec - 1
      }
      return(dec)
    }
  }
)
