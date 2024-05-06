EMJMCMC2016$methods(
  hashing = function(bit) # a hash function to find where to place the key in the hash
  {
    n <- length(bit)
    if (n < hash.length) {
      return(bittodec.alt(bit))
    }
    return(bittodec.alt(bit[(n - hash.length + 1):n]))
  }
)
