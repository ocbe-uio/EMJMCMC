EMJMCMC2016$methods(
    dehash = function(dec) {
    if (exists("hash.keys1")) {
      return(hash.keys1[dec, ])
    }
    if (exists("hash.keys")) {
      return(hash.keys[dec, ])
    }
    return(NULL)
  }
)
