EMJMCMC2016$methods(
  # lambda function for mtmcmc
  lambda = function(c, alpha, g1, g2, g.domain.pos) # simmetric choice driving function
  {
    if ((c != 0)) {
      res <- ((((1 / c) * (1 + (g1 + g2) * (g.domain.pos)) / (1 - (g1 + g2) * (1 - g.domain.pos)))^alpha))
      # !#if(printable.opt)print(res)
      return(res)
    } else {
      # !#if(printable.opt)print(1)
      return(1)
    }
  }
)
