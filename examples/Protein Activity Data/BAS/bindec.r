 bindec = function(x) # x is the binary no.(entered as a vector) to be converted 
    {
      p = length(x)
      base = rep(0,p)
      for(i in 1:p) base[i] =2^(i-1)
      base = rev(base)
      dec = sum(x*base)
      return(dec)
    }
