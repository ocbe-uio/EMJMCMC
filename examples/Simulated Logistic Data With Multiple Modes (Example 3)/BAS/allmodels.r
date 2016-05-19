# trying to create a matrix which lists all 2^p possible models in terms of indicators
allmodels=function(p=3)
  {
    mat=matrix(0,2^p,p)
     binary=function (x, dim)
{
    if (x == 0) {
        pos <- 1
    }
    else {
        pos <- floor(log(x, 2)) + 1
    }
    if (!missing(dim)) {
        if (pos <= dim) {
            pos <- dim
        }
        else {
            warning("the value of `dim` is too small")
        }
    }
    bin <- rep(0, pos)
    dicotomy <- rep(FALSE, pos)
    for (i in pos:1) {
        bin[i] <- floor(x/2^(i - 1))
        dicotomy[i] <- bin[i] == 1
        x <- x - ((2^(i - 1)) * bin[i])
    }
    bin=rev(bin)
    return( bin)
}
for(i in 1:(2^p))
  mat[i,]=binary((i-1),p)
return(mat)
  }
