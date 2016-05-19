repeatsamp.bas <- function(k)
  {
   library(MASS)
   trueremprob = matrix(NA,k,10)
   for(j in 1:k)
     {
       print(j)
 multi.bas = bas.lm(y~.,data=multidata,n.models=2^15,prior="g-prior",update=500,alpha=100,initprobs="eplogp")
        ssgrid = (1:10)*2^15/100 # this vector has the elements 1,2,...10% of 2^15
        ssgrid = round(ssgrid,0)# 328  655  983 1311 1638 1966 2294 2621 2949 3277
        marg = exp(multi.bas$logmarg)
        cummarg = cumsum(marg)
        truesummarg = cummarg[length(cummarg)]

        
   for(i in 1:length(ssgrid))
        {
          trueremprob[j,i] = 1-(cummarg[ssgrid[i]]/truesummarg)
        }
    }
     return( list(trueprob=trueremprob))
  }

repeatsamp.srs <- function(k)
  {
   library(MASS)
   trueremprob = matrix(NA,k,10)
   for(j in 1:k)
     {
       print(j)
 multi.bas = bas.lm(y~.,data=multidata,n.models=2^15,prior="g-prior",update=2^15,
   alpha=100,initprobs="Uniform")
        ssgrid = (1:10)*2^15/100 # this vector has the elements 1,2,...10% of 2^15
        ssgrid = round(ssgrid,0)# 328  655  983 1311 1638 1966 2294 2621 2949 3277
        marg = exp(multi.bas$logmarg)
        cummarg = cumsum(marg)
        truesummarg = cummarg[length(cummarg)]
        
   for(i in 1:length(ssgrid))
        {
          trueremprob[j,i] = 1-(cummarg[ssgrid[i]]/truesummarg)
        }
    }
     return( list(trueprob=trueremprob))
  }
