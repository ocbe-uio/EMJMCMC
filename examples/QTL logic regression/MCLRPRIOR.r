#install.packages("partitions")
#install.packages("matrixStats")
library(partitions)
library(matrixStats)


N_MCLR_prior <- function(s, m,  K = 15,  M = 5){
  
  # number of leaves in model (should become input variable of function)
  # number of markers
  
  #K Maximal number of trees (default 15)
  #M Maximal number of leaves per tree (default 5)
  
  
  if  (s<1 | s > K*M) {
    N = 0
    warning("Not a valid number of leaves s")
  }
  
  
  else if (s==1) {
    N = m
  } 
  
  
  
  else {
    Tab = restrictedparts(s,K)    #Only works for s > 1
    
    colmax = apply(Tab,2,max)
    
    Tab<-Tab[,which(colmax<=M)]
    
    #If figures are getting too large one might work with sums of logarithms already here
    
    if (s > M*K - 2){
      C1 =  prod(factorial(Tab))      # Combinatorial factor from multinomial coefficient (s_1! * ... * s_k!)
      
      C2  = prod(2^(2*Tab-2*(Tab>0)))  # Factor corresponding to number of trees of specific size (2^(2*s_j-2)) 
      
    }
    else {
      C1 =  colProds(factorial(Tab))      # Combinatorial factor from multinomial coefficient (s_1! * ... * s_k!)
      
      C2  = colProds(2^(2*Tab-2*(Tab>0)))  # Factor corresponding to number of trees of specific size (2^(2*s_j-2)) 
    }
    
    N = sum(C2/C1) * prod((m-s+1):m) 
    
    logN =  log(sum(C2/C1)) + sum(log((m-s+1):m))
    
    
  }
  
  return(N)
  
}

m = 100

N_MCLR_prior(0,m)

N_MCLR_prior(1,m)

c(N_MCLR_prior(2,m), 6 * choose(m,2))

c(N_MCLR_prior(3,m),34 * choose(m,3))

c(N_MCLR_prior(4,m),296 * choose(m,4))

c(N_MCLR_prior(5,m), (120 + 4*60 + 16*50 + 64*15 + 256) * choose(m,5))



N_MCLR_prior(25,m)

N_MCLR_prior(50,m)

N_MCLR_prior(74,m)

N_MCLR_prior(75,m)

N_MCLR_prior(90,m,K = 100,M = 5)