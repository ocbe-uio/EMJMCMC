library(hash)
library(stringi)
#setwd("/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/simulations/simulations")


cosi<-function(x)cos(x/180*pi)
sini<-function(x)sin(x/180*pi)
expi<-function(x)
{
  r<-exp(x)
  if(r==Inf)
    return(10000000)
  else
    return(r)
}

InvX<-function(x)
{
  if(x==0)
    return(10000000)
  else
    return(1/x)
  
}
troot<-function(x)abs(x)^(1/3)
sigmoid<-function(x)exp(-x)


#experiment i
temp = list.files(pattern="postJA3_*")
myfiles = lapply(FUN = read.csv,X = temp,stringsAsFactors=F)

#X<- as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))
length(myfiles)

X<-read.csv("exa1.csv")

simplifyposteriors<-function(X,posteriors,th=0.0001,thf=0.3)
{
  dddel<-which(posteriors[,1]<th)
  if(length(dddel)>0)
    posteriors<-posteriors[-dddel,]
  rhash<-hash()
  for(i in 1:length(posteriors[,2]))
  {
    expr<-posteriors[i,2]
    print(expr)
    res<-model.matrix(data=X,object = as.formula(paste0("RadiusJpt~",expr)))
    ress<-c(stri_flatten(res[,1],collapse = ""),stri_flatten(res[,2],collapse = ""),posteriors[i,1],expr)
    if(!((ress[2] %in% values(rhash))))
      rhash[[ress[1]]]<-ress
    else
    {
      if(ress[1] %in% keys(rhash))
      {
        rhash[[ress[1]]][3]<- (as.numeric(rhash[[ress[1]]][3]) + as.numeric(ress[3]))
        if(stri_length(rhash[[ress[1]]][4])>stri_length(expr))
          rhash[[ress[1]]][4]<-expr
      }
      else
      {
        rhash[[ress[2]]][3]<- (as.numeric(rhash[[ress[2]]][3]) + as.numeric(ress[3]))
        if(stri_length(rhash[[ress[2]]][4])>stri_length(expr))
          rhash[[ress[2]]][4]<-expr
      }
    }
    
  }
  #View(values(rhash))
  res<-as.data.frame(t(values(rhash)[c(3,4),]))
  print(res$V1)
  res$V1<-as.numeric(as.character(res$V1))
  res<-res[which(res$V1>thf),]
  res<-res[order(res$V1, decreasing = T),]
  clear(rhash)
  rm(rhash)
  res[which(res[,1]>1),1]<-1
  colnames(res)<-c("posterior","tree")
  return(res)
}


for(i in 1:length(myfiles))
{
  print(i)
  myfiles[i][[1]]<-simplifyposteriors(X=X,posteriors=myfiles[i][[1]],th=0.0001,thf=0.3)
}





rhash<-hash()

N<-length(myfiles)
alpha<-0.5
clear(rhash)




for(i in 1:min(100,N))
{
  for(j in 1:length(myfiles[[i]]$posterior))
  {
    if(myfiles[[i]]$posterior[j]>=alpha)
    {
      expr<-as.character(myfiles[[i]]$tree[j])
      print(expr)
      res<-abs(model.matrix(data=X,object = as.formula(paste0("RadiusJpt~",expr))))
      ress<-c(stri_flatten(res[,2],collapse = ""),stri_flatten(res[,1],collapse = ""),1,expr)
      if(!(ress[1] %in% values(rhash)))
        rhash[[ress[1]]]<-ress
      else
      {
        if(ress[1] %in% keys(rhash))
        {
          rhash[[ress[1]]][3]<- (as.numeric(rhash[[ress[1]]][3])) + as.numeric(1)
          if(stri_length(rhash[[ress[1]]][4])>stri_length(expr))
            rhash[[ress[1]]][4]<-expr
        }
        else
        {
          rhash[[ress[2]]][3]<- (as.numeric(rhash[[ress[2]]][3])) + as.numeric(1)
          if(stri_length(rhash[[ress[2]]][4])>stri_length(expr))
            rhash[[ress[2]]][4]<-expr
        }
      }
    }
    
  }
  
}


write.csv(x = values(rhash),file = "expJA322.csv")





