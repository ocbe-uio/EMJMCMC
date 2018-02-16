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
temp = list.files(pattern="posteriorsLog1SingRMI*")
myfiles = lapply(FUN = read.csv,X = temp,stringsAsFactors=F)

#X<- as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))
length(myfiles)

X<-read.csv("exa1.csv")


simplifyposteriors<-function(X,posteriors,th=0.0001,thf=0.3)
{
  ord<-order(posteriors[,2],decreasing = T)
  th<-posteriors[ord[4],2]
  posteriors<-posteriors[ord,]
  print(th)
  posteriors<-posteriors[-which(posteriors[,2]<th),]
  
  res<-cbind(posteriors[,2],posteriors[,1])
  res<-res[which(res[,1]>thf),]
  colnames(res)<-c("posterior","tree")
  return(res)
}
for(i in 1:length(myfiles))
{
  print(i)
  myfiles[i][[1]]<-simplifyposteriors(X=X,posteriors=myfiles[i][[1]],th=0.0001,thf=0.1)
}


for(i in 1:length(myfiles))
{
  write.csv(x= myfiles[i],file = paste0("postJA32_",temp[i]))
}

res=round(runif(4,1,96))
for(i in res)
{
  write.csv(x= myfiles[i],file = paste0("postJA32_",temp[i]))
}

clear(rhash)
rm(rhash)


rhash<-hash()

N<-length(myfiles)
alpha<-0.3
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


write.csv(x = values(rhash),file = "expJA33.csv")





