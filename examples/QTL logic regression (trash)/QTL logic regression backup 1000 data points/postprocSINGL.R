library(hash)
library(stringi)
#setwd("/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/simulations/simulations")


simplifyposteriors<-function(X,posteriors,th=0.0001,thf=0.5)
{
  lw<-which(posteriors[,2]<th)
  if(length(lw)>0)
    posteriors<-posteriors[-lw,]
  rhash<-hash()
  for(i in 1:length(posteriors[,1]))
  {
    expr<-posteriors[i,1]
    print(expr)
    res<-model.matrix(data=X,object = as.formula(paste0("Y1~",expr)))
    res[,1]<-res[,1]-res[,2]
    ress<-c(stri_flatten(res[,1],collapse = ""),stri_flatten(res[,2],collapse = ""),posteriors[i,2],expr)
    if(!(ress[1] %in% values(rhash)||(ress[2] %in% values(rhash))))
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
  res<-as.data.frame(t(values(rhash)[c(3,4),]))
  res$V1<-as.numeric(as.character(res$V1))
  res$V2<-as.character(res$V2)
  res<-res[which(res$V1>thf),]
  res<-res[order(res$V1, decreasing = T),]
  clear(rhash)
  rm(rhash)
  res[which(res[,1]>1),1]<-1
  colnames(res)<-c("posterior","tree")
  return(res)
}

X1<- as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = 0.3),dim = c(1000,50)))
Y1=-0.45+0.6*((1-X1$V1)*(X1$V4)) + 0.6*(X1$V8*X1$V11)+0.6*(X1$V5*X1$V9)
X1$Y1<-round(1.0/(1.0+exp(-Y1)))


#experiment i
temp = list.files(pattern="posteriorsLog3SingMCI*")
myfiles = lapply(FUN = read.csv,X = temp,stringsAsFactors=F)[1:100]


#X<- as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))
length(myfiles)

for(i in 1:length(myfiles))
{
  print(i)
  myfiles[i][[1]]<-simplifyposteriors(X=X1,posteriors=myfiles[i][[1]],th=0.0001,thf=0.01)
}


for(i in 1:length(myfiles))
{
  myfiles[i]$tree<-as.character(myfiles[i]$tree)
  write.csv(x= myfiles[i],file = paste0("postLogSing_",temp[i]))
}


N<-length(myfiles)
X<-X1
alpha<-0.5
rhash<-hash()
clear(rhash)
for(i in 1:min(100,N))
{
  for(j in 1:length(myfiles[[i]]$posterior))
  {
    if(myfiles[[i]]$posterior[j]>=alpha)
    {
      expr<-as.character(myfiles[[i]]$tree[j])
      print(expr)
      res<-model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
      res[,1]<-res[,1]-res[,2]
      ress<-c(stri_flatten(res[,1],collapse = ""),stri_flatten(res[,2],collapse = ""),1,expr)
      if(!(ress[1] %in% values(rhash)||(ress[2] %in% values(rhash))))
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

res<-as.data.frame(t(values(rhash)[c(4,3),]))
res$V1<-as.numeric(as.character(res$V1))
res<-res[order(res$V1, decreasing = T),]
colnames(res)<-c("posterior","tree")
write.csv(x = values(rhash),file = "exp1Sing3MC.csv")
clear(rhash)
rm(rhash)



