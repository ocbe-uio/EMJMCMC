library(hash)
library(stringi)
#setwd("/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/simulations/simulations")


cosi<-function(x)cos(x/180*pi)
sini<-function(x)sin(x/180*pi)
m<-function(x,y)(x*y)
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
temp = list.files(pattern="posteriorsJA3_*")
myfiles = lapply(FUN = read.csv,X = temp,stringsAsFactors=F)

details = file.info(list.files(pattern="*posteriorsJA1*"))
details = details[with(details, order(as.POSIXct(mtime),decreasing = T)), ]
files = rownames(details)[1:47]
ids<-NULL
nms<-NULL
i<-0
for(file in files)
{
  i<-i+1
  tmp<-strsplit(x = file,fixed = T,split = c("_","."))[[1]][2]
  tmp<-strsplit(x = tmp,fixed = T,split = ".")[[1]][1]
  if(as.integer(tmp)<=100)
  {
    ids<-c(ids,i)
    nms<-c(nms,tmp)
  }
}
temp<-files[ids]
myfiles = lapply(FUN = read.csv,X = temp,stringsAsFactors=F)

#X<- as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))
length(myfiles)

X<-read.csv("exa1.csv")

simplifyposteriors<-function(X,posteriors,th=0.0001,thf=0.3)
{
  dddel<-which(posteriors[,2]<th)
  if(length(dddel)>0)
    posteriors<-posteriors[-dddel,]
  rhash<-hash()
  for(i in 1:length(posteriors[,1]))
  {
    expr<-posteriors[i,1]
    print(expr)
    res<-model.matrix(data=X,object = as.formula(paste0("RadiusJpt~",expr)))
    ress<-c(stri_flatten(res[,1],collapse = ""),stri_flatten(res[,2],collapse = ""),posteriors[i,2],expr)
    if(!((ress[2] %in% values(rhash))))
    {
      rhash[[ress[2]]]<-ress
      print(ress[3])
    }
    else
    {
      if(ress[2] %in% keys(rhash))
      {
        rhash[[ress[2]]][3]<- (as.numeric(rhash[[ress[2]]][3]) + as.numeric(ress[3]))
        if(stri_length(rhash[[ress[2]]][4])>stri_length(expr))
          rhash[[ress[2]]][4]<-expr
      }

    }

  }
  #View(values(rhash))
  res<-as.data.frame(t(values(rhash)[c(3,4),]))
  row.names(res)<-1:dim(res)[1]
  #print(res$V1)
  res[,1]<-as.numeric(as.character(res[,1]))
  res<-res[which(res[,1]>thf),]
  res<-res[order(res[,1], decreasing = T),]
  clear(rhash)
  rm(rhash)
  res[which(res[,1]>1),1]<-1
  colnames(res)<-c("posterior","tree")
  return(res)
}


for(i in 1:length(myfiles))
{
  print(i)
  write.csv(x = simplifyposteriors(X=X,posteriors=as.matrix(myfiles[i][[1]]),th=0.0001,thf=0.1),file =  paste0("postJA32_",nms[i],".csv"),row.names = F)
}
