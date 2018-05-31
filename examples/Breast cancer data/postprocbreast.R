library(hash)
library(stringi)
setwd("/mn/sarpanitu/ansatte-u2/aliaksah/abeldata/breast cancer/parres/")



details = file.info(list.files(pattern="posteriorsbreast_*"))
details = details[with(details, order(as.POSIXct(mtime),decreasing = T)), ]
files = rownames(details)


myfiles = lapply(FUN = read.csv,X = files,stringsAsFactors=F)


length(myfiles)

source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")
troot<-function(x)abs(x)^(1/3)
sini<-function(x)sin(x/180*pi)
logi<-function(x)log(abs(x+0.1))
gfquar<-function(x)as.integer(x<quantile(x,probs = 0.25))
glquar<-function(x)as.integer(x>quantile(x,probs = 0.75))
gmedi<-function(x)as.integer(x>median(x))
cosi<-function(x)cos(x/180*pi)
gmean<-function(x)as.integer(x>mean(x))
gone<-function(x)as.integer(x>0)
gthird<-function(x)(abs(x)^(1/3))
gfifth<-function(x)(abs(x)^(1/5))
grelu<-function(x)(x*(x>0))
contrelu<-function(x)log(1+exp(x))
gauss<-function(x)exp(-x*x)
simplifyposteriors<-function(X,posteriors,th=0.0001,thf=0.1)
{
  posteriors<-posteriors[-which(posteriors[,2]<th),]
  rhash<-hash()
  for(i in 1:length(posteriors[,1]))
  {


    expr<-posteriors[i,1]
    print(expr)
    res<-model.matrix(data=X,object = as.formula(paste0("X~",expr)))
    ress<-c(stri_flatten(round(res[,2],digits = 4),collapse = ""),stri_flatten(res[,2],collapse = ""),posteriors[i,2],expr)
    if(!((ress[1] %in% values(rhash))))
      rhash[[ress[1]]]<-ress
    else
    {
      if(ress[1] %in% keys(rhash))
      {
        rhash[[ress[1]]][3]<- (as.numeric(rhash[[ress[1]]][3]) + as.numeric(ress[3]))
        if(stri_length(rhash[[ress[1]]][4])>stri_length(expr))
          rhash[[ress[1]]][4]<-expr
      }
    }

  }
  res<-as.data.frame(t(values(rhash)[c(3,4),]))
  res$V1<-as.numeric(as.character(res$V1))
  res<-res[which(res$V1>thf),]
  res<-res[order(res$V1, decreasing = T),]
  clear(rhash)
  rm(rhash)
  res[which(res[,1]>1),1]<-1
  colnames(res)<-c("posterior","tree")

  row.names(res) = 1:length(res$posterior)
  return(res)
}



X<-read.csv("/mn/sarpanitu/ansatte-u2/aliaksah/abeldata/breast cancer/train.csv")

featgmj = hash()

for(i in 1:100)
{
  tmp = simplifyposteriors(X = X,posteriors = myfiles[[i]])
  for(feat in tmp$tree)
  {
    if(!has.key(hash = featgmj,key =  feat ))
      {
        featgmj[[feat]] = as.numeric(1)
      } else{

        featgmj[[feat]] =as.numeric(featgmj[[feat]]) + 1
      }
  }
}


write.csv(x = cbind(keys(featgmj),values(featgmj)),file = paste0("breastparfeatrgmj",".csv"))
