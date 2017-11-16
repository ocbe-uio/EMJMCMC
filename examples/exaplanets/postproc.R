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
#temp = list.files(pattern="posteriorsJA3_*")
#myfiles = lapply(FUN = read.csv,X = temp,stringsAsFactors=F)

details = file.info(list.files(pattern="postJA1_*"))
details = details[with(details, order(as.POSIXct(mtime),decreasing = T)), ]
files = rownames(details)



ids<-NULL
nms<-NULL
i<-0
for(file in files)
{
  i<-i+1
  tmp<-strsplit(x = file,fixed = T,split = c("_","."))[[1]][2]
  tmp<-strsplit(x = tmp,fixed = T,split = ".")[[1]][1]
  if(as.integer(tmp)<=150&&stri_count_fixed(str = file,pattern = "new")[[1]]==0&&stri_count_fixed(str = file,pattern = "REV")[[1]]==0&&stri_count_fixed(str = file,pattern = "S")[[1]]==0)
  {
    ids<-c(ids,i)
    nms<-c(nms,tmp)
  }
}
temp<-files[ids]
myfiles = lapply(FUN = read.csv,X = temp,stringsAsFactors=F)[1:100]

X4<- as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))
Y4<-rnorm(n = 1000,mean = 1+7*(X4$V4*X4$V17*X4$V30*X4$V10)+7*(((X4$V50*X4$V19*X4$V13*X4$V11)>0)) + 9*(X4$V37*X4$V20*X4$V12)+ 7*(X4$V1*X4$V27*X4$V3)
          +3.5*(X4$V9*X4$V2) + 6.6*(X4$V21*X4$V18) + 1.5*X4$V7 + 1.5*X4$V8,sd = 1)
X4$Y4<-Y4

length(myfiles)

X<-read.csv("exa1.csv")

aggreg<-NULL
for(i in 1:length(myfiles))
{
  print(i)
  aggreg <- rbind(aggreg,myfiles[i][[1]])
  #write.csv(x = simplifyposteriors(X=X,posteriors=as.matrix(myfiles[i][[1]]),th=0.0001,thf=0.1),file =  paste0("postJA32_",nms[i],".csv"),row.names = F)
}


#xxx<-simplifyposteriors(X=X,posteriors=as.matrix(myfiles[i][[1]]),th=0.0001,thf=0.3)


rhash<-hash()

N<-length(myfiles)
alpha<-0.25
clear(rhash)




for(i in 1:min(100,N))
{
  for(j in 1:length(myfiles[[i]]$posterior))
  {
    if(myfiles[[i]]$posterior[j]>=alpha)
    {
      expr<-as.character(myfiles[[i]]$tree[j])
      print(expr)
      res<-model.matrix(data=X,object = as.formula(paste0("PeriodDays~",expr)))
      ress<-c(stri_flatten(round(res[,2],digits = 4),collapse = ""),stri_flatten(res[,1],collapse = ""),1,expr)
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


write.csv(x = t(values(rhash)[c(3,4),]),file = "exppaap.csv",row.names = F,col.names = F)




rhash<-hash()

N<-length(myfiles)
alpha<-0.25
clear(rhash)

TPS<-c(stri_flatten(round(model.matrix(data=X,object = as.formula(paste0("RadiusJpt~","I(troot(I(I(I(PeriodDays)*I(PeriodDays))*I(HostStarMassSlrMass))))")))[,2],digits = 4),collapse = ""),stri_flatten(round(model.matrix(data=X,object = as.formula(paste0("RadiusJpt~","I(troot(I(I(I(HostStarRadiusSlrRad)*I(PeriodDays))*I(PeriodDays))))")))[,2],digits = 4),collapse = ""),stri_flatten(round(model.matrix(data=X,object = as.formula(paste0("RadiusJpt~","I(troot(I(I(I(PeriodDays)*I(PeriodDays))*I(HostStarTempK))))")))[,2],digits = 4),collapse = ""))

for(i in 1:min(100,N))
{
  j=1
  while(j <= length(myfiles[[i]]$posterior))
  {
    if(myfiles[[i]]$posterior[j]>=alpha)
    {
      expr<-as.character(myfiles[[i]]$tree[j])
      print(expr)
      res<-model.matrix(data=X,object = as.formula(paste0("RadiusJpt~",expr)))
      ress<-c(stri_flatten(round(res[,2],digits = 4),collapse = ""),stri_flatten(res[,1],collapse = ""),1,expr)
      if(ress[1] %in% TPS)
      {
        j = length(myfiles[[i]]$posterior)+1
        #print(ress[1])
      }
      
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
    j = j + 1
  }

}
write.csv(x = t(values(rhash)[c(3,4),]),file = "expJA1p2.csv",row.names = F,col.names = F)
