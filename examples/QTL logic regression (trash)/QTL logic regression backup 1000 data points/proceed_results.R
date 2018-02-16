library(hash)
library(stringi)

simplifyposteriors<-function(X,posteriors,th=0.0001,thf=0.5)
{
  posteriors<-posteriors[-which(posteriors[,2]<th),]
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
  res<-res[which(res$V1>thf),]
  res<-res[order(res$V1, decreasing = T),]
  clear(rhash)
  rm(rhash)
  res[which(res[,1]>1),1]<-1
  colnames(res)<-c("posterior","tree")
  return(res)
}

# proceed experiment 1
X<-read.csv(file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/X1.csv")[,-1]
posteriors<-read.csv(file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/sim1post.csv",stringsAsFactors = F)
th<-(10)^(-5)
thf<-0.05
res1<-simplifyposteriors(X,posteriors, th,thf)
write.csv(x =res1,row.names = F,file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/post1eta.csv")

# proceed experiment 2
X<-read.csv(file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/X2.csv")[,-1]
posteriors<-read.csv(file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/sim2post.csv",stringsAsFactors = F)
th<-(10)^(-5)
thf<-0.05
res1<-simplifyposteriors(X,posteriors, th,thf)
write.csv(x =res1,row.names = F,file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/post2eta.csv")

# proceed experiment 3
X<-read.csv(file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/X3.csv")[,-1]
posteriors<-read.csv(file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/sim3post.csv",stringsAsFactors = F)
th<-(10)^(-5)
thf<-0.05
res1<-simplifyposteriors(X,posteriors, th,thf)
write.csv(x =res1,row.names = F,file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/post3eta.csv")

# proceed experiment 4
X<-read.csv(file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/X4.csv")[,-1]
posteriors<-read.csv(file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/sim4post.csv",stringsAsFactors = F)
th<-(10)^(-5)
thf<-0.05
res1<-simplifyposteriors(X,posteriors, th,thf)
write.csv(x =res1,row.names = F,file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/post4eta.csv")

# proceed experiment 5
X<-read.csv(file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/X5.csv")[,-1]
posteriors<-read.csv(file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/sim5post.csv",stringsAsFactors = F)
th<-(10)^(-5)
thf<-0.05
res1<-simplifyposteriors(X,posteriors, th,thf)
write.csv(x =res1,row.names = F,file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/post5eta.csv")

