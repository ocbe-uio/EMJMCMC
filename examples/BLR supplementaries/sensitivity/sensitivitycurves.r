library("stringi")
library("parallel")

#dead in the files for sensitivity analysis of the 4 types from sensitivity R and the two priors
temp = list.files(pattern="*_T1J_*")
temp = temp[order(x = unlist(mclapply(FUN = function(str){
  as.integer(stri_split_fixed(str = stri_split_fixed(str = str,pattern = c("_"))[[1]][3],pattern = ".")[[1]][1]
  )},X = temp)),decreasing = F)]
files1J = mclapply(FUN = read.csv,X = temp)

temp = list.files(pattern="*_T1G_*")
temp = temp[order(x = unlist(mclapply(FUN = function(str){
  as.integer(stri_split_fixed(str = stri_split_fixed(str = str,pattern = c("_"))[[1]][3],pattern = ".")[[1]][1]
  )},X = temp)),decreasing = F)]
files1G = mclapply(FUN = read.csv,X = temp)

temp = list.files(pattern="*_T1GB_*")
temp = temp[order(x = unlist(mclapply(FUN = function(str){
  as.integer(stri_split_fixed(str = stri_split_fixed(str = str,pattern = c("_"))[[1]][3],pattern = ".")[[1]][1]
  )},X = temp)),decreasing = F)]
files1GB = mclapply(FUN = read.csv,X = temp)

temp = list.files(pattern="*_T2J_*")
temp = temp[order(x = unlist(mclapply(FUN = function(str){
  as.integer(stri_split_fixed(str = stri_split_fixed(str = str,pattern = c("_"))[[1]][3],pattern = ".")[[1]][1]
  )},X = temp)),decreasing = F)]
files2J = mclapply(FUN = read.csv,X = temp)

temp = list.files(pattern="*_T2G_*")
temp = temp[order(x = unlist(mclapply(FUN = function(str){
  as.integer(stri_split_fixed(str = stri_split_fixed(str = str,pattern = c("_"))[[1]][3],pattern = ".")[[1]][1]
  )},X = temp)),decreasing = F)][-c(21:29)][-c(51,62,73,84,95)]
files2G = mclapply(FUN = read.csv,X = temp)

temp = list.files(pattern="*_T3J_*")
temp = temp[order(x = unlist(mclapply(FUN = function(str){
  as.integer(stri_split_fixed(str = stri_split_fixed(str = str,pattern = c("_"))[[1]][3],pattern = ".")[[1]][1]
  )},X = temp)),decreasing = F)]
files3J = mclapply(FUN = read.csv,X = temp)

temp = list.files(pattern="*_T3G_*")
temp = temp[order(x = unlist(mclapply(FUN = function(str){
  as.integer(stri_split_fixed(str = stri_split_fixed(str = str,pattern = c("_"))[[1]][3],pattern = ".")[[1]][1]
  )},X = temp)),decreasing = F)]
files3G = mclapply(FUN = read.csv,X = temp)


temp = list.files(pattern="*_T4J_*")
temp = temp[order(x = unlist(mclapply(FUN = function(str){
  as.integer(stri_split_fixed(str = stri_split_fixed(str = str,pattern = c("_"))[[1]][3],pattern = ".")[[1]][1]
  )},X = temp)),decreasing = F)]
files4J = mclapply(FUN = read.csv,X = temp)

temp = list.files(pattern="*_T4G_*")
temp = temp[order(x = unlist(mclapply(FUN = function(str){
  as.integer(stri_split_fixed(str = stri_split_fixed(str = str,pattern = c("_"))[[1]][3],pattern = ".")[[1]][1]
  )},X = temp)),decreasing = F)]
files4G = mclapply(FUN = read.csv,X = temp)

#create a dummy matrix of covariates
X= as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))

#create dummy variables corresponding to the true positives under the given scenarios
res.tpl.4=model.matrix(data=X,object = as.formula(paste0("V1~","I(V4*V17*V30*V10)")))
res.tpl.3=model.matrix(data=X,object = as.formula(paste0("V1~","I(V7*V12*V20)")))
res.tpl.2=model.matrix(data=X,object = as.formula(paste0("V1~","I(V2*V9)")))
res.tpl.1=model.matrix(data=X,object = as.formula(paste0("V1~","I(V37)")))

#make a plot for varying effect sizes for Jeffrey's and robust g prior
post.b4 = array(0,100)
post.b3 = array(0,100)
post.b2 = array(0,100)
post.b1 = array(0,100)
for(i in 1:100)
{
  for(j in 1:length(files1J[[i]]$posterior))
  {
    
    expr=as.character(files1J[[i]]$tree[j])
    res=model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
    if(sum(abs(res[,2]-res.tpl.4[,2]))==0)
    {
      post.b4[i] = as.numeric(as.character(files1J[[i]]$posterior[j]))
    }else if(sum(abs(res[,2]-res.tpl.3[,2]))==0)
    {
      post.b3[i] = as.numeric(as.character(files1J[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.2[,2]))==0)
    {
      post.b2[i] = as.numeric(as.character(files1J[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.1[,2]))==0)
    {
      post.b1[i] = as.numeric(as.character(files1J[[i]]$posterior[j]))
    }
  }
}


means.4 = array(0,10)
means.3 = array(0,10)
means.2 = array(0,10)
means.1 = array(0,10)
for(i in 1:10)
{
  if(i<10) k=i else k=0
  ids = which(1:100%%10==k)
  means.4[i]=mean(post.b4[ids])
  means.3[i]=mean(post.b3[ids])
  means.2[i]=mean(post.b2[ids])
  means.1[i]=mean(post.b1[ids])
}


post.b4 = array(0,100)
post.b3 = array(0,100)
post.b2 = array(0,100)
post.b1 = array(0,100)
for(i in 1:100)
{
  for(j in 1:length(files1G[[i]]$posterior))
  {
    
    expr=as.character(files1G[[i]]$tree[j])
    res=model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
    if(sum(abs(res[,2]-res.tpl.4[,2]))==0)
    {
      post.b4[i] = as.numeric(as.character(files1G[[i]]$posterior[j]))
    }else if(sum(abs(res[,2]-res.tpl.3[,2]))==0)
    {
      post.b3[i] = as.numeric(as.character(files1G[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.2[,2]))==0)
    {
      post.b2[i] = as.numeric(as.character(files1G[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.1[,2]))==0)
    {
      post.b1[i] = as.numeric(as.character(files1G[[i]]$posterior[j]))
    }
  }
}


meansa.4 = array(0,10)
meansa.3 = array(0,10)
meansa.2 = array(0,10)
meansa.1 = array(0,10)
for(i in 1:10)
{
  if(i<10) k=i else k=0
  ids = which(1:100%%10==k)
  meansa.4[i]=mean(post.b4[ids])
  meansa.3[i]=mean(post.b3[ids])
  meansa.2[i]=mean(post.b2[ids])
  meansa.1[i]=mean(post.b1[ids])
}
png(filename="plot14.png")
plot(ylim = c(0,1),x=c(1:10)*0.1*7,y = means.4,type = "l", col = 2, xlab =  expression('??'[4]),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=c(1:10)*0.1*7,y = meansa.4,type = "l", col = 4, xlab =   expression('??'[4]),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()
png(filename="plot13.png")
plot(ylim = c(0,1),x=c(1:10)*0.1*9,y = means.3,type = "l", col = 2, xlab =  expression('??'[3]),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=c(1:10)*0.1*9,y = meansa.3,type = "l", col = 4, xlab =   expression('??'[3]),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()
png(filename="plot12.png")
plot(ylim = c(0,1),x=c(1:10)*0.1*3.5,y = means.2,type = "l", col = 2, xlab =  expression('??'[2]),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=c(1:10)*0.1*3.5,y = meansa.2,type = "l", col = 4, xlab =   expression('??'[2]),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()
png(filename="plot11.png")
plot(ylim = c(0,1),x=c(1:10)*0.1*1.5,y = means.1,type = "l", col = 2, xlab =  expression('??'[1]),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=c(1:10)*0.1*1.5,y = meansa.1,type = "l", col = 4, xlab =   expression('??'[1]),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()

#make a plot for varying effect sizes for robust g prior for a fixed signal sizes and n = 1000 and reduced by sqrt(10) signal sizes and n = 10000
post.b4 = array(0,100)
post.b3 = array(0,100)
post.b2 = array(0,100)
post.b1 = array(0,100)
for(i in 1:100)
{
  for(j in 1:length(files1G[[i]]$posterior))
  {
    
    expr=as.character(files1G[[i]]$tree[j])
    res=model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
    if(sum(abs(res[,2]-res.tpl.4[,2]))==0)
    {
      post.b4[i] = as.numeric(as.character(files1G[[i]]$posterior[j]))
    }else if(sum(abs(res[,2]-res.tpl.3[,2]))==0)
    {
      post.b3[i] = as.numeric(as.character(files1G[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.2[,2]))==0)
    {
      post.b2[i] = as.numeric(as.character(files1G[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.1[,2]))==0)
    {
      post.b1[i] = as.numeric(as.character(files1G[[i]]$posterior[j]))
    }
  }
}


means.4 = array(0,10)
means.3 = array(0,10)
means.2 = array(0,10)
means.1 = array(0,10)
for(i in 1:10)
{
  if(i<10) k=i else k=0
  ids = which(1:100%%10==k)
  means.4[i]=mean(post.b4[ids])
  means.3[i]=mean(post.b3[ids])
  means.2[i]=mean(post.b2[ids])
  means.1[i]=mean(post.b1[ids])
}


post.b4 = array(0,100)
post.b3 = array(0,100)
post.b2 = array(0,100)
post.b1 = array(0,100)
for(i in 1:100)
{
  for(j in 1:length(files1GB[[i]]$posterior))
  {
    
    expr=as.character(files1GB[[i]]$tree[j])
    res=model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
    if(sum(abs(res[,2]-res.tpl.4[,2]))==0)
    {
      post.b4[i] = as.numeric(as.character(files1GB[[i]]$posterior[j]))
    }else if(sum(abs(res[,2]-res.tpl.3[,2]))==0)
    {
      post.b3[i] = as.numeric(as.character(files1GB[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.2[,2]))==0)
    {
      post.b2[i] = as.numeric(as.character(files1GB[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.1[,2]))==0)
    {
      post.b1[i] = as.numeric(as.character(files1GB[[i]]$posterior[j]))
    }
  }
}


meansa.4 = array(0,10)
meansa.3 = array(0,10)
meansa.2 = array(0,10)
meansa.1 = array(0,10)
for(i in 1:10)
{
  if(i<10) k=i else k=0
  ids = which(1:100%%10==k)
  meansa.4[i]=mean(post.b4[ids])
  meansa.3[i]=mean(post.b3[ids])
  meansa.2[i]=mean(post.b2[ids])
  meansa.1[i]=mean(post.b1[ids])
}
png(filename="plot14B.png")
plot(ylim = c(0,1),x=c(1:10)*0.1,y = means.4,type = "l", col = 4, xlab =  expression('K'),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=c(1:10)*0.1,y = meansa.4,type = "l", col = 4,lty = 2, xlab =   expression('K'),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()
png(filename="plot13B.png")
plot(ylim = c(0,1),x=c(1:10)*0.1,y = means.3,type = "l", col = 4, xlab =  expression('K'),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=c(1:10)*0.1,y = meansa.3,type = "l", col = 4,lty = 2, xlab =  expression('K'),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()
png(filename="plot12B.png")
plot(ylim = c(0,1),x=c(1:10)*0.1,y = means.2,type = "l", col = 4, xlab =  expression('K'),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=c(1:10)*0.1,y = meansa.2,type = "l", col = 4,lty = 2, xlab =  expression('K'),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()
png(filename="plot11B.png")
plot(ylim = c(0,1),x=c(1:10)*0.1,y = means.1,type = "l", col = 4, xlab =  expression('K'),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=c(1:10)*0.1,y = meansa.1,type = "l", col = 4,lty = 2, xlab =  expression('K'),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()


#make a plot for the varying n, fixed effect sizes for both Jeffrey's and robust g priors
post.b4 = array(0,100)
post.b3 = array(0,100)
post.b2 = array(0,100)
post.b1 = array(0,100)
for(i in 1:100)
{
  for(j in 1:length(files2J[[i]]$posterior))
  {
    
    expr=as.character(files2J[[i]]$tree[j])
    res=model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
    if(sum(abs(res[,2]-res.tpl.4[,2]))==0)
    {
      post.b4[i] = as.numeric(as.character(files2J[[i]]$posterior[j]))
    }else if(sum(abs(res[,2]-res.tpl.3[,2]))==0)
    {
      post.b3[i] = as.numeric(as.character(files2J[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.2[,2]))==0)
    {
      post.b2[i] = as.numeric(as.character(files2J[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.1[,2]))==0)
    {
      post.b1[i] = as.numeric(as.character(files2J[[i]]$posterior[j]))
    }
  }
}


means.4 = array(0,10)
means.3 = array(0,10)
means.2 = array(0,10)
means.1 = array(0,10)
for(i in 1:10)
{
  if(i<10) k=i else k=0
  ids = which(1:100%%10==k)
  means.4[i]=mean(post.b4[ids])
  means.3[i]=mean(post.b3[ids])
  means.2[i]=mean(post.b2[ids])
  means.1[i]=mean(post.b1[ids])
}



post.b4 = array(0,100)
post.b3 = array(0,100)
post.b2 = array(0,100)
post.b1 = array(0,100)
for(i in 1:100)
{
  for(j in 1:length(files2G[[i]]$posterior))
  {
    
    expr=as.character(files2G[[i]]$tree[j])
    res=model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
    if(sum(abs(res[,2]-res.tpl.4[,2]))==0)
    {
      post.b4[i] = as.numeric(as.character(files2G[[i]]$posterior[j]))
    }else if(sum(abs(res[,2]-res.tpl.3[,2]))==0)
    {
      post.b3[i] = as.numeric(as.character(files2G[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.2[,2]))==0)
    {
      post.b2[i] = as.numeric(as.character(files2G[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.1[,2]))==0)
    {
      post.b1[i] = as.numeric(as.character(files2G[[i]]$posterior[j]))
    }
  }
}


meansa.4 = array(0,10)
meansa.3 = array(0,10)
meansa.2 = array(0,10)
meansa.1 = array(0,10)
for(i in 1:10)
{
  if(i<10) k=i else k=0
  ids = which(1:100%%10==k)
  meansa.4[i]=mean(post.b4[ids])
  meansa.3[i]=mean(post.b3[ids])
  meansa.2[i]=mean(post.b2[ids])
  meansa.1[i]=mean(post.b1[ids])
}
png(filename="plot24.png")
plot(ylim = c(0,1),x=100*(1:10),y = means.4,type = "l", col = 2, xlab = expression('n'),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=100*(1:10),y = meansa.4,type = "l", col = 4, xlab =    expression('n'),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()
png(filename="plot23.png")
plot(ylim = c(0,1),x=100*(1:10),y = means.3,type = "l", col = 2, xlab =  expression('n'),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=100*(1:10),y = meansa.3,type = "l", col = 4, xlab =   expression('n'),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()
png(filename="plot22.png")
plot(ylim = c(0,1),x=100*(1:10),y = means.2,type = "l", col = 2, xlab =  expression('n'),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=100*(1:10),y = meansa.2,type = "l", col = 4, xlab =   expression('n'),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()
png(filename="plot21.png")
plot(ylim = c(0,1),x=100*(1:10),y = means.1,type = "l", col = 2, xlab =  expression('n'),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=100*(1:10),y = meansa.1,type = "l", col = 4, xlab =   expression('n'),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()

#make a plot for varying d, fixed effect sizes and n for both Jeffrey's and robust g priors
post.b4 = array(0,100)
post.b3 = array(0,100)
post.b2 = array(0,100)
post.b1 = array(0,100)
for(i in 1:100)
{
  for(j in 1:length(files3J[[i]]$posterior))
  {
    
    expr=as.character(files3J[[i]]$tree[j])
    res=model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
    if(sum(abs(res[,2]-res.tpl.4[,2]))==0)
    {
      post.b4[i] = as.numeric(as.character(files3J[[i]]$posterior[j]))
    }else if(sum(abs(res[,2]-res.tpl.3[,2]))==0)
    {
      post.b3[i] = as.numeric(as.character(files3J[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.2[,2]))==0)
    {
      post.b2[i] = as.numeric(as.character(files3J[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.1[,2]))==0)
    {
      post.b1[i] = as.numeric(as.character(files3J[[i]]$posterior[j]))
    }
  }
}


means.4 = array(0,10)
means.3 = array(0,10)
means.2 = array(0,10)
means.1 = array(0,10)
for(i in 1:10)
{
  if(i<10) k=i else k=0
  ids = which(1:100%%10==k)
  means.4[i]=mean(post.b4[ids])
  means.3[i]=mean(post.b3[ids])
  means.2[i]=mean(post.b2[ids])
  means.1[i]=mean(post.b1[ids])
}



post.b4 = array(0,100)
post.b3 = array(0,100)
post.b2 = array(0,100)
post.b1 = array(0,100)
for(i in 1:100)
{
  for(j in 1:length(files3G[[i]]$posterior))
  {
    
    expr=as.character(files3G[[i]]$tree[j])
    res=model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
    if(sum(abs(res[,2]-res.tpl.4[,2]))==0)
    {
      post.b4[i] = as.numeric(as.character(files3G[[i]]$posterior[j]))
    }else if(sum(abs(res[,2]-res.tpl.3[,2]))==0)
    {
      post.b3[i] = as.numeric(as.character(files3G[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.2[,2]))==0)
    {
      post.b2[i] = as.numeric(as.character(files3G[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.1[,2]))==0)
    {
      post.b1[i] = as.numeric(as.character(files3G[[i]]$posterior[j]))
    }
  }
}


meansa.4 = array(0,10)
meansa.3 = array(0,10)
meansa.2 = array(0,10)
meansa.1 = array(0,10)
for(i in 1:10)
{
  if(i<10) k=i else k=0
  ids = which(1:100%%10==k)
  meansa.4[i]=mean(post.b4[ids])
  meansa.3[i]=mean(post.b3[ids])
  meansa.2[i]=mean(post.b2[ids])
  meansa.1[i]=mean(post.b1[ids])
}
png(filename="plot34.png")
plot(ylim = c(0,1),x=15*(1:10),y = means.4,type = "l", col = 2, xlab = expression('d'),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=15*(1:10),y = meansa.4,type = "l", col = 4, xlab =    expression('d'),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()
png(filename="plot33.png")
plot(ylim = c(0,1),x=15*(1:10),y = means.3,type = "l", col = 2, xlab =  expression('d'),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=15*(1:10),y = meansa.3,type = "l", col = 4, xlab =   expression('d'),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()
png(filename="plot32.png")
plot(ylim = c(0,1),x=15*(1:10),y = means.2,type = "l", col = 2, xlab =  expression('d'),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=15*(1:10),y = meansa.2,type = "l", col = 4, xlab =   expression('d'),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()
png(filename="plot31.png")
plot(ylim = c(0,1),x=15*(1:10),y = means.1,type = "l", col = 2, xlab =  expression('d'),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=15*(1:10),y = meansa.1,type = "l", col = 4, xlab =   expression('d'),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()


#make a plot for varying r of the misspecifed and missing true leaf, with fixed effect sizes and n, for both Jeffrey's and robust g priors
post.b4 = array(0,100)
post.b3 = array(0,100)
post.b2 = array(0,100)
post.b1 = array(0,100)
for(i in 1:100)
{
  for(j in 1:length(files4J[[i]]$posterior))
  {
    
    expr=as.character(files4J[[i]]$tree[j])
    #print(expr)
    res=model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
    if(sum(abs(res[,2]-res.tpl.4[,2]))==0)
    {
      post.b4[i] = as.numeric(as.character(files4J[[i]]$posterior[j]))
    }else if(sum(abs(res[,2]-res.tpl.3[,2]))==0)
    {
      post.b3[i] = as.numeric(as.character(files4J[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.2[,2]))==0)
    {
      post.b2[i] = as.numeric(as.character(files4J[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.1[,2]))==0)
    {
      post.b1[i] = as.numeric(as.character(files4J[[i]]$posterior[j]))
    }
  }
}


means.4 = array(0,10)
means.3 = array(0,10)
means.2 = array(0,10)
means.1 = array(0,10)
for(i in 1:10)
{
  if(i<10) k=i else k=0
  ids = which(1:100%%10==k)
  means.4[i]=mean(post.b4[ids])
  means.3[i]=mean(post.b3[ids])
  means.2[i]=mean(post.b2[ids])
  means.1[i]=mean(post.b1[ids])
}
png(filename="plot1.png")
plot(ylim = c(0,1),x=c(1:10)*0.1,y = means.4,type = "l", col = 2,xlab =  expression('r'),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power",main="Power detecting a tree against signal strength")
dev.off()


post.b4 = array(0,100)
post.b3 = array(0,100)
post.b2 = array(0,100)
post.b1 = array(0,100)
for(i in 1:100)
{
  for(j in 1:length(files4G[[i]]$posterior))
  {
    
    expr=as.character(files4G[[i]]$tree[j])
    #print(expr)
    res=model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
    if(sum(abs(res[,2]-res.tpl.4[,2]))==0)
    {
      post.b4[i] = as.numeric(as.character(files4G[[i]]$posterior[j]))
    }else if(sum(abs(res[,2]-res.tpl.3[,2]))==0)
    {
      post.b3[i] = as.numeric(as.character(files4G[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.2[,2]))==0)
    {
      post.b2[i] = as.numeric(as.character(files4G[[i]]$posterior[j]))
    }
    else if(sum(abs(res[,2]-res.tpl.1[,2]))==0)
    {
      post.b1[i] = as.numeric(as.character(files4G[[i]]$posterior[j]))
    }
  }
}


meansa.4 = array(0,10)
meansa.3 = array(0,10)
meansa.2 = array(0,10)
meansa.1 = array(0,10)
for(i in 1:10)
{
  if(i<10) k=i else k=0
  ids = which(1:100%%10==k)
  meansa.4[i]=mean(post.b4[ids])
  meansa.3[i]=mean(post.b3[ids])
  meansa.2[i]=mean(post.b2[ids])
  meansa.1[i]=mean(post.b1[ids])
}
png(filename="plot44.png")
plot(ylim = c(0,1),x=0.1*(1:10),y = means.4,type = "l", col = 2, xlab = expression('r'),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=0.1*(1:10),y = meansa.4,type = "l", col = 4, xlab =    expression('r'),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()
png(filename="plot43.png")
plot(ylim = c(0,1),x=0.1*(1:10),y = means.3,type = "l", col = 2, xlab =  expression('r'),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=0.1*(1:10),y = meansa.3,type = "l", col = 4, xlab =   expression('r'),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()
png(filename="plot42.png")
plot(ylim = c(0,1),x=0.1*(1:10),y = means.2,type = "l", col = 2, xlab = expression('r'),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=0.1*(1:10),y = meansa.2,type = "l", col = 4, xlab =   expression('r'),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()
png(filename="plot41.png")
plot(ylim = c(0,1),x=0.1*(1:10),y = means.1,type = "l", col = 2, xlab =  expression('r'),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=0.1*(1:10),y = meansa.1,type = "l", col = 4, xlab =   expression('r'),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()

