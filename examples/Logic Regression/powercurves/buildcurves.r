library("stringi")

#setwd("/mn/sarpanitu/ansatte-u2/aliaksah/abeldata1/logic-g-prior/powercurves")
#setwd("/mn/sarpanitu/ansatte-u2/aliaksah/abeldata1/logic-g-prior/copcurves")

temp = list.files(pattern="*_T1J_*")
temp = temp[order(x = unlist(mclapply(FUN = function(str){
  as.integer(stri_split_fixed(str = stri_split_fixed(str = str,pattern = c("_"))[[1]][3],pattern = ".")[[1]][1]
)},X = temp)),decreasing = F)][-31]
files1J = mclapply(FUN = read.csv,X = temp)

temp = list.files(pattern="*_T1G_*")
temp = temp[order(x = unlist(mclapply(FUN = function(str){
  as.integer(stri_split_fixed(str = stri_split_fixed(str = str,pattern = c("_"))[[1]][3],pattern = ".")[[1]][1]
  )},X = temp)),decreasing = F)]
files1G = mclapply(FUN = read.csv,X = temp)

temp = list.files(pattern="*_T2J_*")
temp = temp[order(x = unlist(mclapply(FUN = function(str){
  as.integer(stri_split_fixed(str = stri_split_fixed(str = str,pattern = c("_"))[[1]][3],pattern = ".")[[1]][1]
  )},X = temp)),decreasing = F)][-c(30:36)]
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

X<- as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))

res.tpl<-model.matrix(data=X,object = as.formula(paste0("V1~","I(V4*V17*V30*V10)")))

post.b5 = array(0,100)
for(i in 1:100)
{
  for(j in 1:length(files1J[[i]]$posterior))
  {

      expr<-as.character(files1J[[i]]$tree[j])
      #print(expr)
      res<-model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
      if(sum(abs(res[,2]-res.tpl[,2]))==0)
      {
        post.b5[i] = as.numeric(as.character(files1J[[i]]$posterior[j]))
      }
  }
}
png(filename="plot1a.png")
thickness = c(20,15,10,5)
plot(ylim = c(0,1),x=c(1:10),y = post.b5[1:10],type = "l", col = 2, lwd=thickness[1],xlab = expression('β'[3]),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Posterior",main="Posterior of L_3 against signal strength (Jeffrey's prior)")
for(i in 2:10)
  lines(x=c(1:10),y = post.b5[c((10*(i-1)+1):(10*(i-1)+10))],type = "l", col = i+1, lwd=thickness[i])
dev.off()

means = array(0,10)
for(i in 1:10)
{
  if(i<10) k=i else k=0
  ids = which(1:100%%10==k)
  means[i]=mean(post.b5[ids])
}
png(filename="plot1.png")
plot(ylim = c(0,1),x=c(1:10),y = means,type = "l", col = 2,xlab =  expression('β'[3]),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power",main="Power detecting of L_3 against signal strength (Jeffrey's prior)")
dev.off()
meansa=means
post.b5 = array(0,100)
for(i in 1:100)
{
  for(j in 1:length(files1G[[i]]$posterior))
  {

    expr<-as.character(files1G[[i]]$tree[j])
    #print(expr)
    res<-model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
    if(sum(abs(res[,2]-res.tpl[,2]))==0)
    {
      post.b5[i] = as.numeric(as.character(files1G[[i]]$posterior[j]))
    }
  }
}
png(filename="plot2a.png")
thickness = c(20,15,10,5)
plot(ylim = c(0,1),x=c(1:10),y = post.b5[1:10],type = "l", col = 2, lwd=thickness[1],xlab =  expression('β'[3]),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Posterior",main="Posterior of L_3 against signal strength (Robust g-prior)")
for(i in 2:4)
  lines(x=c(1:10),y = post.b5[c((10*(i-1)+1):(10*(i-1)+10))],type = "l", col = i+1, lwd=thickness[i])
dev.off()

means = array(0,10)
for(i in 1:10)
{
  if(i<10) k=i else k=0
  ids = which(1:100%%10==k)
  means[i]=mean(post.b5[ids])
}
png(filename="plot2.png")
plot(ylim = c(0,1),x=c(1:10),y = meansa,type = "l", col = 2,xlab =  expression('β'[4]),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=c(1:10),y = means,type = "l", col = 4,xlab =  expression('β'[4]),cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()


post.b5 = array(0,100)
for(i in 1:100)
{
  for(j in 1:length(files2J[[i]]$posterior))
  {

    expr<-as.character(files2J[[i]]$tree[j])
    #print(expr)
    res<-model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
    if(sum(abs(res[,2]-res.tpl[,2]))==0)
    {
      post.b5[i] = as.numeric(as.character(files2J[[i]]$posterior[j]))
    }
  }
}

png(filename="plot3a.png")
thickness = c(20,15,10,5)
plot(ylim = c(0,1),x=c(1:10)*100,y = post.b5[1:10],type = "l", col = 2, lwd=thickness[1],xlab = "n",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Posterior",main="Posterior of L_3 against sample size (Jeffrey's prior)")
for(i in 2:4)
  lines(x=c(1:10)*100,y = post.b5[c((10*(i-1)+1):(10*(i-1)+10))],type = "l", col = i+1, lwd=thickness[i])
dev.off()


means = array(0,10)
for(i in 1:10)
{
  if(i<10) k=i else k=0
  ids = which(1:100%%10==k)
  means[i]=mean(post.b5[ids])
}
png(filename="plot3.png")
plot(ylim = c(0,1),x=c(1:10)*100,y = means,type = "l", col = 2,xlab = "n",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power",main="Power detecting of L_3 against sample size (Jeffrey's prior)")
dev.off()

meansa=means

post.b5 = array(0,100)
for(i in 1:100)
{
  for(j in 1:length(files2G[[i]]$posterior))
  {

    expr<-as.character(files2G[[i]]$tree[j])
    #print(expr)
    res<-model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
    if(sum(abs(res[,2]-res.tpl[,2]))==0)
    {
      post.b5[i] = as.numeric(as.character(files2G[[i]]$posterior[j]))
    }
  }
}

png(filename="plot4a.png")
thickness = c(20,15,10,5)
plot(ylim = c(0,1),x=c(1:10)*100,y = post.b5[1:10],type = "l", col = 2, lwd=thickness[1],xlab = "n",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Posterior",main="Posterior of L_3 against sample size (Robust g-prior)")
for(i in 2:4)
  lines(x=c(1:10)*100,y = post.b5[c((10*(i-1)+1):(10*(i-1)+10))],type = "l", col = i+1, lwd=thickness[i])
dev.off()


means = array(0,10)
for(i in 1:10)
{
  if(i<10) k=i else k=0
  ids = which(1:100%%10==k)
  means[i]=mean(post.b5[ids])
}
png(filename="plot4.png")
plot(ylim = c(0,1),x=c(1:10)*100,y = meansa,type = "l", col = 2,xlab = "n",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=c(1:10)*100,y = means,type = "l", col = 4,xlab = "n",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()


post.b5 = array(0,100)
for(i in 1:100)
{
  for(j in 1:length(files3J[[i]]$posterior))
  {

    expr<-as.character(files3J[[i]]$tree[j])
    #print(expr)
    res<-model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
    if(sum(abs(res[,2]-res.tpl[,2]))==0)
    {
      post.b5[i] = as.numeric(as.character(files3J[[i]]$posterior[j]))
    }
  }
}

png(filename="plot5a.png")
thickness = c(20,15,10,5)
plot(ylim = c(0,1),x=c(1:10)*15,y = post.b5[1:10],type = "l", col = 2, lwd=thickness[1],cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,xlab = "d",ylab="Posterior",main="Posterior of L_3 against pool size (Jeffrey's prior)")
for(i in 2:4)
  lines(x=c(1:10)*15,y = post.b5[c((10*(i-1)+1):(10*(i-1)+10))],type = "l", col = i+1, lwd=thickness[i])
dev.off()


means = array(0,10)
for(i in 1:10)
{
  if(i<10) k=i else k=0
  ids = which(1:100%%10==k)
  means[i]=mean(post.b5[ids])
}
png(filename="plot5.png")
plot(ylim = c(0,1),x=c(1:10)*15,y = means,type = "l", col = 2,xlab = "d",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power",main="Power detecting of L_3 against pool size (Jeffrey's prior)")
dev.off()

meansa=means

post.b5 = array(0,100)
for(i in 1:100)
{
  for(j in 1:length(files3G[[i]]$posterior))
  {

    expr<-as.character(files3G[[i]]$tree[j])
    #print(expr)
    res<-model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
    if(sum(abs(res[,2]-res.tpl[,2]))==0)
    {
      post.b5[i] = as.numeric(as.character(files3G[[i]]$posterior[j]))
    }
  }
}

png(filename="plot6a.png")
thickness = c(20,15,10,5)
plot(ylim = c(0,1),x=c(1:10)*15,y = post.b5[1:10],type = "l",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col = 2, lwd=thickness[1],xlab = "n",ylab="Posterior",main="Posterior of L_3 against pool size (Robust g-prior)")
for(i in 2:4)
  lines(x=c(1:10)*15,y = post.b5[c((10*(i-1)+1):(10*(i-1)+10))],type = "l", col = i+1, lwd=thickness[i])
dev.off()

means = array(0,10)
for(i in 1:10)
{
  if(i<10) k=i else k=0
  ids = which(1:100%%10==k)
  means[i]=mean(post.b5[ids])
}
png(filename="plot8.png")
plot(ylim = c(0,1),x=c(1:10)*15,y = meansa,type = "l", col = 2,xlab = "d",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
lines(ylim = c(0,1),x=c(1:10)*15,y = means,type = "l", col = 4,xlab = "d",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="Power")
dev.off()

