# example of generating the addressed data set
M=100
resa<-array(data = 0,dim = c(16,M*3))
NM= 1000
post.popul <- array(0,M)
max.popul <- array(0,M)
set.seed(040590)
X1<- as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))
Y1<-rnorm(n = 1000,mean = 1+0.7*(X1$V1*X1$V4) + 0.8896846*(X1$V8*X1$V11)+1.434573*(X1$V5*X1$V9),sd = 1)
X1$Y1<-Y1
write.csv(x =X1,file = "X1.csv")

