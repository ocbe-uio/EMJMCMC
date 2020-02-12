#read the most recent stable version of the package
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")
update.packages("https://github.com/aliaksah/EMJMCMC2016/blob/master/EMJMCMC_1.4.2_R_x86_64-pc-linux-gnu.tar.gz?raw=true", 
                 repos = NULL, type="source")
# load the EMJMCMC package
library(EMJMCMC)
set.seed(040590)
#make sure that you are using Mac Os or Linux (mclapply is currently not supported for Windows unless some mclapply hack function for Windows is preloaded in your R session)



## Construct a binary correlation matrix for M = 50 cariables
M = 50
m = clusterGeneration::rcorrmatrix(M,alphad=2.5) 
#print the highest 10 correlations in the data
print(unique(sort(abs(m),decreasing = T))[1:10])
#print the lowest 10 correlations in the data
print(unique(sort(abs(m),decreasing = F))[1:10])
#simulate 1000 binary variables with a given correlation's structure
X = bindata::rmvbin(1000, margprob = rep(0.5,M), bincorr = m)

print(unique(sort(abs(cor(X)),decreasing = T))[1:10])

melted_cormat = reshape2::melt(cor(X))
ggplot2::ggplot(data = melted_cormat, 
ggplot2::aes(x=Var1, y=Var2, fill=value)) + 
  ggplot2::geom_tile() +  
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.title.y =  ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank())

#simulate Gaussian responses from a model with two-way interactions
#which is considered in S.4 of the paper
X = data.frame(X)
Y=rnorm(n = 1000,mean = 1+0.7*(X$X1*X$X4) + 0.8896846*(X$X8*X$X11)+1.434573*(X$X5*X$X9),sd = 1)
X$Y=Y

#specify the initial formula
formula1 = as.formula(paste(colnames(X)[M+1],"~ 1 +",paste0(colnames(X)[-c(M+1)],collapse = "+")))
df = as.data.frame(X)

#run the inference with robust g prior
res4G = LogicRegr(formula = formula1,data = df ,family = "Gaussian",prior = "G",report.level = 0.5,d = 15,cmax = 2,kmax = 15,p.and = 0.9,p.not = 0.01,p.surv = 0.2,ncores = 32)
print(res4G$feat.stat)

#run the inference with Jeffrey's prior
res4J = LogicRegr(formula = formula1,data = data.example,family = "Gaussian",prior = "J",report.level = 0.5,d = 15,cmax = 2,kmax = 15,p.and = 0.9,p.not = 0.01,p.surv = 0.2,ncores = 32)
print(res4J$feat.stat)


#now use qtl package to generate data in a more realistic scenario

library(qtl)
set.seed(040590)  #you can change this of course

# simulate 5 autosomes, each with 10 equally spaced markers 
# with different chromosomal lengths 
map = sim.map(c(100,90,80,60,40), M/5, include.x=FALSE, eq.spacing=TRUE)
plotMap(map)

# The longer the chromosomal length the less correlated are markers
# (at least theoretically)

# Now simulate data from a backcross design using the sim.cross function
n.ind = 1000   #sample size

simbc = sim.cross(map, type="bc", n.ind=n.ind,
                    model=rbind(c(1,45,1), c(5,20,1), c(5,50,1)))

# The relevant genotype data is in the structure geno 
str(simbc$geno)   #it is a bit complicated

# Get an X matrix for each chromosome
X.list = vector("list", 5)
for (chr in 1:5){
  X.list[[chr]] = pull.geno(simbc, chr)
}

#Check the correlations between markers within the same chromosome
lapply(X.list, cor)


#Create one large X matrix which you can the use to make your own 
# simulations of Y data with a logic regression model
X = cbind(X.list[[1]],X.list[[2]],X.list[[3]],X.list[[4]],X.list[[5]])-1
#permute elements of X

#plot the correlation's structure 
melted_cormat = reshape2::melt(cor(as.matrix(X)))
ggplot2::ggplot(data = melted_cormat, 
                ggplot2::aes(x=Var1, y=Var2, fill=value)) + 
  ggplot2::geom_tile() +  
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.title.y =  ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank())
X = X[,sample.int(n = M,size = M,replace = F)]


X = data.frame(X)
Y=rnorm(n = 1000,mean = 1+0.7*(X[,1]*X[,4]) + 0.8896846*(X[,8]*X[,11])+1.434573*(X[,5]*X[,9]),sd = 1)
X$Y=Y

print(paste0("L1 = ", paste(names(X)[c(1,4)],collapse = ", ",sep = "")))
print(paste0("L2 = ", paste(names(X)[c(8,11)],collapse = ", ",sep = "")))
print(paste0("L3 = ", paste(names(X)[c(5,9)],collapse = ", ",sep = "")))

#specify the initial formula
formula1 = as.formula(paste(colnames(X)[M+1],"~ 1 +",paste0(colnames(X)[-c(M+1)],collapse = "+")))
data.example = as.data.frame(X)

#run the inference with robust g prior
res4G = LogicRegr(formula = formula1,data = data.example,family = "Gaussian",prior = "G",report.level = 0.5,d = 15,cmax = 2,kmax = 15,p.and = 0.9,p.not = 0.01,p.surv = 0.2,ncores = 32)
print(res4G$feat.stat)
#run the inference with Jeffrey's prior
res4J = LogicRegr(formula = formula1,data = data.example,family = "Gaussian",prior = "J",report.level = 0.5,d = 15,cmax = 2,kmax = 15,p.and = 0.9,p.not = 0.01,p.surv = 0.2,ncores = 32)
print(res4J$feat.stat)

