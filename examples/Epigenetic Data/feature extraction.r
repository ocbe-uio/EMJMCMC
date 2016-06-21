#read and preprocess the data
workdir<-"/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/Untitled Folder/"


X <- read.table(paste(workdir,"GSM1085222_mC_calls_Col_0.tsv",sep = "",collapse = ""), header=T)
dim(X)
cc<-diff(X$pos,lag = 1)
cc<-c(0,cc)
X$base_dist <- cc
X$base_dist[which(X$base_dist<0)]<-0
X$CHG <- as.integer((substr(as.character(X$mc_class), 3, 3)=="G" & substr(as.character(X$mc_class), 2, 2)!="G"))
X$CG <- as.integer(substr(as.character(X$mc_class), 2, 2)=="G")
X$CHH <- as.integer(substr(as.character(X$mc_class), 2, 2)!="G" & substr(as.character(X$mc_class), 3, 3)!="G")
X$DT1 <- as.integer(X$base_dist==1)
X$DT2 <- as.integer(X$base_dist==2)
X$DT3 <- as.integer(X$base_dist==3)
X$DT4 <- as.integer(X$base_dist==4)
X$DT5 <- as.integer(X$base_dist==5)
X$DT6_20 <- as.integer(X$base_dist>=6 & X$base_dist<20)
X$DT20_inf <- as.integer(X$base_dist>=20)
X$unmethylated_bases <-X$total_bases - X$methylated_bases

#read and preprocess the annotation
GGroup <- read.table(paste(workdir,"GeneGroups.csv",sep = "",collapse = ""), header=T,sep = '\t',blank.lines.skip = TRUE,fill = TRUE)

Y <- read.table(paste(workdir,"Genes.csv",sep = "",collapse = ""), header=T,sep = '\t')


Express<- read.table(paste(workdir,"GSM1086840_Col_0_expression2.tsv",sep = "",collapse = ""), header=T,sep = '\t')

for(i in 1:nrow(GGroup))
{
  #print(i)
  if(!(as.character(GGroup$TYPE[[i]]) %in% colnames(X)))
    X[[as.character(GGroup$TYPE[[i]])]]<-0
  set<-which(gsub(" ", "", as.character(Y$Locus.tag), fixed = TRUE) == gsub(" ", "", as.character(GGroup$GENE_ID[[i]]), fixed = TRUE))
  set<-c(set,((which(toupper(gsub(" ", "", as.character(Y$Locus.tag), fixed = TRUE)) == gsub(" ", "", as.character(GGroup$GENE_ID[[i]]), fixed = TRUE)))))
  set<-unique(set)
  if((  length(set)>1 ))
  {

    print(paste(as.character(GGroup$GENE_ID[[i]])," is found with protein products:"))
    print(as.character(Y$Protein.product[set]))
    for(j in set)
    {
      #print(Y$Start[j])
      set4<-which(as.integer(X$pos)>=as.integer(Y$Start[j]) & as.integer(X$pos)<=as.integer(Y$Stop[j]) & X$chrom==as.integer(Y$Replicon.Name[j]))
      if(length(set4)>0)
      {
        print(paste(as.character(Y$Protein.product[set]),"found in epi data with ",length(set4)," base pairs"))
        X[[as.character(GGroup$TYPE[[i]])]][set4]<-1
      }
    }

  }
  else
  {
    print(paste(as.character(GGroup$GENE_ID[[i]])," is not found in EPI dataset"))
  }
}
X$none <- 1
X$none <- X$none - (X$Mα + X$Mβ + X$Mβ + X$Mγ + X$Mδ + X$MIKC)

# now mark all promoters, coding regions and transcriptions
X$coding<-0
X$promoter<-0
X$aftron<-0
#define length of the promoter regions
lprom_length<-250
laftr_length<-250
for(i in 1:nrow(Y))
{
  if(i%%10==0)
    print(paste(i,"genes proceed"))
  set4<-which(as.integer(X$pos)>=as.integer(Y$Start[i]) & as.integer(X$pos)<=as.integer(Y$Stop[i]) & X$chrom==as.integer(Y$Replicon.Name[i]))
  set5<-which(X$chrom==as.integer(Y$Replicon.Name[i]))
  if(length(set4)>0)
  {

  ls4<-set4[length(set4)]
  fs4<-set4[1]
  set6<-(fs4-lprom_length):(fs4-1)
  set6<-intersect(set6, set5)
  set81<-NULL
  set82<-NULL
  if(i!=1){
    set8<-intersect(set6, set7)
    if(length(set8)>0){
      set81<-set8[1:round(length(set8)/2)]
      set82<-set8[(round(length(set8)/2)+1):length(set8)]
      if(length(set82)>0)
        X$aftron[set82]<-0
    }
  }
  set7<-(ls4+1):(ls4+laftr_length)
  set7<-intersect(set7, set5)
  #devide the sets in a clever way

    X$coding[set4]<-1
    X$promoter[set6]<-1
    if(length(set81)>0)
      X$promoter[set81]<-0
    X$aftron[set7]<-1
    X$promoter[set4]<-0
    X$aftron[set4]<-0
  }
}

ss<-which(X$promoter + X$aftron == 2)
length(ss)
X$promoter[ss]<-0
X$aftron[ss]<-0
ss1<-which(X$promoter==1)
ss2<-which(X$aftron==1)
ss3<-which(X$coding==1)
length(ss1)
length(ss2)
length(ss3)

classes<-which(X$coding +  X$aftron + X$promoter > 0 & X$strand == "+")
length(classes)


X$express<-0
for(i in 1:nrow(Express))
{
    if(i%%10==0)
      print(paste(i,"gene expressions proceed"))
    set4<-which(as.integer(X$pos)>=as.integer(Express$start[i]) & as.integer(X$pos)<= as.integer(Express$end[i]) & X$chrom==as.integer(Express$chrom[i]))
    if(length(set4)>0)
    {
      X$express[set4]<-Express$fpkm[i]
    }
}

data<-X[classes,]
View(data[which(data$promoter ==1),])
data$P <-data$methylated_bases/data$total_bases
data$pos <-1:dim(data)[1]
closeAllConnections()

# set up the parameters of the simulation or optimization
M<-5
size<-1
fparam <- colnames(X)[c(9:10,12:17)]
fobserved <- colnames(X)[c(5,19)]

bittodec(c(0,0,0,0,0,0,0,0,0,0,0))
#
# xxx<-statistics1

statistics1 <- big.matrix(nrow = 2 ^(length(fparam)+1), ncol = 14 + length(fparam),init = NA, type = "double")
statistics <- describe(statistics1)

#stat<-big.matrix(nrow = 2 ^(length(fparam)+1), ncol = 6,init = 0, type = "double")

# statistics1[1:2048,9:13]=0
# statistics1[1:2048,3]=0
# #describe(statistics)
# ssastat1 <- statistics1
# ssastat  <- statistics
#
# allstat1 <- statistics1
# allstat  <- statistics
#
# sstat1 <- statistics1
# sstat  <- statistics
#
# statistics1 <- big.matrix(nrow = 2 ^(14), ncol = 13,init = rnorm(n = 2 ^(14)*13,mean = 1000,sd = 100), type = "double")
# statistics <- describe(statistics1)
# statistics1[,1]<-rnorm(n = 2 ^(14),mean = -1000,sd = 30)
# statistics1[,2]<-rnorm(n = 2 ^(14),mean =  -statistics1[,1],sd = 10)
# statistics1[,3]<-as.integer(runif(n = 2 ^(14),min = 0, max = 100*statistics1[,2]))
# statistics1[,4:13]<-as.integer(runif(n = 2 ^(14)*10,min = 0, max = statistics1[,3]))
# statistics1[,1]<-statistics1[,1]/100
# statistics1<-read.big.matrix(filename = "/mn/anatu/ansatte-u3/aliaksah/Desktop/publication 1/plots/bigmatrix 5 .csv")
# statistics <- describe(statistics1)
# statistics1[1,]<-NA

partype <- "PSOCK"
max.N <-5
min.N <-2
ar.p.max<-10
binar<-c(0,1,1,2,2,2,2,2,2)
#binar<-c(0,1,1,2,2,2,2,2,2,3,3,3,3,3)
a.P.prior <- as.integer(runif(n = length(fparam)+1,min = 1,max=1))
paral.type <-"SOCK"
max.cpu.per.run <- 5
deltadata<-1000
opt<-array(data = "",dim = M+1)
opt[1]<-paste(c("MOD_ID","CHROM","LBASE","UBASE","CONST",fparam,"WAIC","MLIK","TIME"),collapse = " ")
vect<-list()
# do the simulations or optimization
# View(X[X$chrom == 1,][6399210:6400710,])
done<-TRUE
iii<-8
withRestarts(tryCatch({

  chrom_id<-round(runif(n = 1,min = 1, max = 5))
  chrom_id<-1
  uub<-dim(X[X$chrom == chrom_id,])[1]
  vect<-list()
  for(iii in 3:M)
  {

  lb<-2073472#runif(n = 1, min = 1, max = uub - deltadata)
  ub<-lb+1500#runif(n = 1, min = lb, max = lb + deltadata)
  tmps<-paste(iii,chrom_id,lb,ub,collapse = ";")
  data<-X[X$chrom == chrom_id,][(lb-50000):(ub-50000),1:25]
  data$P <-data$methylated_bases/data$total_bases
  data$pos <-1:dim(data)[1]

  LocImproveParam <- list(t.min = 0.01, t.init =1000, dt = 6, M = 10, max.cpu  = 5)
  vect[[iii]]<-list(max.cpu = 1,objective = 1,p.prior = runif(n = length(fparam)+1, min = 0.5,max = 0.5),LocImproveParam = LocImproveParam,  changeble = array(data = 0,dim =(length(fparam)+1) ), data = data,switch.type=3,n.size=5, LocImprove = 2, fobserved = fobserved,min.N = 3,max.N=7,isobsbinary = binar,fparam = fparam,p.add = 0.2,ar.p = 1,maxit = 3,trit = 20,eps = 0.001,paral.type = paral.type, range = 0.3, burnin = 1 ,seed = runif(n=1, min = 1,max = 10000), RJMCMCvsWAICMC =0 , beta.mu.prior = 0 ,beta.tau.prior = 0.001, tau.mu.prior = 1, tau.prec.prior = 0.001, psi.mu.prior = 0,psi.prec.prior = 0.001, dic.t = FALSE, mlik.t = TRUE, waic.t = TRUE)
  #print(paste("finish preparing MCMC parameters", i))
  statistics1 <- big.matrix(nrow = 2 ^(length(fparam)+1), ncol = 14 + length(fparam),init = NA, type = "double")
  statistics <- describe(statistics1)
#   statistics1[,1]<-rnorm(n = 2048,mean = -300,sd = 30)
#   statistics1[,2]<- -2*statistics1[,1] + rnorm(n = 2048,mean = 0,sd = 100)
#  statistics1[,3:13]<- 0
#   statistics1[,4:8]<- 1
  vect[[iii]]$changeble <- array(data = 0,dim =(length(fparam)+1) )
  vect[[iii]]$burnin = 100
  vect[[iii]]$min.N = 1
  vect[[iii]]$max.N = 6
  vect[[iii]]$n.size= 10
  vect[[iii]]$p.add <- 0.5
  vect[[iii]]$max.cpu <- 3
  vect[[iii]]$switch.type <- 1
  vect[[iii]]$LocImprove <- c(50,0,0,0,100) # distribution of the search methods to be addressed not need to be normilized
  vect[[iii]]$steps<-1
  vect[[iii]]$min.N.randomize <- 0
  vect[[iii]]$max.N.randomize <- 2
  vect[[iii]]$type.randomize  <- 3
  vect[[iii]]$maxit = 40000
  vect[[iii]]$aa = 0.9
  vect[[iii]]$cc = 0
  vect[[iii]]$seed = runif(n = 1,min = 0,max = 100000)
  LocImproveParam <- list(t.min = 0.01, t.init =100,dt = 3, M = 2, max.cpu  = 4)
  vect[[iii]]$LocImproveParam<-LocImproveParam
  invisible(res<-mcmc_wrap(vect[[iii]]))

  opt[iii+1]<-paste(tmps,paste(res$vars,collapse = " "),res$waic, res$mlik, res$time)
  write.big.matrix(x = statistics1,filename = paste("/mn/anatu/ansatte-u3/aliaksah/Desktop/publication 1/plots/bigmatrix",iii,".csv"))

  ids=which(statistics1[,1]==-100000)
  ids<-c(ids,which(statistics1[,1]==-Inf))
  ids<-c(ids,which(statistics1[,2]==Inf))
#   length(which(statistics1[,1]<0))
#   length(which(statistics1[,2]>0))
#   length(which(statistics1[,3]>0))
#   length(which(statistics1[,9]>0))
#   length(which(statistics1[,10]>0))
#   length(which(statistics1[,11]>0))
#   length(which(statistics1[,12]>0))
#   length(which(statistics1[,13]>0))

  if(length(ids)!=0)
  {
    mlik.lim<-c(min(statistics1[-ids,1],na.rm = TRUE),max(statistics1[-ids,1],na.rm = TRUE))
    waic.lim<-c(min(statistics1[-ids,2],na.rm = TRUE),max(statistics1[-ids,2],na.rm = TRUE))
  }else
  {
    mlik.lim<-c(min(statistics1[,1],na.rm = TRUE),max(statistics1[,1],na.rm = TRUE))
    waic.lim<-c(min(statistics1[,2],na.rm = TRUE),max(statistics1[,2],na.rm = TRUE))
  }

  norm<-0.5*sqrt(sum(statistics1[,3],na.rm = TRUE))

  jpeg(file=paste(workdir,iii,"plot0.jpg",sep = ""))
  plot(xlab = "posterior probability of a visit found", ylab="visit type",c(0.1,0.2,0.3,0.4,0.5,0.6,0.7),c(7,7,7,-1,-1,-1,-1),ylim = c(0,9), xlim=c(0,1),pch=19, col = 7,cex= c(2,2,2,0,0,0,0))
  points(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7),c(6,6,6,-1,-1,-1,-1),pch=8,  col = 5,cex=  c(1,2,3,0,0,0,0))
  points(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7),c(1,1,1,-1,-1,-1,-1),pch=2,  col = 2,cex= c(1,2,3,0,0,0,0))
  points(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7),c(2,2,2,-1,-1,-1,-1),pch=3,  col = 3,cex= c(1,2,3,0,0,0,0))
  points(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7),c(3,3,3,-1,-1,-1,-1),pch=4,  col = 4,cex= c(1,2,3,0,0,0,0))
  points(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7),c(4,4,4,-1,-1,-1,-1),pch=6,  col = 6,cex=c(1,2,3,0,0,0,0))
  points(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7),c(5,5,5,-1,-1,-1,-1),pch=1,  col = 1,cex= c(1,2,3,0,0,0,0))
  points(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7),c(8,8,8,-1,-1,-1,-1),pch=19,  col = 2,cex= c(2,2,2,0,0,0,0))
  legend(x = 0.35,y = 1.5,legend = "- SA 1 local optimization accepted",bty="n" )
  legend(x = 0.35,y = 2.5,legend = "- MTMCMC local optimization accepted",bty="n" )
  legend(x = 0.35,y = 3.5,legend = "- GREEDY local optimization accepted",bty="n" )
  legend(x = 0.35,y = 4.5,legend = "- SA 2 local optimization accepted",bty="n" )
  legend(x = 0.35,y = 5.5,legend = "- NO local optimization accepted",bty="n" )
  legend(x = 0.35,y = 6.5,legend = "- Totally accepted",bty="n" )
  legend(x = 0.35,y = 7.5,legend = "- Totally explored",bty="n" )
  legend(x = 0.35,y = 8.5,legend = "- Bayes formula based posterior",bty="n" )
  dev.off()

  jpeg(file=paste(workdir,iii,"plot1.jpg",sep = ""))
  plot(ylim = mlik.lim, xlab = "model_id", ylab="MLIK" , statistics1[,1],pch=19, col = 7,cex= 1*((statistics1[,4]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
  points(statistics1[,1],pch=8,  col = ifelse(statistics1[,3]>0,5,0),cex= ifelse(statistics1[,3]>0,statistics1[,3]/norm+1,0))
  points(statistics1[,1],pch=2,  col = ifelse(statistics1[,9]>0,2,0),cex= ifelse(statistics1[,9]>0,statistics1[,9]/norm+1,0))
  points(statistics1[,1],pch=3,  col = ifelse(statistics1[,10]>0,3,0),cex= ifelse(statistics1[,10]>0,statistics1[,10]/norm+1,0))
  points(statistics1[,1],pch=4,  col = ifelse(statistics1[,11]>0,4,0),cex= ifelse(statistics1[,11]>0,statistics1[,11]/norm+1,0))
  points(statistics1[,1],pch=6,  col = ifelse(statistics1[,12]>0,6,0),cex= ifelse(statistics1[,12]>0,statistics1[,12]/norm+1,0))
  points(statistics1[,1],pch=1,  col = ifelse(statistics1[,13]>0,1,0),cex= ifelse(statistics1[,13]>0,statistics1[,13]/norm+1,0))
  dev.off()

  jpeg(file=paste(workdir,iii,"plot2.jpg",sep = ""))
  plot(ylim = waic.lim, xlab = "model_id", ylab="WAIC" ,statistics1[,2],pch=19, col = 7,cex= 1*((statistics1[,4]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
  points(statistics1[,2],pch=8,  col = ifelse(statistics1[,3]>0,5,0),cex= ifelse(statistics1[,3]>0,statistics1[,3]/norm+1,0))
  points(statistics1[,2],pch=2,  col = ifelse(statistics1[,9]>0,2,0),cex= ifelse(statistics1[,9]>0,statistics1[,9]/norm+1,0))
  points(statistics1[,2],pch=3,  col = ifelse(statistics1[,10]>0,3,0),cex= ifelse(statistics1[,10]>0,statistics1[,10]/norm+1,0))
  points(statistics1[,2],pch=4,  col = ifelse(statistics1[,11]>0,4,0),cex= ifelse(statistics1[,11]>0,statistics1[,11]/norm+1,0))
  points(statistics1[,2],pch=6,  col = ifelse(statistics1[,12]>0,6,0),cex= ifelse(statistics1[,12]>0,statistics1[,12]/norm+1,0))
  points(statistics1[,2],pch=1,  col = ifelse(statistics1[,13]>0,1,0),cex= ifelse(statistics1[,13]>0,statistics1[,13]/norm+1,0))
  dev.off()

  jpeg(file=paste(workdir,iii,"plot3.jpg",sep = ""))
  plot(ylim = mlik.lim, xlim = waic.lim,xlab = "WAIC", ylab="MLIK" ,statistics1[,2],statistics1[,1], pch=19, col = 7,cex= 1*((statistics1[,4]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
  points(statistics1[,2],statistics1[,1],pch=8,  col = ifelse(statistics1[,3]>0,5,0),cex= ifelse(statistics1[,3]>0,statistics1[,3]/norm+1,0))
  points(statistics1[,2],statistics1[,1],pch=2,  col = ifelse(statistics1[,9]>0,2,0),cex= ifelse(statistics1[,9]>0,statistics1[,9]/norm+1,0))
  points(statistics1[,2],statistics1[,1],pch=3,  col = ifelse(statistics1[,10]>0,3,0),cex= ifelse(statistics1[,10]>0,statistics1[,10]/norm+1,0))
  points(statistics1[,2],statistics1[,1],pch=4,  col = ifelse(statistics1[,11]>0,4,0),cex= ifelse(statistics1[,11]>0,statistics1[,11]/norm+1,0))
  points(statistics1[,2],statistics1[,1],pch=6,  col = ifelse(statistics1[,12]>0,6,0),cex= ifelse(statistics1[,12]>0,statistics1[,12]/norm+1,0))
  points(statistics1[,2],statistics1[,1],pch=1,  col = ifelse(statistics1[,13]>0,1,0),cex= ifelse(statistics1[,13]>0,statistics1[,13]/norm+1,0))
  dev.off()

  norm1<-(sum(statistics1[,3],na.rm = TRUE))

  xyz<-which(statistics1[,1]<0)
  xyz<-intersect(xyz,which(statistics1[,1]!=-10000))
  moddee<-which(statistics1[xyz,1]==max(statistics1[xyz,1]))[1]
  zyx<-array(data = NA,dim = 2 ^(length(fparam)+1))
  nconsum<-sum(exp(-statistics1[moddee,1]+statistics1[xyz,1]))
  if( nconsum > 0)
  {
    zyx[xyz]<-exp(statistics1[xyz,1]-statistics1[moddee,1])/nconsum
    y.post.lim<-c(0,max(zyx[xyz]))
  }else{

    zyx[xyz]<-statistics1[xyz,3]/norm1
  }


  if(is.nan(y.post.lim[2]))
    y.post.lim[2]<-max(statistics1[,3]/norm1,na.rm = TRUE)
  if(is.nan(y.post.lim[2]))
    y.post.lim[2]<-1

  jpeg(file=paste(workdir,iii,"plot4.jpg",sep = ""))
  plot(xlab = "model_id", ylab="Pr(M(model_id)|D)",ylim = y.post.lim, statistics1[,3]/norm1,pch=19, col = 7,cex= 3*((statistics1[,4]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
  points(statistics1[,3]/norm1,pch=8,  col = ifelse(statistics1[,3]>0,5,0),cex= ifelse(statistics1[,3]>0,statistics1[,3]/norm+1,0))
  points(statistics1[,3]/norm1,pch=2,  col = ifelse(statistics1[,9]>0,2,0),cex= ifelse(statistics1[,9]>0,statistics1[,9]/norm+1,0))
  points(statistics1[,3]/norm1,pch=3,  col = ifelse(statistics1[,10]>0,3,0),cex= ifelse(statistics1[,10]>0,statistics1[,10]/norm+1,0))
  points(statistics1[,3]/norm1,pch=4,  col = ifelse(statistics1[,11]>0,4,0),cex= ifelse(statistics1[,11]>0,statistics1[,11]/norm+1,0))
  points(statistics1[,3]/norm1,pch=6,  col = ifelse(statistics1[,12]>0,6,0),cex= ifelse(statistics1[,12]>0,statistics1[,12]/norm+1,0))
  points(statistics1[,3]/norm1,pch=1,  col = ifelse(statistics1[,13]>0,1,0),cex= ifelse(statistics1[,13]>0,statistics1[,13]/norm+1,0))
  points(zyx,pch=20,  col = 2,cex= 2)
  dev.off()



  jpeg(file=paste(workdir,iii,"plot5.jpg",sep = ""))
  plot(xlim = waic.lim, ylim = y.post.lim,xlab = "WAIC(M)", ylab="Pr(M|D)",statistics1[,2],statistics1[,3]/norm1, pch=19, col = 7,cex= 3*((statistics1[,4]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
  points(statistics1[,2],statistics1[,3]/norm1,pch=8,  col = ifelse(statistics1[,3]>0,5,0),cex= ifelse(statistics1[,3]>0,statistics1[,3]/norm+1,0))
  points(statistics1[,2],statistics1[,3]/norm1,pch=2,  col = ifelse(statistics1[,9]>0,2,0),cex= ifelse(statistics1[,9]>0,statistics1[,9]/norm+1,0))
  points(statistics1[,2],statistics1[,3]/norm1,pch=3,  col = ifelse(statistics1[,10]>0,3,0),cex= ifelse(statistics1[,10]>0,statistics1[,10]/norm+1,0))
  points(statistics1[,2],statistics1[,3]/norm1,pch=4,  col = ifelse(statistics1[,11]>0,4,0),cex= ifelse(statistics1[,11]>0,statistics1[,11]/norm+1,0))
  points(statistics1[,2],statistics1[,3]/norm1,pch=6,  col = ifelse(statistics1[,12]>0,6,0),cex= ifelse(statistics1[,12]>0,statistics1[,12]/norm+1,0))
  points(statistics1[,2],statistics1[,3]/norm1,pch=1,  col = ifelse(statistics1[,13]>0,1,0),cex= ifelse(statistics1[,13]>0,statistics1[,13]/norm+1,0))
  points(statistics1[,2],zyx,pch=20,  col = 2,cex= 2)
  dev.off()

  jpeg(file=paste(workdir,iii,"plot6.jpg",sep = ""))
  plot(xlim = mlik.lim,ylim = y.post.lim, xlab = "MLIK", ylab="Pr(M|D)",statistics1[,1],statistics1[,3]/norm1, pch=19, col = 7,cex= 3*((statistics1[,4]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
  points(statistics1[,1],statistics1[,3]/norm1,pch=8,  col = ifelse(statistics1[,3]>0,5,0),cex= ifelse(statistics1[,3]>0,statistics1[,3]/statistics1[,3]*3 +1,0))
  points(statistics1[,1],statistics1[,3]/norm1,pch=2,  col = ifelse(statistics1[,9]>0,2,0),cex= ifelse(statistics1[,9]>0,statistics1[,9]/statistics1[,3]*3 +1,0))
  points(statistics1[,1],statistics1[,3]/norm1,pch=3,  col = ifelse(statistics1[,10]>0,3,0),cex= ifelse(statistics1[,10]>0,statistics1[,10]/statistics1[,3]*3 +1,0))
  points(statistics1[,1],statistics1[,3]/norm1,pch=4,  col = ifelse(statistics1[,11]>0,4,0),cex= ifelse(statistics1[,11]>0,statistics1[,11]/statistics1[,3]*3 +1,0))
  points(statistics1[,1],statistics1[,3]/norm1,pch=6,  col = ifelse(statistics1[,12]>0,6,0),cex= ifelse(statistics1[,12]>0,statistics1[,12]/statistics1[,3]*3 +1,0))
  points(statistics1[,1],statistics1[,3]/norm1,pch=1,  col = ifelse(statistics1[,13]>0,1,0),cex= ifelse(statistics1[,13]>0,statistics1[,13]/statistics1[,3]*3 +1,0))
  points(statistics1[,1],zyx,pch=20,  col = 2,cex= 2)
  dev.off()


  lldd<-2^(length(fparam)+1)

  moddee<-which(zyx==max(zyx,na.rm = TRUE))
  iidr<-which(statistics1[moddee,1]==max(statistics1[moddee,1],na.rm = TRUE))
  iidr<-which(statistics1[moddee[iidr],3]==max(statistics1[moddee[iidr],3],na.rm = TRUE))
  moddee<-moddee[iidr]
  if(length(moddee)>1)
    moddee<-moddee[1]

  vec<-dectobit(moddee-1)
  varcur<-c(array(0,dim = (length(fparam)+1 -length(vec))),vec)
  df = data.frame(varcur)


  for(i in 1:(lldd -1))
  {
    if(i==moddee)
    {
       next

    }else
    {
      vec<-dectobit(i-1)
      varcur<-c(array(0,dim = (length(fparam)+1 -length(vec))),vec)
      df<-cbind(df,varcur)
      #colnames(x = df)[i] <- paste("solution ",i)
    }
  }
  df<-t(df)

  x<-dist(x = df,method = "binary")

  dists<-c(0,x[1:lldd-1])

  #length(dists)
  #which(dists==0)


  jpeg(file=paste(workdir,iii,"plot7.jpg",sep = ""))
  plot(xlab = "x = |M* - M| = distance from the main mode", ylab="MLIK(M)",ylim = mlik.lim,,y=c(statistics1[moddee,1],statistics1[-moddee,1]),x=dists,pch=19, col = 7,cex= 1*((statistics1[,4]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
  points(y=c(statistics1[moddee,1],statistics1[-moddee,1]),x=dists,pch=8,  col = ifelse(c(statistics1[moddee,3],statistics1[-moddee,3])>0,5,0),cex= ifelse(c(statistics1[moddee,3],statistics1[-moddee,3])>0,c(statistics1[moddee,3],statistics1[-moddee,3])/norm+1,0))
  points(y=c(statistics1[moddee,1],statistics1[-moddee,1]),x=dists,pch=2,  col = ifelse(c(statistics1[moddee,9],statistics1[-moddee,9])>0,2,0),cex= ifelse(c(statistics1[moddee,9],statistics1[-moddee,9])>0,c(statistics1[moddee,9],statistics1[-moddee,9])/norm+1,0))
  points(y=c(statistics1[moddee,1],statistics1[-moddee,1]),x=dists,pch=3,  col = ifelse(c(statistics1[moddee,10],statistics1[-moddee,10])>0,3,0),cex= ifelse(c(statistics1[moddee,10],statistics1[-moddee,10])>0,c(statistics1[moddee,10],statistics1[-moddee,10])/norm+1,0))
  points(y=c(statistics1[moddee,1],statistics1[-moddee,1]),x=dists,pch=4,  col = ifelse(c(statistics1[moddee,11],statistics1[-moddee,11])>0,4,0),cex= ifelse(c(statistics1[moddee,11],statistics1[-moddee,11])>0,c(statistics1[moddee,11],statistics1[-moddee,11])/norm+1,0))
  points(y=c(statistics1[moddee,1],statistics1[-moddee,1]),x=dists,pch=6,  col = ifelse(c(statistics1[moddee,12],statistics1[-moddee,12])>0,6,0),cex= ifelse(c(statistics1[moddee,12],statistics1[-moddee,12])>0,c(statistics1[moddee,12],statistics1[-moddee,12])/norm+1,0))
  points(y=c(statistics1[moddee,1],statistics1[-moddee,1]),x=dists,pch=1,  col = ifelse(c(statistics1[moddee,13],statistics1[-moddee,13])>0,1,0),cex= ifelse(c(statistics1[moddee,13],statistics1[-moddee,13])>0,c(statistics1[moddee,13],statistics1[-moddee,13])/norm+1,0))
  dev.off()

  jpeg(file=paste(workdir,iii,"plot8.jpg",sep = ""))
  plot(xlab = "x = |M* - M| = distance from the main mode", ylab="WAIC(M)",ylim = waic.lim,y=c(statistics1[moddee,2],statistics1[-moddee,2]),x=dists,pch=19, col = 7,cex= 1*((statistics1[,4]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
  points(y=c(statistics1[moddee,2],statistics1[-moddee,2]),x=dists,pch=8,  col = ifelse(c(statistics1[moddee,3],statistics1[-moddee,3])>0,5,0),cex= ifelse(c(statistics1[moddee,3],statistics1[-moddee,3])>0,c(statistics1[moddee,3],statistics1[-moddee,3])/norm+1,0))
  points(y=c(statistics1[moddee,2],statistics1[-moddee,2]),x=dists,pch=2,  col = ifelse(c(statistics1[moddee,9],statistics1[-moddee,9])>0,2,0),cex= ifelse(c(statistics1[moddee,9],statistics1[-moddee,9])>0,c(statistics1[moddee,9],statistics1[-moddee,9])/norm+1,0))
  points(y=c(statistics1[moddee,2],statistics1[-moddee,2]),x=dists,pch=3,  col = ifelse(c(statistics1[moddee,10],statistics1[-moddee,10])>0,3,0),cex= ifelse(c(statistics1[moddee,10],statistics1[-moddee,10])>0,c(statistics1[moddee,10],statistics1[-moddee,10])/norm+1,0))
  points(y=c(statistics1[moddee,2],statistics1[-moddee,2]),x=dists,pch=4,  col = ifelse(c(statistics1[moddee,11],statistics1[-moddee,11])>0,4,0),cex= ifelse(c(statistics1[moddee,11],statistics1[-moddee,11])>0,c(statistics1[moddee,11],statistics1[-moddee,11])/norm+1,0))
  points(y=c(statistics1[moddee,2],statistics1[-moddee,2]),x=dists,pch=6,  col = ifelse(c(statistics1[moddee,12],statistics1[-moddee,12])>0,6,0),cex= ifelse(c(statistics1[moddee,12],statistics1[-moddee,12])>0,c(statistics1[moddee,12],statistics1[-moddee,12])/norm+1,0))
  points(y=c(statistics1[moddee,2],statistics1[-moddee,2]),x=dists,pch=1,  col = ifelse(c(statistics1[moddee,13],statistics1[-moddee,13])>0,1,0),cex= ifelse(c(statistics1[moddee,13],statistics1[-moddee,13])>0,c(statistics1[moddee,13],statistics1[-moddee,13])/norm+1,0))
  dev.off()

  jpeg(file=paste(workdir,iii,"plot9.jpg",sep = ""))
  plot(xlab = "x = |M* - M| = distance from the main mode",ylim = y.post.lim, ylab="Pr(M|D)",y=c(statistics1[moddee,3],statistics1[-moddee,3])/norm1,x=dists,pch=19, col = 7,cex= 3*((statistics1[,4]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
  points(y=c(statistics1[moddee,3],statistics1[-moddee,3])/norm1,x=dists,pch=8,  col = ifelse(c(statistics1[moddee,3],statistics1[-moddee,3])>0,5,0),cex= ifelse(c(statistics1[moddee,3],statistics1[-moddee,3])>0,c(statistics1[moddee,3],statistics1[-moddee,3])/norm+1,0))
  points(y=c(statistics1[moddee,3],statistics1[-moddee,3])/norm1,x=dists,pch=2,  col = ifelse(c(statistics1[moddee,9],statistics1[-moddee,9])>0,2,0),cex= ifelse(c(statistics1[moddee,9],statistics1[-moddee,9])>0,c(statistics1[moddee,9],statistics1[-moddee,9])/norm+1,0))
  points(y=c(statistics1[moddee,3],statistics1[-moddee,3])/norm1,x=dists,pch=3,  col = ifelse(c(statistics1[moddee,10],statistics1[-moddee,10])>0,3,0),cex= ifelse(c(statistics1[moddee,10],statistics1[-moddee,10])>0,c(statistics1[moddee,10],statistics1[-moddee,10])/norm+1,0))
  points(y=c(statistics1[moddee,3],statistics1[-moddee,3])/norm1,x=dists,pch=4,  col = ifelse(c(statistics1[moddee,11],statistics1[-moddee,11])>0,4,0),cex= ifelse(c(statistics1[moddee,11],statistics1[-moddee,11])>0,c(statistics1[moddee,11],statistics1[-moddee,11])/norm+1,0))
  points(y=c(statistics1[moddee,3],statistics1[-moddee,3])/norm1,x=dists,pch=6,  col = ifelse(c(statistics1[moddee,12],statistics1[-moddee,12])>0,6,0),cex= ifelse(c(statistics1[moddee,12],statistics1[-moddee,12])>0,c(statistics1[moddee,12],statistics1[-moddee,12])/norm+1,0))
  points(y=c(statistics1[moddee,3],statistics1[-moddee,3])/norm1,x=dists,pch=1,  col = ifelse(c(statistics1[moddee,13],statistics1[-moddee,13])>0,1,0),cex= ifelse(c(statistics1[moddee,13],statistics1[-moddee,13])>0,c(statistics1[moddee,13],statistics1[-moddee,13])/norm+1,0))
  points(y= c(zyx[moddee],zyx[-moddee]),x=dists,pch=20,  col = 2,cex= 2)
  dev.off()

  # further address subset of the set of the best solution of cardinality 1024
  if(lldd>1024)
  {
    lldd<-1024
    quant<-(sort(statistics1[,1],decreasing = TRUE)[lldd+1])
    indmds<-which(statistics1[,1]>quant)
    length(indmds)

  }else{
    quant<- -Inf
    indmds<-1:lldd
  }





  vec<-dectobit(moddee-1)
  varcur<-c(array(0,dim = (length(fparam)+1 -length(vec))),vec)
  df = data.frame(varcur)



  for(i in 1:(lldd -1))
  {
    if(i==moddee)
    {
      next

    }else
    {
      vec<-dectobit(indmds[i]-1)
      varcur<-c(array(0,dim = (length(fparam)+1 -length(vec))),vec)
      df<-cbind(df,varcur)
      #colnames(x = df)[i] <- paste("solution ",i)
    }
  }
  df<-t(df)

  x<-dist(x = df,method = "binary")

  dists<-c(0,x[1:lldd-1])


  fit.mds <- cmdscale(d = x,eig=FALSE, k=2) # k is the number of dim


  #fit.mds # view results
  x.mds <- fit.mds[,1]
  y.mds <- fit.mds[,2]
  jpeg(file=paste(workdir,iii,"plot10.jpg",sep = ""))
  plot(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])>0,c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])/norm1*100,0))
  points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=2,pch = 19,cex= c(zyx[moddee],zyx[setdiff(indmds, moddee)])*100,0)
  points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])>0,c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])/norm1*100,0))
  points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])>0,c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])/norm1*100,0))
  points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])>0,c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])/norm1*100,0))
  points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])>0,c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])/norm1*100,0))
  points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])>0,c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])/norm1*100,0))
  points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])>0,c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])/norm1*100,0))
  points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])>0,c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])/norm1*100,0))
  points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])>0,c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])/norm1*100,0))
  points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])>0,c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])/norm1*100,0))
  points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])>0,c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])/norm1*100,0))
  points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])>0,c(statistics1[moddee,3],statistics1[setdiff(indmds, moddee),3])/norm1*100,0))
  points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=7,pch = 19,cex= 0.4)
  points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=1,pch = 19,cex= 0.01)
  dev.off()

  }


  fileConn<-file=paste(workdir,"outresults.csv",sep = "")
  #for(i in 2:M)
  #{
  #  write(opt[i],file=fileConn,append=TRUE)
  #}
  writeLines(opt, fileConn,con = ",",sep = "\n")
  close(fileConn)


}),abort = function(){options(error=traceback); onerr<-TRUE})

#
#
# tttr1<-res$time
# res$model
#
#
# tttr<-res$time
# for(i in 1:M)
# {
#   if(is.null(res[[i]]$model))
#   {
#     opt[i+1]<-paste("error occured at cpu ",i)
#     next
#   }
#   result<-res[[i]]
#   #print(summary(result$model))
#   covobs <- fparam[which(result$vars[-1] %in% c(1,3))]
#   obsconst<-as.integer(result$vars[1] %in% c(1,3))
#   globobs<-array(data = 0,dim = length(fparam)+1)
#   result$model$summary.fixed$mean
#   tmp1<-result$model$summary.fixed$mean
#   globobs[which(result$vars %in% c(1,3))] <-tmp1
#   tmp2<-summary(result$model)[[4]][[1]]
#   opt[i+1]<-as.character(paste(c(paste(i),globobs,paste(tmp2,collapse = "*"),result$model$waic[[1]],paste(round(result$p.post.pred,digits = 2 ), collapse = "*"),paste(result$p.post, collapse = "*"),result$eps,result$time), collapse = ";"))
#   #print(i)
# }
#
# statistics1[1281,1]
#
# max(statistics1[,1],na.rm = TRUE)
# # write results to the output file
#
#
# jpeg(file=paste(workdir,iii,"plot11.jpg",collapse = "",sep = ""))
# plot(ylim = mlik.lim, xlab = "model_id", ylab="MLIK" , statistics1[,1],pch=19, col = 7,cex= statistics1[,4]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])
# points(statistics1[,1],pch=8,  col = ifelse(statistics1[,3]>0,5,0),cex= ifelse(statistics1[,3]>0,statistics1[,3]/norm+1,0))
# dev.off()
#
# jpeg(file=paste(workdir,iii,"plot12.jpg",collapse = "",sep = ""))
# plot(ylim = mlik.lim, xlab = "model_id", ylab="MLIK" , statistics1[,1],pch=19, col = 7,cex= statistics1[,4]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])
# points(statistics1[,1],pch=19, col = 5,cex= statistics1[,4])
# points(statistics1[,1],pch=2,  col = ifelse(statistics1[,9]>0,2,0),cex= ifelse(statistics1[,9]>0,statistics1[,9]/norm+1,0))
# dev.off()
#
# jpeg(file=paste(workdir,iii,"plot13.jpg",collapse = "",sep = ""))
# plot(ylim = mlik.lim, xlab = "model_id", ylab="MLIK" , statistics1[,1],pch=19, col = 7,cex= statistics1[,4]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])
# points(statistics1[,1],pch=19, col = 5,cex= statistics1[,5])
# points(statistics1[,1],pch=3,  col = ifelse(statistics1[,10]>0,3,0),cex= ifelse(statistics1[,10]>0,statistics1[,10]/norm+1,0))
# dev.off()
# jpeg(file=paste(workdir,iii,"plot14.jpg",collapse = "",sep = ""))
# plot(ylim = mlik.lim, xlab = "model_id", ylab="MLIK" , statistics1[,1],pch=19, col = 7,cex= statistics1[,4]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])
# points(statistics1[,1],pch=19, col = 5,cex= statistics1[,6])
# points(statistics1[,1],pch=4,  col = ifelse(statistics1[,11]>0,4,0),cex= ifelse(statistics1[,11]>0,statistics1[,11]/norm+1,0))
# dev.off()
# jpeg(file=paste(workdir,iii,"plot15.jpg",collapse = "",sep = ""))
# plot(ylim = mlik.lim, xlab = "model_id", ylab="MLIK" , statistics1[,1],pch=19, col = 7,cex= statistics1[,4]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])
# points(statistics1[,1],pch=19, col = 5,cex= statistics1[,7])
# points(statistics1[,1],pch=6,  col = ifelse(statistics1[,12]>0,6,0),cex= ifelse(statistics1[,12]>0,statistics1[,12]/norm+1,0))
# dev.off()
# jpeg(file=paste(workdir,iii,"plot16.jpg",collapse = "",sep = ""))
# plot(ylim = mlik.lim, xlab = "model_id", ylab="MLIK" , statistics1[,1],pch=19, col = 7,cex= statistics1[,4]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])
# points(statistics1[,1],pch=19, col = 5,cex= statistics1[,8])
# points(statistics1[,1],pch=1,  col = ifelse(statistics1[,13]>0,1,0),cex= ifelse(statistics1[,13]>0,statistics1[,13]/norm+1,0))
# dev.off()
#
#
#

