##################################
source("~/MCLR_library_V3.r")
library(parallel)
options(mc.cores=30) # Specify how many cores you use
########################################################################
###############
### DETERMINE THE MODEL 
###############
###		SETTINGS 
PParam=p=p.param=0.5
n.sym=nsym_fix=nsym=200
n.obs=nobs=nrows=1000
n.lvs=nlvs=7
n.trs=ntrs=ntrs_fix=4
n.itr=nitr=nitr_fix=200000
lambda=lambda_fix=2
A=a=a.start=a.end=0.7
nsnp=n.snp=100
### beta's
Bvec=c(1.5*rep(1,6),3.5,6.6,8.9,9,9.8,7) 
tree.ind=Tind=c(39,38,31,25,23,16,10,4,56,43,20,60,100,5,95,1,64,30,14,47,21,7)#
n.tree.mod=12
tree.lvs=c(4,3,3,2,2,2,1,1,1,1,1,1)
eff.sizes= Bvec[length(Bvec):1]
slope=1
op=c("&","&","&","&","&","&","&","&","&","&","&")
########## model 				
M_2=make.model(X=snpmatrix(nsnp=n.snp,nrows=n.obs,p=PParam),n.tree.mod=n.tree.mod,tree.lvs=tree.lvs,eff.sizes= Bvec[length(Bvec):1],slope=slope,tree.ind=Tind,op=op)
object=M_2
### interactions from the model
true.int=model_int(X=snpmatrix(nsnp=n.snp,nrows=n.obs,p=PParam),n.tree.mod=n.tree.mod,tree.lvs=tree.lvs,eff.sizes= Bvec[length(Bvec):1],slope=slope,tree.ind=Tind,op=op)
true.int=true.int$m.int
if(!is.null(true.int)){
 write.table(true.int, file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  true.int1=read.table(file = "M.int.txt")
  fn <- "M.int.txt"
  if (file.exists(fn)) file.remove(fn)
  }
true.int=apply(t(apply(as.matrix(true.int1),1,tree.merge)),1,intpaste)

################## simulation procedure #####################
fpath=file.path(paste("~/RESULTS/M_S4_L2_k4_",sep =''),fsep="")
dir.create(file.path(paste(fpath,'RES',sep =''),fsep=""))
tmp=proc.time()
res<-mclapply(1:nsym,function(i)
{
dir.create(file.path(paste(fpath,i,sep =''),fsep=""))
setwd(file.path(paste(fpath,i,sep =''),fsep=""))
data.sim.proc(object=M_2,
n.sym=nsym,
n.snp=n.snp,
n.obs=n.obs,
n.lvs=nlvs_fix,
n.trs=ntrs_fix,
n.itr=nitr_fix,
a.start=A,
a.end=A,
lambda=lambda_fix,
p.param=PParam,
n.tree.mod=12,tree.lvs=c(4,3,3,2,2,2,1,1,1,1,1,1),eff.sizes= Bvec[length(Bvec):1],slope=1,tree.ind=Tind,op=c("&","&","&","&","&","&","&","&","&","&","&"),mypath=file.path(paste(fpath,'RES',sep =''),fsep="")
)
}
)
time.S=proc.time()-tmp
time.S

#####################	Summaries ##################	
if(!is.null(true.int)){
 write.table(true.int, file = "~/M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  true.int1=read.table(file = "~/M.int.txt")
  fn <- "M.int.txt"
  if (file.exists(fn)) file.remove(fn)
  }
true.int=apply(t(apply(as.matrix(true.int1),1,tree.merge)),1,intpaste)
##########	base name for a resulting files 	
Out.file="~/RESULTS/M_S4_L2_k4_RES/result_p=05"
######### summaries 	
summarize.ILA(paste(Out.file,"_MCLR",sep=""),A=A,p.param=P,n.sym=n.sym,true.int=true.int,object=object,n.snp=n.snp,n.obs=n.obs,n.itr=n.itr,boxplots.all=FALSE)
summarize.ILA(paste(Out.file,"_poiss",sep=""),A=A,p.param=P,n.sym=n.sym,true.int=true.int,object=object,n.snp=n.snp,n.obs=n.obs,n.itr=n.itr,boxplots.all=FALSE)
summarize.ETA(paste(Out.file,"_MCLR_ETA",sep=""),A=A,p.param=P,n.sym=n.sym,true.int=true.int,object=object,n.snp=n.snp,n.obs=n.obs,n.itr=n.itr,boxplots.all=FALSE)				
summarize.ETA(paste(Out.file,"_poiss_ETA",sep=""),A=A,p.param=P,n.sym=n.sym,true.int=true.int,object=object,n.snp=n.snp,n.obs=n.obs,n.itr=n.itr,boxplots.all=FALSE)
summarize.CETA(paste(Out.file,"_poiss_CETA",sep=""),A=A,p.param=P,n.sym=n.sym,true.int=true.int,object=object,n.snp=n.snp,n.obs=n.obs,n.itr=n.itr,boxplots.all=FALSE)
summarize.CETA(paste(Out.file,"_MCLR_CETA",sep=""),A=A,p.param=P,n.sym=n.sym,true.int=true.int,object=object,n.snp=n.snp,n.obs=n.obs,n.itr=n.itr,boxplots.all=FALSE)







