#####################
#####################
rm(list=ls())
library(gtools)
library(LogicReg)

##########################################################################
vec.trees=function(vec,nlvs,ntrs,indic)
## determine non-empty trees 
## vec... vector of trees from a model in a 1000-2000 notation
# nlvs ... maximum allowed no. of  leaves per tree
# ntrs... maximum allowed numbed of trees in a model 
# if indic=1 returns a nonempty trees from a model vec in a matrix form 
# else it returns a number of nonempty trees in the model vec 
{
stats.seq<- embed(as.numeric(vec), nlvs)# subvectors of length nlvs
# check if there are empty trees and count them
# indices for trees 
t=1 # starting index
It=NULL
while(t<ntrs*nlvs){
It<-c(It,t)
t=t+nlvs
}
S=sum(colSums(t(stats.seq[It,]==rep(0,nlvs)))==nlvs) ## number of empty trees in the ith model
if(indic==1){return(stats.seq[It,ncol(stats.seq[It,]):1])} else{return(nmod.trees=ntrs-S)}
}

##########################################################################
remove.levels=function(f)
{
return(f[which(f %in% names(table(f))[table(f)==1]), drop=TRUE])
}
##########################################################################
is.operator=function(node)
## function to determine if the node is an operator (AND --> 1000 or OR---> 2000)
## node ... an element of a tree in the 1000-2000 notation
## 
{
if((node!=1000)&&(node!=2000)) return(FALSE)
return(TRUE)
}

##########################################################################
is.leaf=function(node)
## function to determine if the node is a a leaf 
## node ... an element of a tree in the 1000-2000 notation
{
if((node!=1000)&&(node!=2000)&&(node!=0)){return(TRUE)} else {return(FALSE)}
}
##########################################################################
vec.leaves=function(vec)
## function to determine how many leaves are there in a tree 
## vec... a tree in the 1000-2000 notation
{
return(sum(sapply(as.numeric(vec),is.leaf)))
}

##########################################################################
Lfun=function(lvs,nsnp)
## function to calculate the number of all possible trees with lvs leaves which can be obtain while considering nsnp snps
## THIS APPROACH CONTAINS MANY REPETITIONS, I.E. SOME MODELS ARE COUNTED MORE THAN ONCE (SO IN PRINCIPLE THE TOTAL NUMBER IS TOO LARGE)
{
L=NULL
l=NULL
if(lvs==1) L=2*nsnp
if(lvs>1){
pp=1
for(j in 2:lvs){
l[pp]=2^(lvs-2)*choose(nsnp,j)*4^j
pp=pp+1
}
L=sum(l)
}
return(L)
}

##########################################################################
n.mod<-function(nsnp,tr,nlvs)
## function to  count possible number of models with tr trees and nlvs leaves
# nsnp number of snps
# tr number of trees in the model 
# nlvs number of leaves in the model 
{
maxlvs=(nlvs+1)/2# max number of leaves which can be put on nlvs positions 
s=c(tr:(maxlvs*tr))
al.vec=NULL
res=t(as.matrix(rep(0,4)))
## all permutations
if(tr>0){
al.vec= combinations(n=maxlvs, r=tr, v=1:nlvs, set=TRUE, repeats.allowed=TRUE)
N=vector("list",nrow(al.vec))
Lmat=NULL
for( i in 1:nrow(al.vec)){
N[[i]][1]=Lfun(al.vec[i,1],nsnp)
if(tr>1){
for(k in 2:tr){
if(al.vec[i,k]==al.vec[i,k-1]) N[[i]][k]=N[[i]][k-1]-1
if(al.vec[i,k]!=al.vec[i,k-1]) N[[i]][k]=Lfun(al.vec[i,k],nsnp)
}
}
Lmat=rbind(Lmat,N[[i]])
}
# product of rows : 
res.mat=NULL
res.mat=apply(Lmat, 1, prod)
## sums for equal numbers of leaves 
Is<-vector('list',length(s))
for(i in 1:nrow(al.vec)){ 
is=NULL
for(j in 1:length(s))
if(rowSums(al.vec)[i]==s[j]) Is[[j]]<-c(Is[[j]],i)
}
R=NULL
for(i in 1:length(s)){
R[i]=sum(res.mat[Is[[i]]])
}
res=NULL
res=cbind(rep(nsnp, length(s)),rep(tr, length(s)), s, R)
}
colnames(res)=NULL
#colnames(res)=c('nsnp','no.of trees','no. of leaves', 'no.of models')
return(res)
}
##########################################################################
intpaste=function(vec)
{
paste(vec,collapse=" ")
}
##########################################################################
tree.merge=function(vec)
## function to transform a tree (vec) to a unique representation 
# vec is a vector containing a tree in a 1000-2000 notation
{
new.tree=NULL
vec=as.numeric(vec)
sorted.tree=NULL# a vector for a sorted tree
nlvs=length(vec)
if(is.leaf(vec[1])){ # no interaction in a tree
sorted.tree=vec
}else{if(is.leaf(vec[2])){ # leaf on the second position 
sorted.tree[1]=vec[1]
i=2
pair=c(i,i+1)
pair.end=c(nlvs-1,nlvs)
vv=sign(vec[pair])*sort(abs(vec[pair]),decreasing=TRUE)
sorted.tree[pair]=vv
if(nlvs==7){sorted.tree[pair.end]=c(0,0)
sorted.tree[pair.end-2]=vec[pair.end]}
}else{sorted.tree=vec}}

if(is.list(sorted.tree)){svec=as.vector(unlist(sorted.tree))
						}else{svec=sorted.tree}

n.op=NULL#number of DIFFERENT operators	
op.ind=NULL
op.ind=as.vector(which(unlist(lapply(vec,is.operator))))			
n.op=length(unique(vec[op.ind]))
if(n.op==2){#if different op-s: sort within a prime implicant
for(i in 1:((nlvs-1)/2)){
if(!is.leaf(sorted.tree[i])&&(is.leaf(sorted.tree[2*i]))&&(is.leaf(sorted.tree[2*i+1]))){
								pair=c(2*i,2*i+1)
								svec[pair]=sort(abs(sorted.tree[pair]),decreasing=TRUE)
								}
				sorted.tree=svec
				for(ii in 1:length(which(vec<0)))#add minuses where needed
				sorted.tree[which(abs(vec[which(vec<0)])[ii]==sorted.tree)]=-sorted.tree[which(abs(vec[which(vec<0)])[ii]==sorted.tree)]
				}
			}else{sorted.tree=sort(abs(svec),decreasing=TRUE)
			      for(ii in 1:length(which(vec<0)))#add minuses where needed
				sorted.tree[which(abs(vec[which(vec<0)])[ii]==sorted.tree)]=-sorted.tree[which(abs(vec[which(vec<0)])[ii]==sorted.tree)]}
## reducing equivalent representations
if(nlvs==7){
op.no=NULL
op.no=length(op.ind)# number of ALL operators
if(op.no==3){
			if(n.op==1){if(vec[op.ind][1]==2000) sorted.tree=c(rep(1000,3),-sorted.tree[4:7])
						n.un.var=NULL
						n.un.var=length(unique(sorted.tree[4:7]))
						if(n.un.var<4)
						sorted.tree=c(rep(sorted.tree[op.ind][1],(n.un.var-1)),unique(sorted.tree[4:7]),rep(0,(8-2*n.un.var)))
						}else{
							new.tree=sorted.tree
							if(length(which(vec[op.ind]==c(1000,2000,2000)))==3)
							new.tree=c(2000,rep(1000,2),-sorted.tree[4:7])
							if(length(which(vec[op.ind]==c(1000,2000,1000)))==3)
							new.tree=c(2000,1000,2000,-sorted.tree[4:7])
							if(length(which(vec[op.ind]==c(1000,1000,2000)))==3)
							new.tree=c(2000,2000,1000,-sorted.tree[4:7])
							sorted.tree=new.tree
							##having sorted tree-reduce trees with repeated var-s
							#if(length(which(sorted.tree[4:5]%in%sorted.tree[6:7]))==1)
							if((length(unique(sorted.tree[4:7]))==3)&&(length(unique(sorted.tree[4:5]))==2)&&(length(unique(sorted.tree[6:7]))==2))
								{
								if(length(which(sorted.tree[op.ind]==c(2000,1000,1000)))==3)
								{
								rep.ind=which(sorted.tree[4:5]%in%sorted.tree[6:7])
								rep.ind2=which(sorted.tree[6:7]==sorted.tree[4:5][rep.ind])
								new.tree=c(2000,1000,-sorted.tree[4:5][rep.ind],-sort(c(abs(sorted.tree[4:5][-rep.ind]),abs(sorted.tree[6:7][-rep.ind2])),decreasing=TRUE),rep(0,2))
								}
								if(length(which(sorted.tree[op.ind]==c(2000,1000,2000)))==3)
								new.tree=c(2000,sorted.tree[6:7],rep(0,4))
								if(length(which(sorted.tree[op.ind]==c(2000,2000,1000)))==3)
								new.tree=c(2000,sorted.tree[4:5],rep(0,4))
								}
							if((length(unique(sorted.tree[4:7]))==2)&&(length(unique(sorted.tree[4:5]))==2)&&(length(unique(sorted.tree[6:7]))==2))
								{
								if(length(which(sorted.tree[op.ind]==c(2000,1000,1000)))==3)
								{
								new.tree=c(1000,sorted.tree[4:5],rep(0,4))
								}else{
								new.tree=c(2000,sorted.tree[4:5],rep(0,4))}
								}	
								
							if(length(unique(sorted.tree[4:5]))==1){
								new.tree=c(sorted.tree[1],sorted.tree[3:4],sorted.tree[6:7],rep(0,2))
								if(length(unique(sorted.tree[6:7]))==1)
								new.tree=c(sorted.tree[1],sign(c(sorted.tree[4],sorted.tree[6]))*sort(abs(c(sorted.tree[4],sorted.tree[6])),decreasing=TRUE),rep(0,4))
								sorted.tree<-new.tree
								}else{if(length(unique(sorted.tree[6:7]))==1)
									 new.tree<-c(sorted.tree[1:2],sorted.tree[6],sorted.tree[4:5],rep(0,2))
									 }
							sorted.tree<-new.tree		 
							}
			vec1=NULL				
			vec1=sorted.tree
			n.op=NULL#number of DIFFERENT operators	
			op.ind=NULL
			op.ind=as.vector(which(unlist(lapply(vec1,is.operator))))			
			n.op=length(unique(vec1[op.ind]))
			op.no=NULL
			op.no=length(op.ind)# number of ALL operators
			}

if(op.no==2){
			if(n.op==1){if(sorted.tree[op.ind][1]==2000) sorted.tree=c(rep(1000,2),-sorted.tree[3:5],rep(0,2))
						sorted.tree=c(sorted.tree[1:2],sign(sorted.tree[3:5])*sort(abs(sorted.tree[3:5]),decreasing=TRUE),rep(0,2))
						if(length(unique(sorted.tree[3:5]))==2) sorted.tree=c(1000,sign(unique(sorted.tree[3:5]))*sort(abs(unique(sorted.tree[3:5])),decreasing=TRUE),rep(0,4))
						}else{
						if(length(which(vec[op.ind]==c(1000,2000)))==2)
							sorted.tree=c(2000,1000,-sorted.tree[3:5],rep(0,2))
						if(length(unique(sorted.tree[3:5]))==2){
							where.rep=which(sorted.tree==sorted.tree[op.no+which(duplicated(sorted.tree[3:5]))])
							if((length(which(where.rep==c(3,4)))==2)||(length(which(where.rep==c(3,5)))==2)) sorted.tree=c(sorted.tree[3],rep(0,6))
							if(length(which(where.rep==c(4,5)))==2) sorted.tree=c(2000,sorted.tree[3:4],rep(0,4))
							}
							}
						
			}

##reduction of two-way 2000-complimented expressions
if(sorted.tree[1]==2000){
if(length(which(sign(sorted.tree)[2:3]==c(-1,-1)))==2)
sorted.tree=c(1000,-sorted.tree[2:3],rep(0,4))
}	
}		
return(sorted.tree)
}

##########################################################################		
tree.operators=function(vec)
## returns all operators from a tree (vec) 
{
p=sapply(vec,is.operator)
return(vec[which(p)])
}		

##########################################################################
Int.tree.ind=function(mat,ntrs,nlvs,lev=1)
{
if(lev==1){
return(which(rowSums(matrix(sapply(mat,is.operator),nrow=ntrs,ncol=nlvs))>0))}
else{return(intersect(which(rowSums(matrix(sapply(mat,is.leaf),nrow=ntrs,ncol=nlvs))==1),which(rowSums(matrix(sapply(mat,is.operator),nrow=ntrs,ncol=nlvs))==0)))}#trees in ith model containing main effects
 }

#########################################################################
 take.adds=function(mat,ind)
 {
 return(mat[ind,1])
 }
##########################################################################
 take.ints=function(mat,ind)
 {
 return(mat[ind,])
 } 
##########################################################################
count.rows <- function(x) 
#This function returns a data.frame with columns counts with the row counts, 
# ind containing a list of row indices into the original x and the unique rows of x.
{ 
  if (is.matrix (x) && (dim (x) [2] == 1))
    x <- as.vector (x) 
 
  order.x <- do.call(order,as.data.frame(x))
   
  if (is.vector (x)) {
    equal.to.previous <-
      x[tail(order.x,-1)] == x[head(order.x,-1)]
  } else {
    equal.to.previous <-
      rowSums(x[tail(order.x,-1),] != x[head(order.x,-1),])==0
  }	
 
  indices <-  split (order.x, cumsum (c (TRUE, !equal.to.previous)))
 
  if (is.vector (x)) {
    x <- x [sapply (indices, function (x) x [[1]]), drop = FALSE]
  } else {
    x <- x [sapply (indices, function (x) x [[1]]), ,drop = FALSE]
  }
  
  data.frame (counts = sapply (indices, length) , ind = I(indices), x) 
}

##########################################################################
model.ind=function(nmtrs)
# FUNCTION TO CREATE A VECTOR CONTAINING AN INDEX OF A MODEL FOR A CORRESPONDING TREE 
# nmtrs is a vector containing numbers of trees for all visited models
{
i=1
model.ind=NULL
while(i<=length(nmtrs)){
model.ind=c(model.ind,rep(i,nmtrs[i]))
i=i+1
}
return(model.ind)
}

##########################################################################
model.prior.poiss=function(lambda=lambda,ntrs=4,nlvs=7,nmodtrs,nmodlvs,nsnp)
### MODEL PRIOR --  WITH POISSON DISTRIBUTION
# lambda ... a parameter of a poissson distribution
# ntrs ... maximum no. of trees allowed
# nlvs... maximum no. of leaves allowed
# nmodtrs ... no. of trees in the model
# nmodlvs ... no. of leaves in the model
{
maxlvs=(nlvs+1)/2# max number of leaves which can be put on nlvs positions 
Ck=NULL
Ck=((lambda^nmodtrs)*exp(-lambda))/factorial(nmodtrs) #(1-aParam)/(aParam^nmodtrs-aParam^(maxlvs*nmodtrs+1))
t0=which(nmodtrs==0)
if(length(t0)>0){
		nmod=sapply(unique(nmodtrs[-which(nmodtrs==0)]),n.mod.poiss,nsnp=nsnp,nlvs=nlvs)
		if(!is.list(nmod)){
		pom=matrix(nmod,ncol=4)
		nmod=NULL
		nmod=pom
		}
		Ind=NULL
		N=NULL
		if(length(unique(nmodtrs[-which(nmodtrs==0)]))==1){
		Ind=which(nmodtrs==unique(nmodtrs[-which(nmodtrs==0)]))#save indices of models which have unique(nmodtrs)[i] trees
		rows=nmodlvs[Ind]-nmodtrs[Ind]+1# which rows of nmod should be taken
		N[Ind]=nmod[rows,4]
		}else{
			for( i in 1:length(unique(nmodtrs[-which(nmodtrs==0)]))){
			Ind=which(nmodtrs==unique(nmodtrs[-which(nmodtrs==0)])[i])#save indices of models which have unique(nmodtrs)[i] trees
			rows=nmodlvs[Ind]-nmodtrs[Ind]+1# which rows of nmod should be taken
			N[Ind]=nmod[[i]][rows,4]
			}
			}
			N[t0]=1
}else{
	nmod=sapply(unique(nmodtrs),n.mod.poiss,nsnp=nsnp,nlvs=nlvs)
	if(!is.list(nmod)){
	pom=matrix(nmod,ncol=4)
	nmod=NULL
	nmod=pom
	}
	Ind=NULL
	N=NULL
	if(length(unique(nmodtrs))==1){
	Ind=which(nmodtrs==unique(nmodtrs))#save indices of models which have unique(nmodtrs)[i] trees
	rows=nmodlvs[Ind]-nmodtrs[Ind]+1# which rows of nmod should be taken
	N[Ind]=nmod[rows,4]
	}else{
	for( i in 1:length(unique(nmodtrs))){
	Ind=which(nmodtrs==unique(nmodtrs)[i])#save indices of models which have unique(nmodtrs)[i] trees
	rows=nmodlvs[Ind]-nmodtrs[Ind]+1# which rows of nmod should be taken
	N[Ind]=nmod[[i]][rows,4]
	}
	}
	}
return(Ck/N)
}

##########################################################################
Lfun.poiss=function(lvs,nsnp)
## function to calculate the number of all possible trees with lvs leaves which can be obtained while 
## considering nsnp snps
## THIS APPROACH CONTAINS MANY REPETITIONS, I.E. SOME MODELS ARE COUNTED MORE THAN ONCE (SO IN PRINCIPLE THE TOTAL NUMBER IS TOO LARGE)
{
return((nsnp^lvs)*2^(3*lvs-2)/factorial(lvs))
}

##########################################################################
n.mod.poiss<-function(nsnp,tr,nlvs)
## function to  count possible number of models with nmod.trs trees and nmod.lvs leaves
# tr number of trees in the model 
# nlvs number of leaves in the model 
{
maxlvs=(nlvs+1)/2# max number of leaves which can be put on nlvs positions 
s=c(tr:(maxlvs*tr))
al.vec=NULL
res=t(as.matrix(rep(0,4)))
## all permutations
if(tr>0){
al.vec= combinations(n=maxlvs, r=tr, v=1:nlvs, set=TRUE, repeats.allowed=TRUE)
N=vector("list",nrow(al.vec))
Lmat=NULL
for( i in 1:nrow(al.vec)){
N[[i]][1]=Lfun.poiss(al.vec[i,1],nsnp)
if(tr>1){
for(k in 2:tr){
if(al.vec[i,k]==al.vec[i,k-1]) N[[i]][k]=N[[i]][k-1]-1
if(al.vec[i,k]!=al.vec[i,k-1]) N[[i]][k]=Lfun.poiss(al.vec[i,k],nsnp)
}
}
Lmat=rbind(Lmat,N[[i]])
}
# product of rows : 
res.mat=NULL
res.mat=apply(Lmat, 1, prod)
## sums for equal numbers of leaves 
Is<-vector('list',length(s))
for(i in 1:nrow(al.vec)){ 
is=NULL
for(j in 1:length(s))
if(rowSums(al.vec)[i]==s[j]) Is[[j]]<-c(Is[[j]],i)
}
R=NULL
for(i in 1:length(s)){
R[i]=sum(res.mat[Is[[i]]])
}
res=NULL
res=cbind(rep(nsnp, length(s)),rep(tr, length(s)), s, R)
}
colnames(res)=NULL
#colnames(res)=c('nsnp','no.of trees','no. of leaves', 'no.of models')
return(res)
}

##########################################################################
contain.add<-function(add,vec)
#function to check if a vector contains a particular variable
{
return(!is.na(match(add,vec)))
}
##########################################################################
variables.from.model<-function(vec)
#function to return indices of variables in the model
{
vec0=vec[which(abs(vec)!=0)]
vecOp=vec0[which(vec0<1000)]
return(unique(abs(vecOp)))
}

##########################################################################
models.cont.expr=function(vec, List)
# function to check which elements of List contain all elements of vec
# List... a list of models
# vec ... a tree in a 1000-2000 notation 
{
res=NULL
for(i in 1:length(List)){
if(sum(vec%in%List[[i]])>=2) res<-c(res,i) #==length(unI[[1]])
}
return(res)
}
##########################################################################

##########################################################################
models.cont.add=function(add, List)
# function to check which elements of List contain an add element 
# List... a list of models
# add ... index of an add. variable 
{
res=NULL
for(i in 1:length(List)){
if(sum(add%in%List[[i]])>=1) res<-c(res,i) #==length(unI[[1]])
}
return(res)
}
##########################################################################
subvectors=function(vec)
#function to find all sub-pairs, sub-triples and fours of vector vec indices
{
ind.pair=combinations(length(vec),2)
subpairs=matrix(vec[ind.pair],ncol=2)
subtriples=NULL
fours=NULL
if(length(vec)>3)
{ind.triple=combinations(length(vec),3)
subtriples=matrix(vec[ind.triple],ncol=3)
fours=vec
}
if(length(vec)==3) subtriples=t(as.matrix(vec))
return(list(subpairs=subpairs,subtriples=subtriples,fours=fours))
}
##########################################################################
subv.in.tree=function(vec, List)
## which models from the List contain a tree vec 
# returns a vector of indices of these models
{
res=NULL
for(i in 1:length(List)){
#if is  matrix
if(is.matrix(List[[i]])){
if(nrow(List[[i]])>0){
for(j in 1:nrow(List[[i]]))
if(sum(vec%in%abs(List[[i]][j,]))==length(vec)) res<-c(res,i) 
}
}
if(!is.matrix(List[[i]])){
V=variables.from.model(as.numeric(List[[i]]))
if(sum(vec%in%V)==length(vec)) res<-c(res,i) 
}
}
return(res)
}
##############################################
######################
model.PI.poiss=function(filename="slogiclisting.tmp",nlvs=3,ntrs=4, lambda=4,nrows=1000,outFile="mypost_results",nsnp=nsnp)
# function to get all prime implicants and their posterior
# for each model and for each separate tree in a model
# creates files with lists of expressions and corresponding posteriors
{
outFile1=paste(outFile,"_ILA.txt",sep="")
outFile3=paste(outFile,"_CETA",sep="")
stats=read.table(filename)
scores=stats[,2]
stats=stats[,4:ncol(stats)]
stats.seq=data.frame(t(stats))
stats.seq=as.list(stats.seq)
T=lapply(stats.seq,vec.trees,nlvs=nlvs,ntrs=ntrs,indic=1)
nmlvs=apply(stats,1,vec.leaves)
nmtrs=apply(stats,1,vec.trees,nlvs=nlvs,ntrs=ntrs,indic=0)
nm=length(T)# number of models on the list 
#########################################
# count priors for models 
priors=NULL
priors=model.prior.poiss(lambda,ntrs,nlvs,nmodtrs=nmtrs,nmodlvs=nmlvs,nsnp=nsnp)
#########################################
# newscores
newscores=NULL
newscores=-scores#-0.5*nmtrs*log(nrows)
Enewscores=NULL
Enewscores=exp(newscores)
############################################
#Addind=lapply(T,Int.tree.ind,ntrs=ntrs,nlvs=nlvs,lev=0)# in vhich tree additive effects appear
PIa=apply(stats,1,variables.from.model)#which variables contains the specific model# 
alladds=as.numeric(unlist(PIa))
#############
# list of all different additive effects 
unadds=unique(sort(abs(alladds)))
#########################################
# take an additive effect , take the prior of all models that contain it, take corresponding score and sum the prior times the new score
Ia=vector("list",length(unadds))
for( i in 1:length(unadds)){
Ij=NULL# models which contain a particular add effect
for( j in 1:length(PIa)) 
if((!is.na(match(unadds[i],PIa[[j]])))||(!is.na(match(-unadds[i],PIa[[j]]))))# ith additive effect appears in jth model
Ij<-c(Ij,j)
Ia[[i]]=Ij# save indices of such models 
}
#########################################
## 
addposts=NULL
res.add=NULL
for(i in 1:length(Ia))
addposts[i]=sum(priors[Ia[[i]]]*Enewscores[Ia[[i]]])/sum(priors*Enewscores)
foo=cbind(unadds,round(addposts,digits=6))
if(length(which(foo[,2]>0))>0){
res.add=foo[order(foo[,2],foo[,1],decreasing=TRUE),]
res.add=res.add[which(res.add[,2]>0),]
}
###### interactions 
Iind=lapply(T,Int.tree.ind,ntrs=ntrs,nlvs=nlvs,lev=1)
PIl=vector("list",length(Iind))
for(ww in 1:length(Iind)) PIl[[ww]]=take.ints(mat=T[[ww]],ind=Iind[[ww]]) 
PIi=do.call(rbind, PIl)#unlisted- as a matrix with trees in rows
res.add.p=NULL
res.int2=NULL
res.int3=NULL
res.int4=NULL

if(length(PIi)>0){
PIm=t(apply(PIi,1,tree.merge))# posortuj zawartosc drzew
PI=NULL
if(nrow(PIm)>2){PI=count.rows(PIm)
				Iint=PI[,2]
				m.ind=model.ind(nmtrs)
				Iint=lapply(Iint,function(i) m.ind[i])
				unints=PI[,c(3:ncol(PI))]
				unints=t(apply(unints,1,tree.merge))
				unI=apply(unints,1,variables.from.model)
				}
if(nrow(PIm)==1){
				PI=t(matrix(c(1,1,PIm),nrow=9))
				Iint=PI[,2]
				m.ind=model.ind(nmtrs)
				Iint=lapply(Iint,function(i) m.ind[i])
				unints=PI[,c(3:ncol(PI))]
				unints=tree.merge(unints)
				unI=variables.from.model(unints)
				unI=t(as.matrix(unI))
				}
if(nrow(PIm)==2){PI=unique(PIm)
				Iint=vector("list",2)
				if(nrow(PI)==1) Iint[[1]]=c(1,2)
				if(nrow(PI)==2){Iint[[1]]=1
								Iint[[2]]=2
								}
				m.ind=model.ind(nmtrs)
				Iint=lapply(Iint,function(i) m.ind[i])
				unints=PI
				unints=t(apply(unints,1,tree.merge))
				unI=apply(unints,1,variables.from.model)
				}
# #########################################
# # take an interaction , take the prior of all models that contain it, take corresponding score and sum the prior times the new score
#########################################
if(nrow(PIm)>1){
if(!is.list(unI)){
x=unI
unI=split(x, rep(1:ncol(x), each = nrow(x)))
}
}else{pom=vector("list",1)
		pom[[1]]<-unI
		unI<-pom}

L2=lapply(unI,models.cont.expr,List=PIa)# lista indeksow modeli zawierajacych co najmniej dwie zmienne sposrod unI ---> dalej spr czy sa w tym samym drzewie (tylko dla tych modeli)
IntLength=as.numeric(unlist(lapply(unI,length)))#dlugosc interakcji
unI2=NULL
unI3=NULL
sub.unI3=NULL
un.2=NULL
if(length(which(IntLength==2))>0) unI2=unI[which(IntLength==2)]# dla tych nie liczyc podwektorow
# do unI2 nie liczyc podwektorow, wszystkie maja dlugosc 2 
Pairs3=NULL
if(length(which(IntLength>2))>0){
unI3=unI[which(IntLength>2)]#dla tych znalezc podwektory
sub.unI3=lapply(unI3,subvectors)
# all pairs being subvectors of larger vectors
for(i in 1:length(sub.unI3)) Pairs3<-rbind(Pairs3,sub.unI3[[i]]$subpairs)
Pairs3=t(apply(Pairs3,1,sort))# sort pairs within rows 
}
unPairs3=unique(Pairs3)# take only unique pairs 
mat.unI2=NULL
s.unI2=NULL
un.unI2=NULL 
un.2=NULL
if(!is.null(unI2)){
mat.unI2=t(matrix(unlist(unI2),nrow=2))# pairs from unI2 as matrix
s.unI2=t(apply(mat.unI2,1,sort))
un.unI2=unique(s.unI2)# take only unique pairs 
}
un.2=unique(rbind(unPairs3,un.unI2))# connect subpairs of larger vec and vec length 2 , then take only unique ones--> these are all pairs to be checked 
# which models contain them in one tree
# znalezc modele ktore maja te pary w jednym drzewie 
# subv.in.tree(un.2[1,],PIl)# indices of models containing a specific subvector
Lsub2=apply(un.2,1,subv.in.tree,List=PIl)# a list-indices of models containing a specific subpair (pair)--> for each pair- a vector of indices 
int2posts=NULL
foo=NULL
sig.int2=NULL
res.int2=NULL
if(length(Lsub2)>0){
for(i in 1:length(Lsub2))
int2posts[i]=sum(priors[Lsub2[[i]]]*Enewscores[Lsub2[[i]]])/sum(priors*Enewscores)
foo=cbind(un.2,round(int2posts,digits=6))
res.int2=foo[order(foo[,(ncol(foo))],decreasing=TRUE),]
sig.int2=NULL
if(is.matrix(res.int2)){
if(length(which(res.int2[,ncol(foo)]>0))>0){
is=which(res.int2[,ncol(foo)]>0)
res.int2=res.int2[is,]
}else{res.int2=res.int2[-c(1:nrow(res.int2)),]}
if(is.matrix(res.int2)){
if(nrow(res.int2)>0){
sigint2<-res.int2[,-ncol(res.int2)]
C1=apply(sigint2,1,intpaste)
sig.int2=cbind(C1,res.int2[,ncol(res.int2)])
}
}else{sigint2<-res.int2[-length(res.int2)]
C1=intpaste(sigint2)
sig.int2=c(C1,res.int2[length(res.int2)])}
}

if(!is.matrix(res.int2)){
if(length(res.int2)>0){
sigint2<-res.int2[-length(res.int2)]
C1=intpaste(sigint2)
sig.int2=c(C1,res.int2[length(res.int2)])
}
}
}

sig.int3=NULL
res.int3=NULL
Triples=NULL# 
unTriples=NULL
if(length(sub.unI3)>0){
for(i in 1:length(sub.unI3)) Triples<-rbind(Triples,sub.unI3[[i]]$subtriples)
Triples=t(apply(Triples,1,sort))# sort triples within rows 
unTriples=unique(Triples)# take only unique triples
Lsub3=apply(unTriples,1,subv.in.tree,List=PIl)# a list-indices of models containing a specific subtriple (triple)--> for each triple- a vector of indices 
## data frame saved to file 
int3posts=NULL
foo=NULL
sig.int3=NULL
res.int3=NULL
if(length(Lsub3)>0){
for(i in 1:length(Lsub3))
int3posts[i]=sum(priors[Lsub3[[i]]]*Enewscores[Lsub3[[i]]])/sum(priors*Enewscores)
foo=cbind(unTriples,round(int3posts,digits=6))
res.int3=foo[order(foo[,(ncol(foo))],decreasing=TRUE),]
sig.int3=NULL
if(is.matrix(res.int3)){
if(length(which(res.int3[,ncol(foo)]>0))>0){
is=which(res.int3[,ncol(foo)]>0)
res.int3=res.int3[is,]
} else{res.int3=res.int3[-c(1:nrow(res.int3)),]}
if(is.matrix(res.int3)){
if(nrow(res.int3)>0){
sigint3<-res.int3[,-ncol(res.int3)]
C1=apply(sigint3,1,intpaste)
sig.int3=cbind(C1,res.int3[,ncol(res.int3)])
}}else{sigint3<-res.int3[-length(res.int3)]
C1=intpaste(sigint3)
sig.int3=c(C1,res.int3[length(res.int3)])}
}
if(!is.matrix(res.int3)){
if(length(res.int3)>0){
sigint3<-res.int3[-length(res.int3)]
C1=intpaste(sigint3)
sig.int3=c(C1,res.int3[length(res.int3)])
}
}
}
}
sig.int4=NULL
res.int4=NULL
sigint4=NULL
Fours=NULL# all fours
unFours=NULL
if(length(which(IntLength>3))>0){
if(length(sub.unI3)>0){
for(i in 1:length(sub.unI3)) Fours<-rbind(Fours,sub.unI3[[i]]$fours)
Fours=t(apply(Fours,1,sort))# sort fours within rows 
unFours=unique(Fours)# take only unique fours
Lsub4=apply(unFours,1,subv.in.tree,List=PIl)# a list-indices of models containing a specific subpair (pair)--> for each pair- a vector of indices 
## data frame saved to file 
int4posts=NULL
foo=NULL
if(length(Lsub4)>0){
for(i in 1:length(Lsub4))
int4posts[i]=sum(priors[Lsub4[[i]]]*Enewscores[Lsub4[[i]]])/sum(priors*Enewscores)
foo=cbind(unFours,round(int4posts,digits=6))
res.int4=foo[order(foo[,(ncol(foo))],decreasing=TRUE),]
sig.int4=NULL
if(is.matrix(res.int4)){
if(length(which(res.int4[,ncol(foo)]>0))>0){
is=which(res.int4[,ncol(foo)]>0)
res.int4=res.int4[is,]
} else{res.int4=res.int4[-c(1:nrow(res.int4)),]}
if(is.matrix(res.int4)){
if(nrow(res.int4)>0){
sigint4<-res.int4[,-ncol(res.int4)]
C1=apply(sigint4,1,intpaste)
sig.int4=cbind(C1,res.int4[,ncol(res.int4)])
}}else{sigint4<-res.int4[-length(res.int4)]
C1=intpaste(sigint4)
sig.int4=c(C1,res.int4[length(res.int4)])}
}
if(!is.matrix(res.int4)){
if(length(res.int4)>0){
sigint4<-res.int4[-length(res.int4)]
C1=intpaste(sigint4)
sig.int4=c(C1,res.int4[length(res.int4)])
}
}
}
}
}


sig.int<-rbind(sig.int2,sig.int3,sig.int4)

int2=NULL
if(!is.null(sig.int2)){
if(is.matrix(sig.int2))
int2<-cbind(as.numeric(sigint2[,1]),as.numeric(sigint2[,2])
,as.numeric(sig.int2[,2]))
if(!is.matrix(sig.int2))
int2<-c(as.numeric(sigint2[1]),as.numeric(sigint2[2]),as.numeric(sig.int2[2]))
}
int3=NULL
if(!is.null(sig.int3)){
if(is.matrix(sig.int3))
int3<-cbind(as.numeric(sigint3[,1]),as.numeric(sigint3[,2]),as.numeric(sigint3[,3]),as.numeric(sig.int3[,2]))
if(!is.matrix(sig.int3))
int3<-c(as.numeric(sigint3[1]),as.numeric(sigint3[2]),as.numeric(sigint3[3]),as.numeric(sig.int3[2]))
}
int4=NULL
if(!is.null(sig.int4)){
if(is.matrix(sig.int4))
int4<-cbind(as.numeric(sigint4[,1]),as.numeric(sigint4[,2]),as.numeric(sigint4[,3]),as.numeric(sigint4[,4]),as.numeric(sig.int4[,2]))
if(!is.matrix(sig.int4))
int4<-c(as.numeric(sigint4[1]),as.numeric(sigint4[2]),as.numeric(sigint4[3]),as.numeric(sigint4[4]),as.numeric(sig.int4[2]))
}

res.int2=NULL
if(!is.null(int2)){
if(is.matrix(int2)){
				if(nrow(int2)>1){
				mpom=rbind(matrix(rep(1000,nrow(int2)),ncol=nrow(int2)),t(as.matrix(int2[,1:2])),matrix(rep(0,4*nrow(int2)),ncol=nrow(int2)))
				res.int2=cbind(apply(t(mpom),1,intpaste),int2[,3])
				}else{
				mpom=rbind(matrix(rep(1000,nrow(int2)),ncol=nrow(int2)),as.matrix(int2[,1:2]),matrix(rep(0,4*nrow(int2)),ncol=nrow(int2)))
				res.int2=cbind(apply(t(mpom),1,intpaste),int2[,3])
				}
				}
if(!is.matrix(int2)){
				mpom=c(1000,int2[1:2],rep(0,4))		
				res.int2=c(intpaste(mpom),int2[3])
				res.int2=t(as.matrix(res.int2))
				}						
				}
				
				
res.int3=NULL
if(!is.null(int3)){
if(is.matrix(int3)){
				if(nrow(int3)>1){
				mpom=rbind(matrix(rep(1000,2*nrow(int3)),ncol=nrow(int3)),t(as.matrix(int3[,1:3])),matrix(rep(0,2*nrow(int3)),ncol=nrow(int3)))
				res.int3=cbind(apply(t(mpom),1,intpaste),int3[,4])
				}else{
				mpom=rbind(matrix(rep(1000,2*nrow(int3)),ncol=nrow(int3)),as.matrix(int3[,1:3]),matrix(rep(0,2*nrow(int3)),ncol=nrow(int3)))
				res.int3=cbind(apply(t(mpom),1,intpaste),int3[,4])
				}
				}
if(!is.matrix(int3)){
				mpom=c(1000,1000,int3[1:3],rep(0,2))	
				res.int3=c(intpaste(mpom),int3[4])
				res.int3=t(as.matrix(res.int3))
				}
				}
res.int4=NULL
if(!is.null(int4)){
				if(is.matrix(int4)){
				if(nrow(int4)>1){
				mpom=rbind(matrix(rep(1000,3*nrow(int4)),ncol=nrow(int4)),t(as.matrix(int4[,1:4])))
				res.int4=cbind(apply(t(mpom),1,intpaste),int4[,5])
				}else{
				mpom=rbind(matrix(rep(1000,3*nrow(int4)),ncol=nrow(int4)),as.matrix(int4[,1:4]))
				res.int4=cbind(apply(t(mpom),1,intpaste),int4[,5])
				}
				}
if(!is.matrix(int4)){
				mpom=c(1000,1000,1000,int4[1:4])		
				res.int4=c(intpaste(mpom),int4[5])
				res.int4=t(as.matrix(res.int4))
				}
				}

	res.add.p=NULL
	res.add.p=res.add#posteriors for main effects
#########################################
write.table(rbind(res.add.p,res.int2,res.int3,res.int4), file = outFile1, append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names =FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}
###########################################
res.addp=NULL
if(length(which(res.add.p[,2]>=0.5))>0) res.addp=res.add.p[which(res.add.p[,2]>=0.5),]
resint2=NULL
if(length(which(as.numeric(res.int2[,2])>=0.5))>0) resint2=res.int2[which(as.numeric(res.int2[,2])>=0.5),]
resint3=NULL
if(length(which(as.numeric(res.int3[,2])>=0.5))>0) resint3=res.int3[which(as.numeric(res.int3[,2])>=0.5),]
resint4=NULL
if(length(which(as.numeric(res.int3[,2])>=0.5))>0) resint4=res.int4[which(as.numeric(res.int4[,2])>=0.5),]

TabAdd=NULL
if(!is.null(res.addp)){
					if(nrow(res.addp)>0) TabAdd=cbind(res.addp[,1],rep(1,nrow(res.addp)),res.addp[,2])
					if(nrow(res.addp)==0){res.addp=t(as.matrix(res.addp))
										  TabAdd=cbind(res.addp[,1],1,res.addp[,2])
										}
						}
TabInt=rbind(resint2,resint3,resint4)
if(!is.null(TabInt)) TabInt=cbind(TabInt[,1],rep(1,nrow(TabInt)),TabInt[,2])

tpfp.ILA=tp.fp.ILA(true.int,tab.add.p=TabAdd,tab.int.p=TabInt,thres=0.15)
###########################################
TP.add.ila<-tpfp.ILA$TP.add
TP.int.ila<-tpfp.ILA$TP.int
FP.add.ila<-tpfp.ILA$FP.add
FP.int.ila<-tpfp.ILA$FP.int
FDR.a.ila<-tpfp.ILA$FDR.a
FDR.i.ila<-tpfp.ILA$FDR.i
FDR.ila<-tpfp.ILA$FDR


# ## now calculate posteriors for a specific trees: ETA
###########################################
outFile2=paste(outFile,"_ETA.txt",sep="")
Addind=lapply(T,Int.tree.ind,ntrs=ntrs,nlvs=nlvs,lev=0)
indT.a=which(as.numeric(unlist(lapply(Addind,length)))>0)# which models contain at least one add effect
PIa=vector("list",length(T))
if(length(indT.a)>0){
for(ww in 1:length(indT.a))
PIa[[indT.a[[ww]]]]=take.adds(T[[indT.a[[ww]]]],ind=as.numeric(unlist(Addind[indT.a[[ww]]])))
}
alladds=as.numeric(unlist(PIa))
###########################################
# list of all different additive effects 
unadds=unique(sort(abs(alladds)))
#########################################
# take an additive effect , take the prior of all models that contain it, take corresponding score and sum the prior times the new score
Ia<-NULL
if(length(unadds)>0){
Ia=vector("list",length(unadds))
for( i in 1:length(unadds)){
Ij=NULL# models which contain a particular add effect
for( j in 1:length(PIa)) 
if((!is.na(match(unadds[i],PIa[[j]])))||(!is.na(match(-unadds[i],PIa[[j]]))))# ith additive effect appears in jth model
Ij<-c(Ij,j)
Ia[[i]]=Ij# save indices of such models 
}
#########################################
}
## data frame saved to file 
addposts=NULL
res.add=NULL
if(length(Ia)>0){
for(i in 1:length(Ia))
addposts[i]=sum(priors[Ia[[i]]]*Enewscores[Ia[[i]]])/sum(priors*Enewscores)
foo=cbind(unadds,round(addposts,digits=6))
if(length(which(foo[,2]>0))>0){
res.add=foo[order(foo[,2],foo[,1],decreasing=TRUE),]
if(is.vector(res.add)) res.add=t(as.matrix(res.add))
#res.add=res.add[which(res.add[,2]>0),]
}
}
if(!is.null(res.add)){
write.table("# additive effects & posterior probability ", outFile2,append = TRUE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
if(is.matrix(res.add))			
write.table(res.add, outFile2,append = TRUE, quote = FALSE, sep = "\t ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
			
if(!is.matrix(res.add))			
write.table(t(as.matrix(res.add)), outFile2,append = TRUE, quote = FALSE, sep = "\t ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
			}
#########################################
###### interactions 
Iind=lapply(T,Int.tree.ind,ntrs=ntrs,nlvs=nlvs,lev=1)
PIl=vector("list",length(Iind))
for(ww in 1:length(Iind)) PIl[[ww]]=take.ints(mat=T[[ww]],ind=Iind[[ww]]) 
#PIl= mapply(take.ints,T,ind=Iind)# all trees with interactions
PIi=do.call(rbind, PIl)#unlisted- as a matrix with trees in rows
unints=NULL
intposts=NULL
foo=NULL
sig.int=NULL
res.int=NULL
sigint=NULL
if(length(PIi)>0)
{
PIm=t(apply(PIi,1,tree.merge))# posortuj zawartosc drzew
m.ind=model.ind(nmtrs)
if(nrow(PIm)>2){PI=count.rows(PIm)
				Iint=PI[,2]
				Iint=lapply(Iint,function(i) m.ind[i])
				unints=PI[,c(3:ncol(PI))]
				}else{
				PIm=matrix(PIm,ncol=7)
				PI=cbind(c(1:nrow(PIm)),c(1:nrow(PIm)),PIm)
				Iint=PI[,2]
				Iint=lapply(Iint,function(i) m.ind[i])
				}

# #########################################
# # take an interaction , take the prior of all models that contain it, take corresponding score and sum the prior times the new score
#########################################
## data frame saved to file 
intposts=NULL
foo=NULL
sig.int=NULL
res.int=NULL
sigint=NULL
if(length(Iint)>0){
for(i in 1:length(Iint))
intposts[i]=sum(priors[Iint[[i]]]*Enewscores[Iint[[i]]])/sum(priors*Enewscores)
intposts=pmin(1,intposts)
if(!is.null(unints)){
if(nrow(unints)>1){
foo=cbind(unints,round(intposts,digits=6))
res.int=foo[order(foo[,(ncol(foo))],decreasing=TRUE),]
sig.int=NULL
if(length(which(res.int[,ncol(foo)]>0))>0)
	{
	is=which(res.int[,ncol(foo)]>0)
	res.int=res.int[is,]
	}else{
			res.int=res.int[-c(1:nrow(res.int)),]
			}
if(nrow(res.int)>0){
					sigint<-res.int[,-ncol(res.int)]
					C1=apply(sigint,1,intpaste)
					sig.int=cbind(C1,res.int[,ncol(res.int)])
					}
					}else{
							foo=c(unints,round(intposts,digits=6))
							res.int=foo#[order(foo[,(ncol(foo))],decreasing=TRUE),]
							sig.int=NULL
							sigint<-res.int[-length(res.int)]
							C1=intpaste(sigint)
							sig.int=c(C1,res.int[length(res.int)])
							}
}
}

if(!is.null(sig.int)){
write.table("# interactions  & posterior probability ", outFile2,append = TRUE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
#if(is.matrix(sig.int))			
write.table(sig.int, outFile2,append = TRUE, quote = FALSE, sep = "\t ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
			}
#########################################
}
res.add.ETA<-NULL
if(is.vector(res.add)) res.add=t(as.matrix(res.add))
if(is.vector(sig.int)) sig.int=t(as.matrix(sig.int))

if(length(which(res.add[,2]>=0.5))>0) res.add.ETA<-res.add[which(res.add[,2]>=0.5),]
res.int.ETA<-NULL
if(length(which(as.numeric(sig.int[,2])>=0.5))>0) res.int.ETA<-sig.int[which(as.numeric(sig.int[,2])>=0.5),]
tab.add.p.eta=NULL

if(!is.null(res.add.ETA)) if(is.matrix(res.add.ETA)){ 
											tab.add.p.eta=cbind(res.add.ETA[,1],rep(1,nrow(res.add.ETA)),res.add.ETA[,2])
											}else{
											tab.add.p.eta=c(res.add.ETA[1],1,res.add.ETA[2])
											tab.add.p.eta=t(as.matrix(tab.add.p.eta))}
tab.int.p.eta=NULL
if(!is.null(res.int.ETA)) if(is.matrix(res.int.ETA)){ 
											tab.int.p.eta=cbind(res.int.ETA[,1],rep(1,nrow(res.int.ETA)),res.int.ETA[,2])
											}else{
											tab.int.p.eta=c(res.int.ETA[1],1,res.int.ETA[2])
											tab.int.p.eta=t(as.matrix(tab.int.p.eta))
											}


tpfp.ETA=tp.fp.ETA(true.int,tab.add.p=tab.add.p.eta,tab.int.p=tab.int.p.eta,thres=0.15)

TP.add.eta<-tpfp.ETA$TP.add
TP.int.eta<-tpfp.ETA$TP.int
FP.add.eta<-tpfp.ETA$FP.add
FP.int.eta<-tpfp.ETA$FP.int
FDR.a.eta<-tpfp.ETA$FDR.a
FDR.i.eta<-tpfp.ETA$FDR.i
FDR.eta<-tpfp.ETA$FDR

######## cumulative for trees 
## checking for all expressions containing variables from the expression ( or their compliments) and the same conector

resExpr=NULL

if(!is.null(sigint)){
Ppost=res.int[,ncol(res.int)]
zm=apply(sigint,1,tvariables.from.model)
# ssigint=t(apply(sigint,1,tree.merge))# sorted trees 
# cr=count.rows(ssigint)
# join trees which contain exactly the same variables and operators
Expr=NULL
sigint1=sigint
#for(ex in 1:length(zm)){
ex=1
# ExprList=NULL
while(ex<(nrow(sigint))){
#if(compare.expr(Expr,sigint[ex,])){ex=ex+2}else{ex=ex+1}
Ni=ex
P=Ppost[ex]
UnL.int=NULL
TreeList=NULL
# PostList=vector("list",nrow(sigint))
for(ni in (ex+1):nrow(sigint1)){
if(length(which(zm[[ex]]%in%zm[[ni]]))==length(zm[[ex]])){
Ni=c(Ni,ni)
P=c(P,Ppost[ni])
}
}
UnL.int=sigint1[ex,]
TreeList=sigint1[Ni,]
PostList=P
################################
## porzadkowanie drzew
if(nrow(TreeList)>0){
untv=TreeList#t(apply(TreeList,1,unique))
if(nrow(TreeList)>1){
# if(!is.list(untv)){
# untv=t(untv)
 # untv=split(untv, rep(1:ncol(untv), each = nrow(untv)))
# }
tv=untv#t(apply(untv,1,abs))
stv=tv#t(apply(tv,1,sort,decreasing=TRUE))
#stv# uzupelnic zerami do drzewa
#n0ind=lapply(stv,function(i) which(i>0))
#lapply(stv,stv[n0ind])
# remove 0s
stvn0=apply(stv,1,rem.null)
}# uzupelnic zerami do drzewa
if(nrow(TreeList)==1){
untv=TreeList
tv=TreeList
stv=tv#sort(tv,decreasing=TRUE)
}
}
lmv=NULL
mv=NULL
if(nrow(stv)>1){
#Lnull=lapply(stv,length)
mv=apply(stv,1,variables.from.model)
lmv=lapply(mv,length)#liczba zmiennych 
Lop=lapply(lmv, function(i) i-1)
sumV=lapply(seq_along(Lop),function(i)
         unlist(Lop[[i]])+unlist(lmv[[i]]))# suma dlugosci operatorow i zmiennych 
Lnull=lapply(seq_along(sumV),function(i)
         7-unlist(sumV[[i]]))# liczba zer do dopisania 
}

if(nrow(stv)==1){
#Lnull=lapply(stv,length)
mv=variables.from.model(stv)
lmv=length(mv)#liczba zmiennych 
Lop=lmv-1
sumV=lmv+Lop#lapply(seq_along(Lop),function(i)
      #   unlist(Lop[[i]])+unlist(lmv[[i]]))# suma dlugosci operatorow i zmiennych 
Lnull=7-sumV
stvn0=rem.null(stv)
}

newTrees=NULL
Cr=NULL
PostTab=NULL
if(nrow(stv)==1){newTrees=stv}
if(nrow(stv)>1){
newTrees=stv
#t(mapply(make.tree,lop=Lop,lnull=Lnull,n0v=stvn0,sumv=sumV,lmv=lmv))
}

if(nrow(newTrees)>2){
#if(nrow(unique(newTrees))>2){
Cr=count.rows(newTrees)
eq.trees=Cr[,2]#which are the same
l.eq.trees=lapply(eq.trees,length)
list.ind=which(as.numeric(unlist(l.eq.trees))>1)# which are repeated
if(length(list.ind)>0){
list.ind.p=eq.trees[list.ind]# go to list index and sum up posteriors with indices from this list inex
Post.list=lapply(list.ind.p,function(i) PostList[i])
#take posteriors of elements from this list
nTsumPost=unlist(as.numeric(lapply(Post.list,sum)))# sum up posteriors in elements of a list
#sum(PostList[eq.trees[list.ind]])
lmv=NULL
mv=NULL
PostTab=NULL
PP=NULL
PP[list.ind]=nTsumPost
rem.ind=c(1:nrow(newTrees))[-unlist(eq.trees[list.ind])]
newTrees1=unique(newTrees)
n.list.ind=c(1:nrow(newTrees1))[-list.ind]
PP[n.list.ind]=PostList[rem.ind]
PostTab=cbind(newTrees1,PP)
mv=apply(newTrees1,1,variables.from.model)
lmv=as.numeric(unlist(lapply(mv,length)))#liczba zmiennych 
}
if(length(list.ind)==0){
newTrees1=newTrees
PostTab=cbind(newTrees,P)
mv=apply(newTrees,1,variables.from.model)
lmv=as.numeric(unlist(lapply(mv,length)))#liczba z
}
}
if(nrow(newTrees)==2){
if(nrow(unique(newTrees))==2){
nTposts=P
nTsumpost=P
newTrees1=newTrees
PostTab=cbind(newTrees,nTsumpost)
# mv=variables.from.model(newTrees)
# lmv=length(mv)
}
}


if(nrow(unique(newTrees))==1){
newTrees1=unique(newTrees)
nTposts=sum(P)
nTsumpost=sum(P)
PostTab=cbind(newTrees,nTsumpost)
mv=variables.from.model(newTrees)
lmv=length(mv)
}

newTrees=newTrees1

## now check if the list of trees contains specific Logic expression within the tree 

# LI=apply(newTrees1,1,tvariables.from.model)
# Lln=lapply(LI,length)
b.ex=newTrees[1,]
p.ex=PostTab[1,8]
Which.cont=1
p.i.2=0
p.i.3=0
p.i.4=0
Lln=NULL
Lln=length(as.numeric(tvariables.from.model(b.ex)))
if(Lln==2) p.i.2=PostTab[1,8]
if(Lln==3) p.i.3=PostTab[1,8]
if(Lln==4) p.i.4=PostTab[1,8]
if(nrow(newTrees)>1){
for(ii in 2:nrow(newTrees)){
op.i=which(sapply(newTrees[ii,],is.operator))# indeks operatora
nop=length(op.i)# number of operators
op=as.numeric(newTrees[ii,which(newTrees[ii,]>=1000)])# operatory
if((length(unique(op))==1)&&(unique(op)==1000)){
# expression contains the specific b.ex
p.ex=p.ex+PostTab[ii,8]
Which.cont=c(Which.cont,ii)
Lli=length(as.numeric(tvariables.from.model(newTrees[ii,])))
if(Lli==2) p.i.2=p.i.2+PostTab[ii,8]
if(Lli==3) p.i.3=p.i.3+PostTab[ii,8]
if(Lli==4) p.i.4=p.i.4+PostTab[ii,8]
}
if(length(unique(op))>1){
if(length(op)==2){
if((op[op.i[1]]==2000)&&(op[op.i[2]]==1000)){
poz=NULL
for(hh in 1:length(zm[[ex]])) poz[hh]=which(newTrees[ii,]==zm[[ex]][hh])
if(sum(poz%in%c(4,5))==2){
p.ex=p.ex+PostTab[ii,8]
Which.cont=c(Which.cont,ii)
p.i.3=p.i.3+PostTab[ii,8]
}
}
}

if(length(op)==3){
if((op[op.i[1]]==2000)&&(op[op.i[2]]==1000)&&(op[op.i[3]]==1000)){
poz=NULL
for(hh in 1:length(zm[[ex]])) poz[hh]=which(newTrees[ii,]==zm[[ex]][hh])
if(sum(poz%in%c(4,5))==2){
p.ex=p.ex+PostTab[ii,8]
Which.cont=c(Which.cont,ii)
p.i.4=p.i.4+PostTab[ii,8]
}
if(sum(poz%in%c(6,7))==2){
p.ex=p.ex+PostTab[ii,8]
Which.cont=c(Which.cont,ii)
p.i.4=p.i.4+PostTab[ii,8]
}
}
}

}
}
}

fin.tab=c(intpaste(b.ex),p.i.2,p.i.3,p.i.4,p.ex)
Expr=sigint[ex,]
ResExpr=fin.tab
resExpr=rbind(resExpr,ResExpr)
ex=ex+1
}
}
TP.add.ceta=TP.add.ila



 
resExpr.p=NULL
tab.int.p.ceta=NULL
if(length(which(as.numeric(resExpr[,5])>=0.5))>1){ 
					resExpr.p=resExpr[which(as.numeric(resExpr[,5])>=0.5),]
					tab.int.p.ceta=cbind(resExpr.p[,1],rep(1,nrow(resExpr.p)),resExpr.p[,5])
}
if(length(which(as.numeric(resExpr[,5])>=0.5))==1){ 
					resExpr.p=resExpr[which(as.numeric(resExpr[,5])>=0.5),]
					tab.int.p.ceta=c(resExpr.p[1],1,resExpr.p[5])
					tab.int.p.ceta=t(as.matrix(tab.int.p.ceta))
					}


#### tp.fp.CETA function
tpfp.ceta=tp.fp.CETA(true.int,tab.add.p=TabAdd,tab.int.p=tab.int.p.ceta,thres=0.15)

TP.add.ceta<-tpfp.ceta$TP.add
TP.int.ceta<-tpfp.ceta$TP.int
FP.add.ceta<-tpfp.ceta$FP.add
FP.int.ceta<-tpfp.ceta$FP.int
FDR.a.ceta<-tpfp.ceta$FDR.a
FDR.i.ceta<-tpfp.ceta$FDR.i
FDR.ceta<-tpfp.ceta$FDR

write.table(res.add.p,file=paste(outFile3,"_add.txt"),append = TRUE, quote = FALSE, sep = "\t ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")		
write.table(resExpr,file=paste(outFile3,"_Int.txt"),append = TRUE, quote = FALSE, sep = "\t ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
######################
return(list(res.add=res.add.p,int2=res.int2,int3=res.int3,int4=res.int4,res.add.T=res.add,sig.int.T=sig.int,res.expr=resExpr,TP.add.ila=TP.add.ila,TP.int.ila=TP.int.ila,FP.add.ila=FP.add.ila,FP.int.ila=FP.int.ila,FDR.a.ila=FDR.a.ila,FDR.i.ila=FDR.i.ila,FDR.ila=FDR.ila,TP.add.eta=TP.add.eta,TP.int.eta=TP.int.eta,FP.add.eta=FP.add.eta,FP.int.eta=FP.int.eta,FDR.a.eta=FDR.a.eta,FDR.i.eta=FDR.i.eta,FDR.eta=FDR.eta,TP.add.ceta=TP.add.ceta,TP.int.ceta=TP.int.ceta,FP.add.ceta=FP.add.ceta,FP.int.ceta=FP.int.ceta,FDR.a.ceta=FDR.a.ceta,FDR.i.ceta=FDR.i.ceta,FDR.ceta=FDR.ceta))
}
#####################################################



tp.fp.ETA=function(true.int,tab.add.p,tab.int.p,thres=0.15)
{
##***********************************************************
l.fp.add<-0
fp.add<-NULL
fr.fp.add<-NULL
p.fp.add<-NULL
l.fp.int<-0
fp.int<-NULL
fr.fp.int<-NULL
p.fp.int<-NULL

l.tp.add<-0
tp.add<-NULL
fr.tp.add<-NULL
p.tp.add<-NULL
l.tp.int<-0
tp.int<-NULL
fr.tp.int<-NULL
p.tp.int<-NULL

if(!is.null(true.int)){
 write.table(true.int, file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  true.int1=read.table(file = "M.int.txt")
  fn <- "M.int.txt"
  if (file.exists(fn)) file.remove(fn)
  }

 
where.int=which(unlist(lapply(apply(true.int1,1,tree.operators),length))>0)#which of true int are really interactions
Adds=true.int1[-where.int,1]

#index of true additive effect in tab.add.p (row)
ind.Adds=rep(0,length(Adds))
for(jj in 1:length(Adds)) if(length(which(tab.add.p[,1]==Adds[jj]))>0) ind.Adds[jj]=which(tab.add.p[,1]==Adds[jj])
fr.tp.add=tab.add.p[ind.Adds,2]
p.tp.add=tab.add.p[ind.Adds,3]
l.tp.add=length(which(ind.Adds>0))
tp.add=tab.add.p[ind.Adds,1]

fr.fp.add=tab.add.p[-ind.Adds,2]
p.fp.add=tab.add.p[-ind.Adds,3]
l.fp.add=ifelse(!is.null(nrow(tab.add.p)),nrow(tab.add.p)-l.tp.add,0)
l.fp.add
fp.add=tab.add.p[-ind.Adds,1]

if(!is.null(tab.int.p))
if(nrow(tab.int.p)>0){
write.table(tab.int.p[,1], "int-p.txt",append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
intP=read.table("int-p.txt")
fn <- "int-p.txt"
if (file.exists(fn)) file.remove(fn)
int.P=t(apply(intP,1,tree.merge))
if(length(apply(int.P,1,intpaste))>1){
tab.int.P=cbind(apply(int.P,1,intpaste),tab.int.p[,c(2,3)])}else{tab.int.P=t(as.matrix(c(apply(int.P,1,intpaste),tab.int.p[,c(2,3)])))}
tab.int.p<-tab.int.P
}

Ints=true.int1[where.int,]

#index of true interaction effect in tab.int.p (row)
ind.Ints=rep(0,nrow(Ints))
for(jj in 1:nrow(Ints)) ind.Ints[jj]=match(true.int[jj],tab.int.p[,1],nomatch=FALSE)
if(length(which(ind.Ints>0))>0){
fr.tp.int=as.numeric(tab.int.p[ind.Ints,2])
p.tp.int=as.numeric(tab.int.p[ind.Ints,3])
l.tp.int=length(which(ind.Ints>0))
tp.int=tab.int.p[ind.Ints,1]
############ FP
fr.fp.int=as.numeric(tab.int.p[-ind.Ints,2])
p.fp.int=as.numeric(tab.int.p[-ind.Ints,3])
l.fp.int=nrow(tab.int.p)-l.tp.int
l.fp.int
fp.int=tab.int.p[-ind.Ints,1]	
}else{
	fr.fp.int=as.numeric(tab.int.p[,2])
	p.fp.int=as.numeric(tab.int.p[,3])
	l.fp.int=nrow(tab.int.p)
	l.fp.int
	fp.int=tab.int.p[,1]}
	
p.fp<-NULL
p.fp<-c(p.fp.add,p.fp.int)
fr.fp<-NULL
fr.fp<-c(fr.fp.add,fr.fp.int)
m.fp.int<-NULL
l.fp<-c(l.fp.add,l.fp.int)
if(length(fp.int)>0){
write.table(fp.int, file = paste("fpint.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
fp.int.v<-read.table(file = paste("fpint.txt",sep=""))
fn <- "fpint.txt"
if (file.exists(fn)) file.remove(fn)
m.fp.int<-fp.int.v
}

##*******
# ## POPRAWNE "KLASYFIKACJE" 
##*****

l.tp<-NULL
l.tp<-c(l.tp.add,l.tp.int)

tp<-NULL
tp<-c(tp.add,tp.int)

fr.tp<-NULL
fr.tp<-c(fr.tp.add,fr.tp.int)

p.tp<-NULL
p.tp<-c(p.tp.add,p.tp.int)

m.tp.add<-NULL
if(!is.null(tp.add)){
tp.add<-as.numeric(tp.add)
m.tp.add=tp.add
}

m.tp.int<-NULL
if(length(tp.int)>0){
  write.table(tp.int, file= "tpint.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  tp.int.v<-read.table(file = "tpint.txt")
  fn <- "tpint.txt"
  if (file.exists(fn)) file.remove(fn)
  m.tp.int<-tp.int.v
}
FDR=rep(0,2)
if(length(which((l.fp+l.tp)==0))>0){ 
FDR[which((l.fp+l.tp)==0)]=0
FDR[-which((l.fp+l.tp)==0)]=l.fp[-which((l.fp+l.tp)==0)]/(l.fp[-which((l.fp+l.tp)==0)]+l.tp[-which((l.fp+l.tp)==0)])
}else{FDR=l.fp/(l.fp+l.tp)}
FDR

TP.add=NULL
if(!is.null(tp.add)) TP.add=cbind(tp.add,fr.tp.add,p.tp.add)
TP.int=NULL
if(!is.null(m.tp.int)) TP.int=cbind(m.tp.int,fr.tp.int,p.tp.int)
FP.add=NULL
if(!is.null(fp.add)) FP.add=cbind(fp.add,fr.fp.add,p.fp.add)
FP.int=NULL
if(!is.null(m.fp.int)) FP.int=cbind(m.fp.int,fr.fp.int,p.fp.int)
colnames(TP.add)=colnames(TP.int)=colnames(FP.add)=colnames(FP.int)=NULL

s.a.tp=NULL
s.i.tp=NULL
s.a.fp=NULL
s.i.fp=NULL
if(is.matrix(TP.add)){s.a.tp=sum(TP.add[,2])}else{
												  if(length(TP.add)>0){s.a.tp=sum(TP.add[2])
																		}else{s.a.tp=0}
												    }
if(is.matrix(TP.int)){s.i.tp=sum(TP.int[,8])}else{
													if(length(TP.int)>0){s.i.tp=sum(TP.int[8])
																			}else{s.i.tp=0}
													}
if(is.matrix(FP.add)){s.a.fp=sum(FP.add[,2])}else{
													if(length(FP.add)>0){s.a.fp=sum(FP.add[2])
																			}else{s.a.fp=0}
													}
if(is.matrix(FP.int)){s.i.fp=sum(FP.int[,8])}else{
													if(length(FP.int)>0){s.i.fp=sum(FP.int[8])
																			}else{s.i.fp=0}
													}

l.tp=c(s.a.tp,s.i.tp)
l.fp=c(s.a.fp,s.i.fp)

ltp=l.tp
lfp=l.fp

FDR=rep(0,2)
if(length(which((l.fp+l.tp)==0))>0){ 
FDR[which((l.fp+l.tp)==0)]=0
FDR[-which((l.fp+l.tp)==0)]=l.fp[-which((l.fp+l.tp)==0)]/(l.fp[-which((l.fp+l.tp)==0)]+l.tp[-which((l.fp+l.tp)==0)])
}else{FDR=l.fp/(l.fp+l.tp)}
FDR

FDR.total=0
if(length(which(sum(lfp)+sum(ltp)==0))>0){ 
FDR.total[which((sum(lfp)+sum(ltp))==0)]=0
FDR.total[-which((sum(lfp)+sum(ltp))==0)]=sum(lfp)[-which((sum(lfp)+sum(ltp))==0)]/(sum(lfp)[-which((sum(lfp)+sum(ltp))==0)]+sum(ltp)[-which((sum(lfp)+sum(ltp))==0)])
}else{FDR.total=sum(lfp)/(sum(lfp)+sum(ltp))}
FDR.total

# ## FP depending on the threshold
# FP.int.t=NULL
# FP.add.t=NULL
 # FP.int.t=FP.int[which(FP.int[,8]/n.sym>=thres),]
 # FP.add.t=FP.add[which(FP.add[,2]/n.sym>=thres),]
# TP.int.t=TP.int
# TP.add.t=TP.add
# colnames(TP.add)=colnames(TP.int)=colnames(FP.add)=colnames(FP.int)=NULL

# s.a.tp=NULL
# s.i.tp=NULL
# s.a.fp=NULL
# s.i.fp=NULL
# if(is.matrix(TP.add.t)){s.a.tp=sum(TP.add.t[,2])}else{
												  # if((length(TP.add.t)>0)&&(nrow(TP.add.t)>0)){s.a.tp=sum(TP.add.t[2])
																		# }else{s.a.tp=0}
												    # }
# if(is.matrix(TP.int.t)){s.i.tp=sum(TP.int.t[,8])}else{
													# if((length(TP.int.t)>0)&&(nrow(TP.int.t)>0)){s.i.tp=sum(TP.int.t[8])
																			# }else{s.i.tp=0}
													# }
# if(is.matrix(FP.add.t)){s.a.fp=sum(FP.add.t[,2])}else{
													# if((length(FP.add.t)>0)&&(nrow(FP.add.t)>0)){s.a.fp=sum(FP.add.t[2])
																			# }else{s.a.fp=0}
													# }
# if(is.matrix(FP.int.t)){s.i.fp=sum(FP.int.t[,8])}else{
													# if((length(FP.int.t)>0)&&(nrow(FP.int.t)>0)){s.i.fp=sum(FP.int.t[8])
																			# }else{s.i.fp=0}
													# }

# l.tp=c(s.a.tp,s.i.tp)
# l.fp=c(s.a.fp,s.i.fp)


# FDR.t=rep(0,2)
# if(length(which((l.fp+l.tp)==0))==1){ 
# #FDR.t[which((l.fp+l.tp)==0)]=0
# FDR.t[-which((l.fp+l.tp)==0)]=l.fp[-which((l.fp+l.tp)==0)]/(l.fp[-which((l.fp+l.tp)==0)]+l.tp[-which((l.fp+l.tp)==0)])
# }#else{FDR.t=l.fp/(l.fp+l.tp)}
# FDR.t

return(list(TP.add=TP.add,TP.int=TP.int,FP.add=FP.add,FP.int=FP.int,FDR.a=FDR[1],FDR.i=FDR[2],FDR=FDR.total))
#,TP.add.t=TP.add.t,TP.int.t=TP.int.t,FP.add.t=FP.add.t,FP.int.t=FP.int.t,FDR.a.t=FDR.t[1],FDR.i.t=FDR.t[2],FDR.t=sum(l.fp)/(sum(l.fp)+sum(l.tp))))
}

##

tp.fp.ETA.plot=function(true.int,tab.add.p,tab.int.p,max.thres=max.thres,OutFile=OutFile){
FDR=NULL
FDR.a.t=NULL
FDR.i.t=NULL
for(hh in 0:10){
res=tp.fp.ETA(true.int,tab.add.p,tab.int.p,thres=max.thres*hh*0.1)
FDR.a.t=c(FDR.a.t,res$FDR.a.t)
FDR.i.t=c(FDR.i.t,res$FDR.i.t)
}
jpeg(file=file.path(OutFile,"_ETA_FDR-t.jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
plot(c(0:10)*0.1*max.thres,FDR.a.t,ylim=c(0,1),xlab="Threshold",ylab="FDR", main="FDR for ETA",type="l",col="red",lty=2,lwd=2)
points(c(0:10)*0.1*max.thres,FDR.i.t,type="l",col="blue",lty=1,lwd=2)
legend(x=0.3,y=0.2,c("FDR add", "FDR int"),lty=c(2,1),col=c("red","blue"),inset = .01,cex=0.6,bty="n", title.adj=.55)
dev.off()
return(list(FDR.a.t=FDR.a.t,FDR.i.t=FDR.i.t))
}


all.subtrees=function(vec)
{
vec=as.numeric(vec)
if(is.leaf(vec[1])){
					subtrees=t(as.matrix(as.numeric(vec)))
					#colnames(subtrees)=NULL
					}else{
						# leaves=variables.from.model(vec)
						leaves=sign(vec[which(unlist(lapply(vec,is.leaf)))])*variables.from.model(vec)
						subtrees=t(rbind(leaves,matrix(rep(0,length(leaves)*6),ncol=length(leaves))))
						colnames(subtrees)=NULL
						intlev=length(tree.operators(vec))+1#level of interaction
						diff.op=length(unique(tree.operators(vec)))# no of different operators (1 or 2)
						if(intlev==2){
										subtrees=rbind(subtrees,vec)
										rownames(subtrees)=NULL
										colnames(subtrees)=NULL
										}
						if(intlev==3){
										if(diff.op==1){
													subpairs=subvectors(leaves)$subpairs
													sub2=cbind(rep(vec[1],nrow(subpairs)),subpairs[,2:1],matrix(rep(0,nrow(subpairs)*4),ncol=4))
													colnames(sub2)=NULL
													subtrees=rbind(subtrees,sub2)
													rownames(subtrees)=NULL
													colnames(subtrees)=NULL	
														}
										if(diff.op==2){
													sub2=rbind(c(vec[c(2,4,5)],rep(0,4)),
													c(vec[1],sort(vec[3:4],decreasing=TRUE),rep(0,4)),
													c(vec[1],sort(vec[c(3,5)],decreasing=TRUE),rep(0,4)))
													colnames(sub2)=NULL	
													subtrees=rbind(subtrees,sub2)
													rownames(subtrees)=NULL
													colnames(subtrees)=NULL
													}
										}
						if(intlev==4){
										if(diff.op==1){
													subpairs=subvectors(leaves)$subpairs
													sub2=cbind(rep(vec[1],nrow(subpairs)),subpairs[,2:1],matrix(rep(0,nrow(subpairs)*4),ncol=4))
													colnames(sub2)=NULL
													subtriples=subvectors(leaves)$subtriples
													sub3=cbind(matrix(rep(vec[1],nrow(subtriples)*2),ncol=2),t(apply(subtriples,1,sort,decreasing=TRUE)),matrix(rep(0,nrow(subtriples)*2),ncol=2))
													colnames(sub3)=NULL
													subtrees=rbind(subtrees,sub2,sub3)
													rownames(subtrees)=NULL
													colnames(subtrees)=NULL		
														}
										if(diff.op==2){
													if(which(duplicated(vec[1:3]))==3){# it means 2000 1000 1000 tree
														m1=c(rep(1000,2),rep(2000,4))
														m21=rbind(vec[6:7],vec[4:5])
														m22=matrix(c(sort(vec[c(5,6)],decreasing=TRUE),sort(vec[c(4,7)],decreasing=TRUE),sort(vec[c(4,6)],decreasing=TRUE),sort(vec[c(5,7)],decreasing=TRUE)),ncol=2)
														m2=rbind(m21,m22)
														m3=matrix(rep(0,4*length(m1)),ncol=4)
														sub2=cbind(m1,m2,m3)
														colnames(sub2)=NULL
														subtrees=rbind(subtrees,sub2)
														m1=matrix(rep(c(2000,1000),each=4),ncol=2)
														m2=rbind(vec[5:7],vec[c(4,6,7)],vec[4:6],vec[c(4,5,7)])
														m3=matrix(rep(0,8),ncol=2)
														sub3=cbind(m1,m2,m3)
														colnames(sub3)=NULL
														subtrees=rbind(subtrees,sub3)
														rownames(subtrees)=NULL
														colnames(subtrees)=NULL
														}
													if(which(duplicated(vec[1:3]))==2){# it means 2000 2000 1000 tree
														m1=c(1000,rep(2000,5))
														m21=rbind(vec[6:7],vec[4:5])
														m22=matrix(c(sort(vec[c(4,6)],decreasing=TRUE),sort(vec[c(4,7)],decreasing=TRUE),sort(vec[c(5,6)],decreasing=TRUE),sort(vec[c(5,7)],decreasing=TRUE)),ncol=2)
														m2=rbind(m21,m22)
														m3=matrix(rep(0,4*length(m1)),ncol=4)
														sub2=cbind(m1,m2,m3)
														colnames(sub2)=NULL
														subtrees=rbind(subtrees,sub2)
														m1=matrix(c(rep(2000,4),rep(1000,2),rep(2000,2)),ncol=2)
														m2=rbind(vec[c(6,7,5)],vec[c(6,7,4)],vec[c(6,4,5)],vec[c(7,4,5)])
														m3=matrix(rep(0,8),ncol=2)
														sub3=cbind(m1,m2,m3)
														colnames(sub3)=NULL
														subtrees=rbind(subtrees,sub3)
														rownames(subtrees)=NULL
														colnames(subtrees)=NULL
														}	
														
														}
										}
						}
subtrees=t(apply(subtrees,1,tree.merge))						
return(subtrees)
}

tp.fp.CETA=function(true.int,tab.add.p,tab.int.p,thres=0.15)
{
##***********************************************************
l.fp.add<-0
fp.add<-NULL
fr.fp.add<-NULL
p.fp.add<-NULL
l.fp.int<-0
fp.int<-NULL
fr.fp.int<-NULL
p.fp.int<-NULL

l.tp.add<-0
tp.add<-NULL
fr.tp.add<-NULL
p.tp.add<-NULL
l.tp.int<-0
tp.int<-NULL
fr.tp.int<-NULL
p.tp.int<-NULL

if(!is.null(true.int)){
 write.table(true.int, file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  true.int1=read.table(file = "M.int.txt")
  fn <- "M.int.txt"
  if (file.exists(fn)) file.remove(fn)
  }

 
where.int=which(unlist(lapply(apply(true.int1,1,tree.operators),length))>0)#which of true int are really interactions
#Adds=true.int1[-where.int,1]
res=vector("list",nrow(true.int1))
for(i in 1:nrow(true.int1))
res[[i]]=all.subtrees(true.int1[i,])
res

true.subtrees=res#apply(true.int1,1,all.subtrees)
true.subtrees
true.sub.add=NULL
for(kk in 1:length(true.subtrees)) true.sub.add=c(true.sub.add,true.subtrees[[kk]][which(rowSums(abs(true.subtrees[[kk]]))<1000),1])

Adds=true.sub.add
#index of true additive effect in tab.add.p (row)
ind.Adds=0
for(jj in 1:length(Adds))
				{
				if(length(which(tab.add.p[,1]==Adds[jj]))>0)
				ind.Adds<-c(ind.Adds,which(tab.add.p[,1]==Adds[jj]))
				}
fr.tp.add=tab.add.p[ind.Adds,2]
p.tp.add=tab.add.p[ind.Adds,3]
l.tp.add=length(which(ind.Adds>0))
tp.add=tab.add.p[ind.Adds,1]

fr.fp.add=tab.add.p[-ind.Adds,2]
p.fp.add=tab.add.p[-ind.Adds,3]
l.fp.add=nrow(tab.add.p)-l.tp.add
l.fp.add
fp.add=tab.add.p[-ind.Adds,1]


## 1) trees which coincide with a true tree --- procedure as in ETA 
int.P=NULL
if(!is.null(tab.int.p)){
if(nrow(tab.int.p)>0){
write.table(tab.int.p[,1], "int-p.txt",append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
intP=read.table("int-p.txt")
fn <- "int-p.txt"
if (file.exists(fn)) file.remove(fn)
int.P=t(apply(intP,1,tree.merge))
if(length(apply(int.P,1,intpaste))>1){
tab.int.P=cbind(apply(int.P,1,intpaste),tab.int.p[,c(2,3)])}else{tab.int.P=t(as.matrix(c(apply(int.P,1,intpaste),tab.int.p[,c(2,3)])))}
tab.int.p<-tab.int.P
}
}
Ints=true.int1[where.int,]

#index of true interaction effect in tab.int.p (row)
ind.Ints=rep(0,nrow(Ints))
for(jj in 1:nrow(Ints)) ind.Ints[jj]=match(true.int[jj],tab.int.p[,1],nomatch=FALSE)
fr.tp.int=as.numeric(tab.int.p[ind.Ints,2])
p.tp.int=as.numeric(tab.int.p[ind.Ints,3])
l.tp.int=length(which(ind.Ints>0))
tp.int=tab.int.p[ind.Ints,1]
ceta.tp.int=NULL
if(sum(ind.Ints)>0){if(length(which(ind.Ints>0))>1){
					ceta.tp.int=cbind(as.matrix(int.P[ind.Ints,]),matrix(as.numeric(tab.int.p[ind.Ints,2:3]),ncol=2))
													}else{
													ceta.tp.int=c(int.P[ind.Ints,],as.numeric(tab.int.p[ind.Ints,2:3]))
													ceta.tp.int=t(as.matrix(ceta.tp.int))
														}
					}
 

ceta.tp.int
int.PP=int.P
tab.int.pp=tab.int.p
if(sum(ind.Ints)>0){
int.PP=int.P[-ind.Ints,]
tab.int.pp=tab.int.p[-ind.Ints,]
}

if(!is.null(tab.int.pp)) if(!is.matrix(tab.int.pp)) tab.int.pp=t(as.matrix(tab.int.pp))
if(!is.null(int.PP)) if(!is.matrix(int.PP)) int.PP=t(as.matrix(int.PP))

ceta.fp.int=NULL
if(!is.null(tab.int.pp))
ceta.fp.int=cbind(int.PP,matrix(as.numeric(tab.int.pp[,2:3]),ncol=2))

## 2) subtrees of the true interactions
#true.subtrees=apply(true.int1,1,all.subtrees)
true.subtrees
ind.sub=NULL
if(!is.null(int.PP)){
for(gg in 1:length(true.int))
ind.sub=c(ind.sub,which(apply(int.PP,1,compare.expr, exprList=true.subtrees[[gg]])))##indices of rows from int.P which are subtrees of all true int
if(length(ind.sub)>0){
if(!is.null(tab.int.pp)) ceta.tp.int=rbind(ceta.tp.int,cbind(int.PP[ind.sub,],matrix(as.numeric(tab.int.pp[ind.sub,2:3]),ncol=2)))
int.PP=int.PP[-ind.sub,]
tab.int.pp=tab.int.pp[-ind.sub,]
if(!is.null(tab.int.pp)) ceta.fp.int=cbind(int.PP,matrix(as.numeric(tab.int.pp[,2:3]),ncol=2))
}
}

## 3) expressions which contain a true tree as a subtree
#subtrees of all expressions 
sub.All=NULL
if(!is.null(int.PP)){
if(nrow(int.PP)>0){ 
res=vector("list",nrow(int.PP))
for(i in 1:nrow(int.PP))
res[[i]]=all.subtrees(int.PP[i,])
res
sub.All=res
}
}

ind.c.sub=NULL #indices of rows from int.P which contain true int as a subtree 
if(!is.null(sub.All)){
for(gg in 1:length(sub.All)){
if(sum(apply(true.int1,1,compare.expr, exprList=sub.All[[gg]]))>0)
ind.c.sub=c(ind.c.sub,gg)
}
}

if(!is.null(ind.c.sub)){
if(length(ind.c.sub)>0)
{
if(length(ind.c.sub)==1){ceta.tp.int=rbind(ceta.tp.int,cbind(t(as.matrix(int.PP[ind.c.sub,])),matrix(as.numeric(tab.int.pp[ind.c.sub,2:3]),ncol=2)))}else{
ceta.tp.int=rbind(ceta.tp.int,cbind(int.PP[ind.c.sub,],matrix(as.numeric(tab.int.pp[ind.c.sub,2:3]),ncol=2)))}
int.PP=int.PP[-ind.c.sub,]
tab.int.pp=tab.int.pp[-ind.c.sub,]
}
}

if(!is.null(tab.int.pp)){
if(is.matrix(tab.int.pp)){
				ceta.fp.int=cbind(int.PP,matrix(as.numeric(tab.int.pp[,2:3]),ncol=2))
				}else{
				ceta.fp.int=c(int.PP,as.numeric(tab.int.pp[2:3]))
				ceta.fp.int=t(as.matrix(ceta.fp.int))
				}
				}


m.tp.int=NULL
m.fp.int=NULL
if(is.null(ceta.tp.int)){l.tp.int=0}else{l.tp.int=nrow(ceta.tp.int)}# number of TP
if(is.null(ceta.fp.int)){l.fp.int=0}else{l.fp.int=nrow(ceta.fp.int)}
#l.fp.int=nrow(ceta.fp.int)# number of FP
if(l.tp.int>0){
fr.tp.int=as.numeric(ceta.tp.int[,8])#frequencies
m.tp.int=ceta.tp.int[,1:7]
p.tp.int=as.numeric(ceta.tp.int[,9])
}
if(l.fp.int>0){
p.fp.int=as.numeric(ceta.fp.int[,9])
fr.fp.int=as.numeric(ceta.fp.int[,8])
m.fp.int=ceta.fp.int[,1:7]
}
p.fp<-NULL
p.fp<-c(p.fp.add,p.fp.int)
fr.fp<-NULL
fr.fp<-c(fr.fp.add,fr.fp.int)
l.fp<-c(l.fp.add,l.fp.int)


##*******
# TP's
##*****

l.tp<-NULL
l.tp<-c(l.tp.add,l.tp.int)

tp<-NULL
tp<-c(tp.add,tp.int)

fr.tp<-NULL
fr.tp<-c(fr.tp.add,fr.tp.int)

p.tp<-NULL
p.tp<-c(p.tp.add,p.tp.int)

m.tp.add<-NULL
if(!is.null(tp.add)){
tp.add<-as.numeric(tp.add)
m.tp.add=tp.add
}


FDR=NULL
if(length(which((l.fp+l.tp)==0))>0){ 
FDR[which((l.fp+l.tp)==0)]=0
FDR[-which((l.fp+l.tp)==0)]=l.fp[-which((l.fp+l.tp)==0)]/(l.fp[-which((l.fp+l.tp)==0)]+l.tp[-which((l.fp+l.tp)==0)])
}else{FDR=l.fp/(l.fp+l.tp)}
FDR

if(!is.null(m.tp.int)) if(!is.matrix(m.tp.int)) m.tp.int=t(as.matrix(m.tp.int))
if(!is.null(m.fp.int)) if(!is.matrix(m.fp.int)) m.fp.int=t(as.matrix(m.fp.int))


TP.add=NULL
if(!is.null(tp.add)) TP.add=cbind(tp.add,fr.tp.add,p.tp.add)
TP.int=NULL
if(!is.null(m.tp.int)) TP.int=cbind(m.tp.int,fr.tp.int,p.tp.int)
FP.add=NULL
if(!is.null(fp.add)) FP.add=cbind(fp.add,fr.fp.add,p.fp.add)
FP.int=NULL
if(!is.null(m.fp.int)) FP.int=cbind(m.fp.int,fr.fp.int,p.fp.int)
colnames(TP.add)=colnames(TP.int)=colnames(FP.add)=colnames(FP.int)=NULL

# TP.add=cbind(tp.add,fr.tp.add,p.tp.add)
# TP.int=cbind(m.tp.int,fr.tp.int,p.tp.int)
# FP.add=cbind(fp.add,fr.fp.add,p.fp.add)
# FP.int=cbind(m.fp.int,fr.fp.int,p.fp.int)
# colnames(TP.add)=colnames(TP.int)=colnames(FP.add)=colnames(FP.int)=NULL

s.a.tp=NULL
s.i.tp=NULL
s.a.fp=NULL
s.i.fp=NULL
if(is.matrix(TP.add)){s.a.tp=sum(TP.add[,2])}else{
												  if(length(TP.add)>0){s.a.tp=sum(TP.add[2])
																		}else{s.a.tp=0}
												    }
if(is.matrix(TP.int)){s.i.tp=sum(TP.int[,8])}else{
													if(length(TP.int)>0){s.i.tp=sum(TP.int[8])
																			}else{s.i.tp=0}
													}
if(is.matrix(FP.add)){s.a.fp=sum(FP.add[,2])}else{
													if(length(FP.add)>0){s.a.fp=sum(FP.add[2])
																			}else{s.a.fp=0}
													}
if(is.matrix(FP.int)){s.i.fp=sum(FP.int[,8])}else{
													if(length(FP.int)>0){s.i.fp=sum(FP.int[8])
																			}else{s.i.fp=0}
													}

l.tp=c(s.a.tp,s.i.tp)
l.fp=c(s.a.fp,s.i.fp)

ltp=l.tp
lfp=l.fp

FDR=NULL
if(length(which((l.fp+l.tp)==0))>0){ 
FDR[which((l.fp+l.tp)==0)]=0
FDR[-which((l.fp+l.tp)==0)]=l.fp[-which((l.fp+l.tp)==0)]/(l.fp[-which((l.fp+l.tp)==0)]+l.tp[-which((l.fp+l.tp)==0)])
}else{FDR=l.fp/(l.fp+l.tp)}
FDR

FDR.total=0
if(length(which(sum(lfp)+sum(ltp)==0))>0){ 
FDR.total[which((sum(lfp)+sum(ltp))==0)]=0
FDR.total[-which((sum(lfp)+sum(ltp))==0)]=sum(lfp)[-which((sum(lfp)+sum(ltp))==0)]/(sum(lfp)[-which((sum(lfp)+sum(ltp))==0)]+sum(ltp)[-which((sum(lfp)+sum(ltp))==0)])
}else{FDR.total=sum(lfp)/(sum(lfp)+sum(ltp))}
FDR.total

# ## FP depending on the threshold
# FP.int.t=NULL
# FP.add.t=NULL
 # FP.int.t=FP.int[which(FP.int[,8]/n.sym>=thres),]
 # FP.add.t=FP.add[which(FP.add[,2]/n.sym>=thres),]
# TP.int.t=TP.int
# TP.add.t=TP.add
# colnames(TP.add)=colnames(TP.int)=colnames(FP.add)=colnames(FP.int)=NULL

# s.a.tp=NULL
# s.i.tp=NULL
# s.a.fp=NULL
# s.i.fp=NULL
# if(is.matrix(TP.add.t)){s.a.tp=sum(TP.add.t[,2])}else{
												  # if(length(TP.add.t)>0){s.a.tp=sum(TP.add.t[2])
																		# }else{s.a.tp=0}
												    # }
# if(is.matrix(TP.int.t)){s.i.tp=sum(TP.int.t[,8])}else{
													# if(length(TP.int.t)>0){s.i.tp=sum(TP.int.t[8])
																			# }else{s.i.tp=0}
													# }
# if(is.matrix(FP.add.t)){s.a.fp=sum(FP.add.t[,2])}else{
													# if(length(FP.add.t)>0){s.a.fp=sum(FP.add.t[2])
																			# }else{s.a.fp=0}
													# }
# if(is.matrix(FP.int.t)){s.i.fp=sum(FP.int.t[,8])}else{
													# if(length(FP.int.t)>0){s.i.fp=sum(FP.int.t[8])
																			# }else{s.i.fp=0}
													# }

# l.tp=c(s.a.tp,s.i.tp)
# l.fp=c(s.a.fp,s.i.fp)


# FDR.t=NULL
# if(length(which((l.fp+l.tp)==0))>0){ 
# FDR.t[which((l.fp+l.tp)==0)]=0
# FDR.t[-which((l.fp+l.tp)==0)]=l.fp[-which((l.fp+l.tp)==0)]/(l.fp[-which((l.fp+l.tp)==0)]+l.tp[-which((l.fp+l.tp)==0)])
# }else{FDR.t=l.fp/(l.fp+l.tp)}
# FDR.t

return(list(TP.add=TP.add,TP.int=TP.int,FP.add=FP.add,FP.int=FP.int,FDR.a=FDR[1],FDR.i=FDR[2],FDR=FDR.total))
}

##
tp.fp.ILA=function(true.int,tab.add.p,tab.int.p,thres=0.15)
### true and false positives for ILA approach
{
##***********************************************************
l.fp.add<-0
fp.add<-NULL
fr.fp.add<-NULL
p.fp.add<-NULL
l.fp.int<-0
fp.int<-NULL
fr.fp.int<-NULL
p.fp.int<-NULL

l.tp.add<-0
tp.add<-NULL
fr.tp.add<-NULL
p.tp.add<-NULL
l.tp.int<-0
tp.int<-NULL
fr.tp.int<-NULL
p.tp.int<-NULL

if(!is.null(true.int)){
 write.table(true.int, file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  true.int1=read.table(file = "M.int.txt")
  fn <- "M.int.txt"
  if (file.exists(fn)) file.remove(fn)
  }

 
where.int=which(unlist(lapply(apply(true.int1,1,tree.operators),length))>0)#which of true int are really interactions

true.leaves=sort(as.numeric(unlist(apply(true.int1,1,variables.from.model))))
Adds=true.leaves
#index of true additive effect in tab.add.p (row)
# ind.Adds=rep(0,length(Adds))
ind.Adds=NULL
for(jj in 1:length(Adds)){
				if(length(which(tab.add.p[,1]==Adds[jj]))>0)
				ind.Adds<-c(ind.Adds,which(tab.add.p[,1]==Adds[jj]))
				}
				
fr.tp.add=tab.add.p[ind.Adds,2]
p.tp.add=tab.add.p[ind.Adds,3]
l.tp.add=length(which(ind.Adds>0))
tp.add=tab.add.p[ind.Adds,1]
if(!is.null(ind.Adds)){
fr.fp.add=tab.add.p[-ind.Adds,2]
p.fp.add=tab.add.p[-ind.Adds,3]
l.fp.add=nrow(tab.add.p)-l.tp.add
l.fp.add
fp.add=tab.add.p[-ind.Adds,1]
}

## 
int.P=NULL
if(!is.null(tab.int.p))
if(nrow(tab.int.p)>0){
write.table(tab.int.p[,1], "int-p.txt",append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
intP=read.table("int-p.txt")
fn <- "int-p.txt"
if (file.exists(fn)) file.remove(fn)
int.P=t(apply(intP,1,tree.merge))
if(length(apply(int.P,1,intpaste))>1){
tab.int.P=cbind(apply(int.P,1,intpaste),tab.int.p[,c(2,3)])}else{tab.int.P=t(as.matrix(c(apply(int.P,1,intpaste),tab.int.p[,c(2,3)])))}
tab.int.p<-tab.int.P
}

Ints=true.int1[where.int,]
# all true leaves

##################################################
leaves.list=apply(true.int1[where.int,],1,variables.from.model)
if(!is.list(leaves.list)){
 leaves.list=split(leaves.list, rep(1:ncol(leaves.list), each = nrow(leaves.list)))
}
vars.int=NULL
if(!is.null(int.P)){
vars.int=apply(int.P,1,variables.from.model)
if(!is.list(vars.int)){
 vars.int=split(vars.int, rep(1:ncol(vars.int), each = nrow(vars.int)))
}
}
# lista indeksow interakcji ktore zawieraja true ints jako podzbiory
ind.tp.ILA=lapply(leaves.list,subv.in.tree,List=vars.int)
ind.tp.ILA

for(kk in 1:length(ind.tp.ILA)){
if(length(leaves.list[[kk]])>1){
## subpairs 
if(!is.null(subvectors(leaves.list[[kk]])$subpairs)){
ind.subs=NULL
ind.subs=apply(subvectors(leaves.list[[kk]])$subpairs,1,subv.in.tree,List=vars.int)## all which contain subpairs
ind.sub.p=NULL
ind.sub.p=which(unlist(lapply(vars.int[unlist(ind.subs)],length))==2)
ind.p=0
if(length(ind.sub.p)>0) ind.p=unlist(ind.subs)[ind.sub.p]
ind.tp.ILA[[kk]]=c(ind.tp.ILA[[kk]],ind.p)
}
## subtriples
if(!is.null(subvectors(leaves.list[[kk]])$subtriples)){
ind.subs=apply(subvectors(leaves.list[[kk]])$subtriples,1,subv.in.tree,List=vars.int)## all which contain subpairs
ind.sub.p=which(unlist(lapply(vars.int[unlist(ind.subs)],length))==3)
ind.p=0
if(length(ind.sub.p)>0) ind.p=unlist(ind.subs)[ind.sub.p]
ind.tp.ILA[[kk]]=c(ind.tp.ILA[[kk]],ind.p)
}
ind.tp.ILA[[kk]]=unique(ind.tp.ILA[[kk]])
}
ind.tp.ILA[[kk]]=ind.tp.ILA[[kk]][which(ind.tp.ILA[[kk]]>0)]
}


# ind.tp.ILA=which(rowSums(t(apply(int.P,1,contain.add,vec=true.leaves)))>0)## where the ILA TPs are



ila.tp.i=vector("list",length(leaves.list))
ila.tp.i.mat=NULL
for(kk in 1:length(ind.tp.ILA)){
ila.tp.int=NULL
if(length(ind.tp.ILA[[kk]])>1){
				ila.tp.int=cbind(int.P[ind.tp.ILA[[kk]],],matrix(as.numeric(tab.int.p[ind.tp.ILA[[kk]],2:3]),ncol=2))}
if(length(ind.tp.ILA[[kk]])==1){
				ila.tp.int=t(as.matrix(c(round(int.P[ind.tp.ILA[[kk]],],digits=0),matrix(as.numeric(tab.int.p[ind.tp.ILA[[kk]],2:3]),ncol=2))))}
ila.tp.i[[kk]]=ila.tp.int
if(length(ila.tp.int)==9) ila.tp.i.mat=rbind(ila.tp.i.mat,ila.tp.int)
}
# ila.tp.i.mat
ila.fp.int=NULL
ind.tp.ILA=as.vector(unlist(ind.tp.ILA))
if(length(ind.tp.ILA)>0){INT=int.P[-ind.tp.ILA,]}else{INT=int.P}
if(!is.null(INT)){
if(is.matrix(INT)){
if(length(ind.tp.ILA)>0){
ila.fp.int=cbind(INT,matrix(as.numeric(tab.int.p[-ind.tp.ILA,2:3]),ncol=2))
}
if(length(ind.tp.ILA)==0){
ila.fp.int=cbind(INT,matrix(as.numeric(tab.int.p[,2:3]),ncol=2))
}
}
if(!is.matrix(INT)){
INT=t(as.matrix(INT))
if(length(ind.tp.ILA)>0){
ila.fp.int=cbind(INT,matrix(as.numeric(tab.int.p[-ind.tp.ILA,2:3]),ncol=2))
}
if(length(ind.tp.ILA)==0){
ila.fp.int=cbind(INT,matrix(as.numeric(tab.int.p[,2:3]),ncol=2))
}
}
}
# l.tp.int=nrow(ila.tp.i.mat)#sum(unlist(lapply(ila.tp.i,nrow)))#nrow(ila.tp.int)# number of TP
# l.fp.int=nrow(ila.fp.int)# number of FP

# fr.tp.int=as.numeric(ila.tp.i.mat[,8])#frequencies
# p.tp.int=as.numeric(ila.tp.i.mat[,9])
# p.fp.int=as.numeric(ila.fp.int[,9])
# fr.fp.int=as.numeric(ila.fp.int[,8])

# m.tp.int=ila.tp.i.mat[,1:7]
# m.fp.int=ila.fp.int[,1:7]

# p.fp<-NULL
# p.fp<-c(p.fp.add,p.fp.int)
# fr.fp<-NULL
# fr.fp<-c(fr.fp.add,fr.fp.int)
# l.fp<-c(l.fp.add,l.fp.int)



m.tp.int=NULL
m.fp.int=NULL
if(is.null(ila.tp.int)){l.tp.int=0}else{l.tp.int=nrow(ila.tp.int)}# number of TP
if(is.null(ila.fp.int)){l.fp.int=0}else{l.fp.int=nrow(ila.fp.int)}
#l.fp.int=nrow(ila.fp.int)# number of FP
if(l.tp.int>0){
fr.tp.int=as.numeric(ila.tp.int[,8])#frequencies
m.tp.int=ila.tp.int[,1:7]
p.tp.int=as.numeric(ila.tp.int[,9])
}
if(l.fp.int>0){
p.fp.int=as.numeric(ila.fp.int[,9])
fr.fp.int=as.numeric(ila.fp.int[,8])
m.fp.int=ila.fp.int[,1:7]
}
p.fp<-NULL
p.fp<-c(p.fp.add,p.fp.int)
fr.fp<-NULL
fr.fp<-c(fr.fp.add,fr.fp.int)
l.fp<-c(l.fp.add,l.fp.int)



##*******
# TP's
# ## saving the level-related tp's to file
##*****

l.tp<-NULL
l.tp<-c(l.tp.add,l.tp.int)



fr.tp<-NULL
fr.tp<-c(fr.tp.add,fr.tp.int)

p.tp<-NULL
p.tp<-c(p.tp.add,p.tp.int)

m.tp.add<-NULL
if(!is.null(tp.add)){
tp.add<-as.numeric(tp.add)
m.tp.add=tp.add
}


# FDR=NULL
# if(length(which((l.fp+l.tp)==0))>0){ 
# FDR[which((l.fp+l.tp)==0)]=0
# FDR[-which((l.fp+l.tp)==0)]=l.fp[-which((l.fp+l.tp)==0)]/(l.fp[-which((l.fp+l.tp)==0)]+l.tp[-which((l.fp+l.tp)==0)])
# }else{FDR=l.fp/(l.fp+l.tp)}
# FDR
if(!is.null(m.tp.int)){
					if(!is.matrix(m.tp.int)) m.tp.int=t(as.matrix(m.tp.int))}
if(!is.null(m.fp.int)){
					if(!is.matrix(m.fp.int)) m.fp.int=t(as.matrix(m.fp.int))
					}

TP.add=cbind(tp.add,fr.tp.add,p.tp.add)
TP.int=cbind(m.tp.int,fr.tp.int,p.tp.int)
FP.add=cbind(fp.add,fr.fp.add,p.fp.add)
FP.int=cbind(m.fp.int,fr.fp.int,p.fp.int)
colnames(TP.add)=colnames(TP.int)=colnames(FP.add)=colnames(FP.int)=NULL

s.a.tp=NULL
s.i.tp=NULL
s.a.fp=NULL
s.i.fp=NULL
if(is.matrix(TP.add)){if(dim(TP.add)[1]>0){s.a.tp=sum(TP.add[,2])}else{s.a.tp=0}}else{
												  if(length(TP.add)>0){s.a.tp=sum(TP.add[2])
																		}else{s.a.tp=0}
												    }
if(is.matrix(TP.int)){if(dim(TP.int)[1]>0){s.i.tp=sum(TP.int[,8])}else{s.i.tp=0}}else{
													if(length(TP.int)>0){s.i.tp=sum(TP.int[8])
																			}else{s.i.tp=0}
													}
if(is.matrix(FP.add)){if(dim(FP.add)[1]>0){s.a.fp=sum(FP.add[,2])}else{s.a.fp=0}}else{
													if(length(FP.add)>0){s.a.fp=sum(FP.add[2])
																			}else{s.a.fp=0}
													}
if(is.matrix(FP.int)){if(dim(FP.int)[1]>0){s.i.fp=sum(FP.int[,8])}else{s.i.fp=0}}else{
													if(length(FP.int)>0){s.i.fp=sum(FP.int[8])
																			}else{s.i.fp=0}
													}

l.tp=c(s.a.tp,s.i.tp)
l.fp=c(s.a.fp,s.i.fp)

ltp=l.tp
lfp=l.fp

FDR=rep(0,2)
if(length(which((l.fp+l.tp)==0))>0){ 
FDR[which((l.fp+l.tp)==0)]=0
FDR[-which((l.fp+l.tp)==0)]=l.fp[-which((l.fp+l.tp)==0)]/(l.fp[-which((l.fp+l.tp)==0)]+l.tp[-which((l.fp+l.tp)==0)])
}else{FDR=l.fp/(l.fp+l.tp)}
FDR


FDR.total=0
if(length(which(sum(lfp)+sum(ltp)==0))>0){ 
FDR.total[which((sum(lfp)+sum(ltp))==0)]=0
FDR.total[-which((sum(lfp)+sum(ltp))==0)]=sum(lfp)[-which((sum(lfp)+sum(ltp))==0)]/(sum(lfp)[-which((sum(lfp)+sum(ltp))==0)]+sum(ltp)[-which((sum(lfp)+sum(ltp))==0)])
}else{FDR.total=sum(lfp)/(sum(lfp)+sum(ltp))}
FDR.total

# ## FP depending on the threshold
# FP.int.t=NULL
# FP.add.t=NULL
 # FP.int.t=FP.int[which(FP.int[,8]/n.sym>=thres),]
 # FP.add.t=FP.add[which(FP.add[,2]/n.sym>=thres),]
# TP.int.t=TP.int
# TP.add.t=TP.add
# colnames(TP.add)=colnames(TP.int)=colnames(FP.add)=colnames(FP.int)=NULL

# s.a.tp=NULL
# s.i.tp=NULL
# s.a.fp=NULL
# s.i.fp=NULL
# if(is.matrix(TP.add.t)){s.a.tp=sum(TP.add.t[,2])}else{
												  # if(length(TP.add.t)>0){s.a.tp=sum(TP.add.t[2])
																		# }else{s.a.tp=0}
												    # }
# if(is.matrix(TP.int.t)){s.i.tp=sum(TP.int.t[,8])}else{
													# if(length(TP.int.t)>0){s.i.tp=sum(TP.int.t[8])
																			# }else{s.i.tp=0}
													# }
# if(is.matrix(FP.add.t)){s.a.fp=sum(FP.add.t[,2])}else{
													# if(length(FP.add.t)>0){s.a.fp=sum(FP.add.t[2])
																			# }else{s.a.fp=0}
													# }
# if(is.matrix(FP.int.t)){s.i.fp=sum(FP.int.t[,8])}else{
													# if(length(FP.int.t)>0){s.i.fp=sum(FP.int.t[8])
																			# }else{s.i.fp=0}
													# }

# l.tp=c(s.a.tp,s.i.tp)
# l.fp=c(s.a.fp,s.i.fp)


# FDR.t=rep(0,2)
# if(length(which((l.fp+l.tp)==0))>0){ 
# FDR.t[which((l.fp+l.tp)==0)]=0
# FDR.t[-which((l.fp+l.tp)==0)]=l.fp[-which((l.fp+l.tp)==0)]/(l.fp[-which((l.fp+l.tp)==0)]+l.tp[-which((l.fp+l.tp)==0)])
# }else{FDR.t=l.fp/(l.fp+l.tp)}
# FDR.t

return(list(TP.add=TP.add,TP.int=TP.int,FP.add=FP.add,FP.int=FP.int,FDR.a=FDR[1],FDR.i=FDR[2],FDR=FDR.total))
#,TP.add.t=TP.add.t,TP.int.t=TP.int.t,FP.add.t=FP.add.t,FP.int.t=FP.int.t,FDR.a.t=FDR.t[1],FDR.i.t=FDR.t[2],FDR.t=sum(l.fp)/(sum(l.fp)+sum(l.tp)),ila.tp.i=ila.tp.i))
}

##

tp.fp.ILA.plot=function(true.int,tab.add.p,tab.int.p,max.thres=max.thres,OutFile=OutFile)
{
FDR=NULL
FDR.a.t=NULL
FDR.i.t=NULL
for(hh in 0:10){
res=tp.fp.ILA(true.int,tab.add.p,tab.int.p,thres=max.thres*hh*0.1)
FDR.a.t=c(FDR.a.t,res$FDR.a.t)
FDR.i.t=c(FDR.i.t,res$FDR.i.t)
}
jpeg(file=file.path(OutFile,"_ILA_FDR-t.jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
plot(c(0:10)*0.1*max.thres,FDR.a.t,ylim=c(0,1),xlab="Threshold",ylab="FDR", main="FDR for ILA",type="l",col="red",lty=2,lwd=2)
points(c(0:10)*0.1*max.thres,FDR.i.t,type="l",col="blue",lty=1,lwd=2)
legend(x=0.3,y=0.2,c("FDR add", "FDR int"),lty=c(2,1),col=c("red","blue"),inset = .01,cex=0.6,bty="n", title.adj=.55)
dev.off()
return(list(FDR.a.t=FDR.a.t,FDR.i.t=FDR.i.t))
}



tp.fp=function(true.int,tab.add.p,tab.int.p)
{
##***********************************************************
l.fp.add<-0
fp.add<-NULL
fr.fp.add<-NULL
p.fp.add<-NULL
l.fp.int<-0
fp.int<-NULL
fr.fp.int<-NULL
p.fp.int<-NULL

if(!is.null(true.int)){
 write.table(true.int, file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  true.int1=read.table(file = "M.int.txt")
  fn <- "M.int.txt"
  if (file.exists(fn)) file.remove(fn)
  }

  
 # tab.add.p[which(tab.add.p[,1]==true.int1[-which(unlist(lapply(apply(true.int1,1,tree.operators),length))>0),1][1]),]
 
where.int=which(unlist(lapply(apply(true.int1,1,tree.operators),length))>0)#which of true int are really interactions
Adds=true.int1[-where.int,1]

#index of true additive effect in tab.add.p (row)
ind.Adds=rep(0,length(Adds))
for(jj in 1:length(Adds)) ind.Adds[jj]=which(tab.add.p[,1]==Adds[jj])
fr.tp.add=tab.add.p[ind.Adds,2]
p.tp.add=tab.add.p[ind.Adds,3]
l.tp.add=length(which(ind.Adds>0))

fr.fp.add=tab.add.p[-ind.Adds,2]
p.fp.add=tab.add.p[-ind.Adds,3]
l.fp.add=nrow(tab.add.p)-l.tp.add
l.fp.add

if(!is.null(true.int)){
if(!is.null(tab.add.p)){
if(is.matrix(tab.add.p)){
d<-matrix(nrow=length(true.int), ncol=length(tab.add.p[,1]))
if(nrow(tab.add.p)>0){
for(s in 1:length(true.int)) for(k in 1:length(tab.add.p[,1])) d[s,k]<-c(tab.add.p[,1][k]!=true.int[s])
for(k in 1:length(tab.add.p[,1]))
if(sum(d[,k])==length(true.int) ) { l.fp.add<-l.fp.add+1
					fp.add<-c(fp.add,tab.add.p[,1][k])
					fr.fp.add<-c(fr.fp.add,tab.add.p[,2][k])
					p.fp.add<-c(p.fp.add,tab.add.p[,3][k])
									}
}
}


if(!is.matrix(tab.add.p)){
d<-matrix(nrow=length(true.int), ncol=length(tab.add.p[1]))
if(length(tab.add.p)>0){
for(s in 1:length(true.int))  d[s,1]<-c(tab.add.p[1]!=true.int[s])
if(sum(d[,1])==length(true.int) ) { l.fp.add<-l.fp.add+1
					fp.add<-c(fp.add,tab.add.p[1])
					fr.fp.add<-c(fr.fp.add,tab.add.p[2])
					p.fp.add<-c(p.fp.add,tab.add.p[3])
									}
}
}
}

if(!is.null(tab.int.p))
if(nrow(tab.int.p)>0){
write.table(tab.int.p[,1], "int-p.txt",append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
intP=read.table("int-p.txt")
fn <- "int-p.txt"
if (file.exists(fn)) file.remove(fn)
int.P=t(apply(intP,1,tree.merge))
if(length(apply(int.P,1,intpaste))>1){
tab.int.P=cbind(apply(int.P,1,intpaste),tab.int.p[,c(2,3)])}else{tab.int.P=t(as.matrix(c(apply(int.P,1,intpaste),tab.int.p[,c(2,3)])))}
tab.int.p<-tab.int.P
}

Ints=true.int1[where.int,]

#index of true interaction effect in tab.int.p (row)
ind.Ints=rep(0,nrow(Ints))
for(jj in 1:nrow(Ints)) ind.Ints[jj]=match(true.int[jj],tab.int.p[,1],nomatch=FALSE)
fr.tp.int=as.numeric(tab.int.p[ind.Ints,2])
p.tp.int=as.numeric(tab.int.p[ind.Ints,3])
l.tp.int=length(which(ind.Ints>0))

fr.fp.int=as.numeric(tab.int.p[-ind.Ints,2])
p.fp.int=as.numeric(tab.int.p[-ind.Ints,3])
l.fp.int=nrow(tab.int.p)-l.tp.int
l.fp.int





if(!is.null(tab.int.p)){
if(!is.null(nrow(tab.int.p))){
d<-matrix(nrow=length(true.int), ncol=length(tab.int.p[,1]))
if(nrow(tab.int.p)>0){
for(s in 1:length(true.int)) for(k in 1:length(tab.int.p[,1])) d[s,k]<-c(tab.int.p[,1][k]!=true.int[s])
for(k in 1:length(tab.int.p[,1]))
if(sum(d[,k])==length(true.int) ) { l.fp.int<-l.fp.int+1
					fp.int<-c(fp.int,tab.int.p[,1][k])
					fr.fp.int<-c(fr.fp.int,tab.int.p[,2][k])
					p.fp.int<-c(p.fp.int,tab.int.p[,3][k])
									}
}
}
}

if(is.null(nrow(tab.int.p))){
d<-matrix(nrow=length(true.int), ncol=length(tab.int.p[1]))
if(length(tab.int.p)>0){
for(s in 1:length(true.int))  d[s,1]<-c(tab.int.p[1]!=true.int[s])
if(sum(d[,1])==length(true.int) ) { l.fp.int<-l.fp.int+1
					fp.int<-c(fp.int,tab.int.p[1])
					fr.fp.int<-c(fr.fp.int,tab.int.p[2])
					p.fp.int<-c(p.fp.int,tab.int.p[3])
									}
}
}

l.fp<-NULL
l.fp<-c(l.fp.add,l.fp.int)
fp.add<-as.numeric(fp.add)
}else{if(is.matrix(tab.add.p)){
		l.fp.add<-nrow(tab.add.p)
		fp.add<-tab.add.p[,1]
		fr.fp.add<-tab.add.p[,2]
		p.fp.add<-tab.add.p[,3]
		}else{
		l.fp.add<-1
		fp.add<-tab.add.p[1]
		fr.fp.add<-tab.add.p[2]
		p.fp.add<-tab.add.p[3]
		}
	if(is.matrix(tab.int.p)){
		l.fp.int<-nrow(tab.int.p)
		fp.int<-tab.int.p[,1]
		fr.fp.int<-tab.int.p[,2]
		p.fp.int<-as.numeric(tab.int.p[,3])
		}else{
		l.fp.int<-1
		fp.int<-tab.int.p[1]
		fr.fp.int<-tab.int.p[2]
		p.fp.int<-as.numeric(tab.int.p[3])
		}
		}

p.fp<-NULL
p.fp<-c(p.fp.add,p.fp.int)
fr.fp<-NULL
fr.fp<-c(fr.fp.add,fr.fp.int)
m.fp.int<-NULL
l.fp<-c(l.fp.add,l.fp.int)
if(length(fp.int)>0){
write.table(fp.int, file = paste("fpint.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
fp.int.v<-read.table(file = paste("fpint.txt",sep=""))
fn <- "fpint.txt"
if (file.exists(fn)) file.remove(fn)
m.fp.int<-fp.int.v
}

##*******
# ## POPRAWNE "KLASYFIKACJE" 
##*****
tp.add<-NULL
fr.tp.add<-NULL
p.tp.add<-NULL
l.tp.add<-0
tp.int<-NULL
fr.tp.int<-NULL
p.tp.int<-NULL
l.tp.int<-0

if(!is.null(true.int))
{
if(is.matrix(tab.add.p)){
if(nrow(tab.add.p)>0){
d<-matrix(nrow=length(true.int), ncol=length(tab.add.p[,1]))
if(nrow(tab.add.p)>0){
for(s in 1:length(true.int)) for(k in 1:length(tab.add.p[,1])) d[s,k]<-c(tab.add.p[,1][k]==true.int[s])
for(s in 1:length(true.int)) for(k in 1:length(tab.add.p[,1]))
if(d[s,k]=='TRUE'){l.tp.add<-l.tp.add+1
				tp.add<-c(tp.add,tab.add.p[,1][[k]])
				fr.tp.add<-c(fr.tp.add,tab.add.p[,2][[k]])
				p.tp.add<-c(p.tp.add,tab.add.p[,3][k])
				}
}
}
}

if(!is.matrix(tab.add.p)){
if(length(tab.add.p)>0){
d<-matrix(rep(0,length(true.int)*length(tab.add.p[1])),ncol=length(tab.add.p[1]))
for(s in 1:length(true.int)) d[s,1]<-c(tab.add.p[1]==true.int[s])
for(s in 1:length(true.int)) 
if(d[s,1]=='TRUE'){l.tp.add<-l.tp.add+1
				tp.add<-c(tp.add,tab.add.p[1])
				fr.tp.add<-c(fr.tp.add,tab.add.p[2])
				p.tp.add<-c(p.tp.add,tab.add.p[3])
				}
}
}


if(is.matrix(tab.int.p)){
if(nrow(tab.int.p)>0){
d<-matrix(nrow=length(true.int), ncol=length(tab.int.p[,1]))
for(s in 1:length(true.int)) for(k in 1:length(tab.int.p[,1])) d[s,k]<-c(tab.int.p[,1][k]==true.int[s])
for(s in 1:length(true.int)) for(k in 1:length(tab.int.p[,1]))
if(d[s,k]=='TRUE'){l.tp.int<-l.tp.int+1
				tp.int<-c(tp.int,tab.int.p[,1][[k]])
				fr.tp.int<-c(fr.tp.int,tab.int.p[,2][[k]])
				p.tp.int<-c(p.tp.int,tab.int.p[,3][k])
				}
}
}

if(!is.matrix(tab.int.p)){
if(length(tab.int.p)>0){
d<-matrix(nrow=length(true.int), ncol=length(tab.int.p[1]))
for(s in 1:length(true.int)) d[s,1]<-c(tab.int.p[1]==true.int[s])
for(s in 1:length(true.int)) 
if(d[s,1]=='TRUE'){l.tp.int<-l.tp.int+1
				tp.int<-c(tp.int,tab.int.p[1])
				fr.tp.int<-c(fr.tp.int,tab.int.p[2])
				p.tp.int<-c(p.tp.int,tab.int.p[3])
				}
}
}




}


l.tp<-NULL
l.tp<-c(l.tp.add,l.tp.int)

tp<-NULL
tp<-c(tp.add,tp.int)

fr.tp<-NULL
fr.tp<-c(fr.tp.add,fr.tp.int)

p.tp<-NULL
p.tp<-c(p.tp.add,p.tp.int)

m.tp.add<-NULL
if(!is.null(tp.add)){
tp.add<-as.numeric(tp.add)
m.tp.add=tp.add
}

m.tp.int<-NULL
if(!is.null(tp.int)){
  write.table(tp.int, file= "tpint.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  tp.int.v<-read.table(file = "tpint.txt")
  fn <- "tpint.txt"
  if (file.exists(fn)) file.remove(fn)
  m.tp.int<-tp.int.v
}
FDR=NULL
if(length(which((l.fp+l.tp)==0))>0){ 
FDR[which((l.fp+l.tp)==0)]=0
FDR[-which((l.fp+l.tp)==0)]=l.fp[-which((l.fp+l.tp)==0)]/(l.fp[-which((l.fp+l.tp)==0)]+l.tp[-which((l.fp+l.tp)==0)])
}else{FDR=l.fp/(l.fp+l.tp)}
TP.add=cbind(tp.add,fr.tp.add,p.tp.add)
TP.int=cbind(m.tp.int,fr.tp.int,p.tp.int)
FP.add=cbind(fp.add,fr.fp.add,p.fp.add)
FP.int=cbind(m.fp.int,fr.fp.int,p.fp.int)
return(list(TP.add=TP.add,TP.int=TP.int,FP.add=FP.add,FP.int=FP.int,FDR.a=FDR[1],FDR.i=FDR[2],FDR=sum(l.fp)/(sum(l.fp)+sum(l.tp))))
}

#######################################################################
find.max.fr=function(vec,nitr)
## function to determine max frequencies for corresponding expressions
# returns a vector of frequencies larger than 0.5 and corresponding indices of interactions
{
post.fr=NULL
Ind=which((vec/nitr)>0)#1/nitr)
post.fr=vec[Ind]/nitr
return(cbind(post.fr=post.fr,Maxi=Ind))
}

##########################################################################
model.PIMCLR=function(res,nitr,outFile="mypost_results.txt",filename="slogiclisting.tmp",nlvs=7,ntrs=5)
# function to get all prime implicants and their posterior freqencies
# this corresponds to the original MCLR method.
{
 # res=fitmodel(
				 # nsnp, nrows,
				 # nitr,
				 # ntrs,	nlvs,
				 # p=0.5,a=0.7,n.tree.mod=12,tree.lvs=c(4,3,3,2,2,2,1,1,1,1,1,1),eff.sizes= Bvec[length(Bvec):1],slope=1,tree.ind=Tind,op=c("&","&","&","&","&","&","||","&","&","&","&"))
outFile3=paste(outFile,"_MCLR_CETA",sep="")		
add=NULL
PP=find.max.fr(res$single,nitr)#single
if(!is.null(PP)) add=PP[,ncol(PP):1]
TabA=NULL
TabAdd=NULL
if(length(which(add[,2]>=0.5))>0) TabA=add[which(add[,2]>=0.5),]
if(!is.null(TabA)) TabAdd=cbind(TabA[,1],rep(1,nrow(TabA)),TabA[,2])

DD=apply(res$double,2,find.max.fr,nitr=nitr)#doubles
DDn=lapply(DD,is.null)
LV=which(unlist(DDn)==FALSE)## which variables interact 
int2=NULL
if(length(LV)>0){
for( i in 1:length(LV)){
pom=cbind(rep(LV[i],nrow(DD[[LV[i]]])),as.numeric(DD[[LV[i]]][,2]),as.numeric(DD[[LV[i]]][,1]))
int2=rbind(int2,as.matrix(pom))
}
}
TabInt=NULL
if(!is.null(int2)){
TabInt=t(mapply(make.tree,rep(1,nrow(int2)),lnull=4,n0v=1000,lmv=2))
TabInt[,2:3]<-int2[,1:2]
TabInt=cbind(TabInt,int2[,3])
Ind=which(as.numeric(TabInt[,8])>=0.5)
if(length(Ind)>1) TabInt=TabInt[which(as.numeric(TabInt[,8])>=0.5),]
if(length(Ind)==1){
				TabInt=TabInt[which(as.numeric(TabInt[,8])>=0.5),]
				TabInt=t(as.matrix(TabInt))
				}
if(length(Ind)==0) TabInt=NULL
}
TabInt

int3=NULL
for(i in 1:dim(res$triple)[3]){
Ti=res$triple[,,i]
TTi=apply(Ti,2,find.max.fr,nitr=nitr)#triples
TTn=lapply(TTi,is.null)
LV=which(unlist(TTn)==FALSE)## which variables interact 
if(length(LV)>0){
for(kk in 1:length(LV)){
pom=cbind(rep(i,nrow(TTi[[LV[kk]]])),rep(LV[kk],nrow(TTi[[LV[kk]]])),as.numeric(TTi[[LV[kk]]][,2]),as.numeric(TTi[[LV[kk]]][,1]))
int3=rbind(int3,as.matrix(pom))
			}
		}
}

TabInt3=NULL
if(!is.null(int3)){
TabInt3=t(mapply(make.tree,rep(2,nrow(int3)),lnull=3,n0v=1000,lmv=2))
TabInt3[,3:5]<-int3[,1:3]
TabInt3=cbind(TabInt3,int3[,4])
Ind=which(as.numeric(TabInt3[,8])>=0.5)
if(length(Ind)>1) TabInt3=TabInt3[which(as.numeric(TabInt3[,8])>=0.5),]
if(length(Ind)==1){
				TabInt3=TabInt3[which(as.numeric(TabInt3[,8])>=0.5),]
				TabInt3=t(as.matrix(TabInt3))
				}
if(length(Ind)==0) TabInt3=NULL
}
TabInt3
TabInt<-rbind(TabInt,TabInt3)


res.int2=NULL
if(!is.null(int2)){
				if(nrow(int2)>1){
				mpom=rbind(matrix(rep(1000,nrow(int2)),ncol=nrow(int2)),t(as.matrix(int2[,1:2])),matrix(rep(0,4*nrow(int2)),ncol=nrow(int2)))
				res.int2=cbind(apply(t(mpom),1,intpaste),int2[,3])
				}else{
				mpom=rbind(matrix(rep(1000,nrow(int2)),ncol=nrow(int2)),as.matrix(int2[,1:2]),matrix(rep(0,4*nrow(int2)),ncol=nrow(int2)))
				res.int2=cbind(apply(t(mpom),1,intpaste),int2[,3])
				}
				}
				
res.int3=NULL
if(!is.null(int3)){
				if(nrow(int3)>1){
				mpom=rbind(matrix(rep(1000,2*nrow(int3)),ncol=nrow(int3)),t(as.matrix(int3[,1:3])),matrix(rep(0,2*nrow(int3)),ncol=nrow(int3)))
				res.int3=cbind(apply(t(mpom),1,intpaste),int3[,4])
				}else{
				mpom=rbind(matrix(rep(1000,2*nrow(int3)),ncol=nrow(int3)),as.matrix(int3[,1:3]),matrix(rep(0,2*nrow(int3)),ncol=nrow(int3)))
				res.int3=cbind(apply(t(mpom),1,intpaste),int3[,4])
				}
				}


			
stats=read.table(filename)
scores=stats[,2]
stats=stats[,4:ncol(stats)]
stats.seq=data.frame(t(stats))
stats.seq=as.list(stats.seq)
T=lapply(stats.seq,vec.trees,nlvs=nlvs,ntrs=ntrs,indic=1)
nmlvs=apply(stats,1,vec.leaves)
nmtrs=apply(stats,1,vec.trees,nlvs=nlvs,ntrs=ntrs,indic=0)
nm=length(T)# number of models on the list 
#######		
# newscores
newscores=NULL
newscores=-scores#-0.5*nmtrs*log(nrows)
Enewscores=NULL
Enewscores=exp(newscores)
#############################
##############################
priors=rep(1,nm)
###########
##########################################
## now calculate posteriors for a specific trees
###########################################
outFile2=paste(outFile,"_ETA.txt",sep="")
Addind=lapply(T,Int.tree.ind,ntrs=ntrs,nlvs=nlvs,lev=0)
indT.a=which(as.numeric(unlist(lapply(Addind,length)))>0)# which models contain at least one add effect
PIa=vector("list",length(T))
if(length(indT.a)>0){
for(ww in 1:length(indT.a))
PIa[[indT.a[[ww]]]]=take.adds(T[[indT.a[[ww]]]],ind=as.numeric(unlist(Addind[indT.a[[ww]]])))
}
alladds=as.numeric(unlist(PIa))
#############
# list of all different additive effects 
unadds=unique(sort(abs(alladds)))
#########################################
# take an additive effect , take the prior of all models that contain it, take corresponding score and sum the prior times the new score
Ia<-NULL
if(length(unadds)>0){
Ia=vector("list",length(unadds))
for( i in 1:length(unadds)){
Ij=NULL# models which contain a particular add effect
for( j in 1:length(PIa)) 
if((!is.na(match(unadds[i],PIa[[j]])))||(!is.na(match(-unadds[i],PIa[[j]]))))# ith additive effect appears in jth model
Ij<-c(Ij,j)
Ia[[i]]=Ij# save indices of such models 
}
#########################################
}
## data frame saved to file 
addposts=NULL
res.add=NULL
if(length(Ia)>0){
for(i in 1:length(Ia))
addposts[i]=sum(priors[Ia[[i]]]*Enewscores[Ia[[i]]])/sum(priors*Enewscores)
foo=cbind(unadds,round(addposts,digits=6))
if(length(which(foo[,2]>0))>=0){
res.add=foo[order(foo[,2],foo[,1],decreasing=TRUE),]
if(is.vector(res.add)) res.add=t(as.matrix(res.add))
#res.add=res.add[which(res.add[,2]>=0),]
}
}
#########################################
###### interactions 
Iind=lapply(T,Int.tree.ind,ntrs=ntrs,nlvs=nlvs,lev=1)
PIl=vector("list",length(Iind))
for(ww in 1:length(Iind)) PIl[[ww]]=take.ints(mat=T[[ww]],ind=Iind[[ww]]) 
#PIl= mapply(take.ints,T,ind=Iind)# all trees with interactions
PIi=do.call(rbind, PIl)#unlisted- as a matrix with trees in rows
unints=NULL
intposts=NULL
foo=NULL
sig.int=NULL
res.int=NULL
sigint=NULL
if(length(PIi)>0)
{
PIm=t(apply(PIi,1,tree.merge))# posortuj zawartosc drzew
m.ind=model.ind(nmtrs)
if(nrow(PIm)>2){PI=count.rows(PIm)
				Iint=PI[,2]
				Iint=lapply(Iint,function(i) m.ind[i])
				unints=PI[,c(3:ncol(PI))]
				}else{
				PIm=matrix(PIm,ncol=7)
				PI=cbind(c(1:nrow(PIm)),c(1:nrow(PIm)),PIm)
				Iint=PI[,2]
				Iint=lapply(Iint,function(i) m.ind[i])
				}

# #########################################
# # take an interaction , take the prior of all models that contain it, take corresponding score and sum the prior times the new score
#########################################
## data frame saved to file 
intposts=NULL
foo=NULL
sig.int=NULL
res.int=NULL
sigint=NULL
if(length(Iint)>0){
for(i in 1:length(Iint))
intposts[i]=sum(priors[Iint[[i]]]*Enewscores[Iint[[i]]])/sum(priors*Enewscores)
intposts=pmin(1,intposts)
if(!is.null(unints)){
if(nrow(unints)>1){
foo=cbind(unints,round(intposts,digits=6))
res.int=foo[order(foo[,(ncol(foo))],decreasing=TRUE),]
sig.int=NULL
if(length(which(res.int[,ncol(foo)]>0))>0)
	{
	is=which(res.int[,ncol(foo)]>0)
	res.int=res.int[is,]
	}else{
			res.int=res.int[-c(1:nrow(res.int)),]
			}
if(nrow(res.int)>0){
					sigint<-res.int[,-ncol(res.int)]
					C1=apply(sigint,1,intpaste)
					sig.int=cbind(C1,res.int[,ncol(res.int)])
					}
					}else{
							foo=c(unints,round(intposts,digits=6))
							res.int=foo#[order(foo[,(ncol(foo))],decreasing=TRUE),]
							sig.int=NULL
							sigint<-res.int[-length(res.int)]
							C1=intpaste(sigint)
							sig.int=c(C1,res.int[length(res.int)])
							}
}
}
if(!is.null(sig.int)){
write.table("# interactions  & posterior probability ", outFile2,append = TRUE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
#if(is.matrix(sig.int))			
write.table(sig.int, outFile2,append = TRUE, quote = FALSE, sep = "\t ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
			}
#########################################
}


if(is.vector(res.add)) sig.int=t(as.matrix(res.add))
if(is.vector(sig.int)) sig.int=t(as.matrix(sig.int))

res.add.ETA<-NULL
if(length(which(res.add[,2]>=0.5))>0) res.add.ETA<-res.add[which(res.add[,2]>=0.5),]
res.int.ETA<-NULL
if(length(which(as.numeric(sig.int[,2])>=0.5))>0) res.int.ETA<-sig.int[which(as.numeric(sig.int[,2])>=0.5),]
tab.add.p.eta=NULL
if(!is.null(res.add.ETA)) if(is.matrix(res.add.ETA)){ 
											tab.add.p.eta=cbind(res.add.ETA[,1],rep(1,nrow(res.add.ETA)),res.add.ETA[,2])
											}else{
											tab.add.p.eta=c(res.add.ETA[1],1,res.add.ETA[2])
											tab.add.p.eta=t(as.matrix(tab.add.p.eta))}
tab.int.p.eta=NULL
if(!is.null(res.int.ETA)) if(is.matrix(res.int.ETA)){ 
											tab.int.p.eta=cbind(res.int.ETA[,1],rep(1,nrow(res.int.ETA)),res.int.ETA[,2])
											}else{
											tab.int.p.eta=c(res.int.ETA[1],1,res.int.ETA[2])
											tab.int.p.eta=t(as.matrix(tab.int.p.eta))
											}

tpfp.ETA=tp.fp.ETA(true.int,tab.add.p=tab.add.p.eta,tab.int.p=tab.int.p.eta,thres=0.15)

TP.add.eta<-tpfp.ETA$TP.add
TP.int.eta<-tpfp.ETA$TP.int
FP.add.eta<-tpfp.ETA$FP.add
FP.int.eta<-tpfp.ETA$FP.int
FDR.a.eta<-tpfp.ETA$FDR.a
FDR.i.eta<-tpfp.ETA$FDR.i
FDR.eta<-tpfp.ETA$FDR
######## cumulative for trees 

resExpr=NULL
if(is.vector(res.int)) res.int=t(as.matrix(res.int))

if(!is.null(sigint)){
Ppost=res.int[,ncol(res.int)]
zm=apply(sigint,1,tvariables.from.model)
# ssigint=t(apply(sigint,1,tree.merge))# sorted trees 
# cr=count.rows(ssigint)
# join trees which contain exactly the same variables and operators
Expr=NULL
sigint1=sigint
#for(ex in 1:length(zm)){
ex=1
# ExprList=NULL
while(ex<(nrow(sigint))){
#if(compare.expr(Expr,sigint[ex,])){ex=ex+2}else{ex=ex+1}
Ni=ex
P=Ppost[ex]
UnL.int=NULL
TreeList=NULL
# PostList=vector("list",nrow(sigint))
for(ni in (ex+1):nrow(sigint1)){
if(length(which(zm[[ex]]%in%zm[[ni]]))==length(zm[[ex]])){
Ni=c(Ni,ni)
P=c(P,Ppost[ni])
}
}
UnL.int=sigint1[ex,]
TreeList=sigint1[Ni,]
PostList=P
################################
## porzadkowanie drzew
if(nrow(TreeList)>0){
untv=TreeList#t(apply(TreeList,1,unique))
if(nrow(TreeList)>1){
# if(!is.list(untv)){
# untv=t(untv)
 # untv=split(untv, rep(1:ncol(untv), each = nrow(untv)))
# }
tv=untv#t(apply(untv,1,abs))
stv=tv#t(apply(tv,1,sort,decreasing=TRUE))
#stv# uzupelnic zerami do drzewa
#n0ind=lapply(stv,function(i) which(i>0))
#lapply(stv,stv[n0ind])
# remove 0s
stvn0=apply(stv,1,rem.null)
}# uzupelnic zerami do drzewa
if(nrow(TreeList)==1){
untv=TreeList
tv=TreeList
stv=tv#sort(tv,decreasing=TRUE)
}
}
lmv=NULL
mv=NULL
if(nrow(stv)>1){
#Lnull=lapply(stv,length)
mv=apply(stv,1,variables.from.model)
lmv=lapply(mv,length)#liczba zmiennych 
Lop=lapply(lmv, function(i) i-1)
sumV=lapply(seq_along(Lop),function(i)
         unlist(Lop[[i]])+unlist(lmv[[i]]))# suma dlugosci operatorow i zmiennych 
Lnull=lapply(seq_along(sumV),function(i)
         7-unlist(sumV[[i]]))# liczba zer do dopisania 
}

if(nrow(stv)==1){
#Lnull=lapply(stv,length)
mv=variables.from.model(stv)
lmv=length(mv)#liczba zmiennych 
Lop=lmv-1
sumV=lmv+Lop#lapply(seq_along(Lop),function(i)
      #   unlist(Lop[[i]])+unlist(lmv[[i]]))# suma dlugosci operatorow i zmiennych 
Lnull=7-sumV
stvn0=rem.null(stv)
}

newTrees=NULL
Cr=NULL
PostTab=NULL
if(nrow(stv)==1){newTrees=stv}
if(nrow(stv)>1){
newTrees=stv
#t(mapply(make.tree,lop=Lop,lnull=Lnull,n0v=stvn0,sumv=sumV,lmv=lmv))
}

if(nrow(newTrees)>2){
#if(nrow(unique(newTrees))>2){
Cr=count.rows(newTrees)
eq.trees=Cr[,2]#which are the same
l.eq.trees=lapply(eq.trees,length)
list.ind=which(as.numeric(unlist(l.eq.trees))>1)# which are repeated
if(length(list.ind)>0){
list.ind.p=eq.trees[list.ind]# go to list index and sum up posteriors with indices from this list inex
Post.list=lapply(list.ind.p,function(i) PostList[i])
#take posteriors of elements from this list
nTsumPost=unlist(as.numeric(lapply(Post.list,sum)))# sum up posteriors in elements of a list
#sum(PostList[eq.trees[list.ind]])
lmv=NULL
mv=NULL
PostTab=NULL
PP=NULL
PP[list.ind]=nTsumPost
rem.ind=c(1:nrow(newTrees))[-unlist(eq.trees[list.ind])]
newTrees1=unique(newTrees)
n.list.ind=c(1:nrow(newTrees1))[-list.ind]
PP[n.list.ind]=PostList[rem.ind]
PostTab=cbind(newTrees1,PP)
mv=apply(newTrees1,1,variables.from.model)
lmv=as.numeric(unlist(lapply(mv,length)))#liczba zmiennych 
}
if(length(list.ind)==0){
newTrees1=newTrees
PostTab=cbind(newTrees,P)
mv=apply(newTrees,1,variables.from.model)
lmv=as.numeric(unlist(lapply(mv,length)))#liczba z
}
}
if(nrow(newTrees)==2){
if(nrow(unique(newTrees))==2){
nTposts=P
nTsumpost=P
newTrees1=newTrees
PostTab=cbind(newTrees,nTsumpost)
# mv=variables.from.model(newTrees)
# lmv=length(mv)
}
}


if(nrow(unique(newTrees))==1){
newTrees1=unique(newTrees)
nTposts=sum(P)
nTsumpost=sum(P)
PostTab=cbind(newTrees,nTsumpost)
mv=variables.from.model(newTrees)
lmv=length(mv)
}

newTrees=newTrees1

## now check if the list of trees contains specific Logic expression within the tree 

# LI=apply(newTrees1,1,tvariables.from.model)
# Lln=lapply(LI,length)
b.ex=newTrees[1,]
p.ex=PostTab[1,8]
Which.cont=1
p.i.2=0
p.i.3=0
p.i.4=0
Lln=NULL
Lln=length(as.numeric(tvariables.from.model(b.ex)))
if(Lln==2) p.i.2=PostTab[1,8]
if(Lln==3) p.i.3=PostTab[1,8]
if(Lln==4) p.i.4=PostTab[1,8]
if(nrow(newTrees)>1){
for(ii in 2:nrow(newTrees)){
op.i=which(sapply(newTrees[ii,],is.operator))# indeks operatora
nop=length(op.i)# number of operators
op=as.numeric(newTrees[ii,which(newTrees[ii,]>=1000)])# operatory
if((length(unique(op))==1)&&(unique(op)==1000)){
# expression contains the specific b.ex
p.ex=p.ex+PostTab[ii,8]
Which.cont=c(Which.cont,ii)
Lli=length(as.numeric(tvariables.from.model(newTrees[ii,])))
if(Lli==2) p.i.2=p.i.2+PostTab[ii,8]
if(Lli==3) p.i.3=p.i.3+PostTab[ii,8]
if(Lli==4) p.i.4=p.i.4+PostTab[ii,8]
}
if(length(unique(op))>1){
if(length(op)==2){
if((op[op.i[1]]==2000)&&(op[op.i[2]]==1000)){
poz=NULL
for(hh in 1:length(zm[[ex]])) poz[hh]=which(newTrees[ii,]==zm[[ex]][hh])
if(sum(poz%in%c(4,5))==2){
p.ex=p.ex+PostTab[ii,8]
Which.cont=c(Which.cont,ii)
p.i.3=p.i.3+PostTab[ii,8]
}
}
}

if(length(op)==3){
if((op[op.i[1]]==2000)&&(op[op.i[2]]==1000)&&(op[op.i[3]]==1000)){
poz=NULL
for(hh in 1:length(zm[[ex]])) poz[hh]=which(newTrees[ii,]==zm[[ex]][hh])
if(sum(poz%in%c(4,5))==2){
p.ex=p.ex+PostTab[ii,8]
Which.cont=c(Which.cont,ii)
p.i.4=p.i.4+PostTab[ii,8]
}
if(sum(poz%in%c(6,7))==2){
p.ex=p.ex+PostTab[ii,8]
Which.cont=c(Which.cont,ii)
p.i.4=p.i.4+PostTab[ii,8]
}
}
}

}
}
}

fin.tab=c(intpaste(b.ex),p.i.2,p.i.3,p.i.4,p.ex)
Expr=sigint[ex,]
ResExpr=fin.tab
resExpr=rbind(resExpr,ResExpr)
ex=ex+1
}
}


resExpr.p=NULL
tab.int.p.ceta=NULL
if(length(which(as.numeric(resExpr[,5])>=0.5))>1){ 
					resExpr.p=resExpr[which(as.numeric(resExpr[,5])>=0.5),]
					tab.int.p.ceta=cbind(resExpr.p[,1],rep(1,nrow(resExpr.p)),resExpr.p[,5])
}
if(length(which(as.numeric(resExpr[,5])>=0.5))==1){ 
					resExpr.p=resExpr[which(as.numeric(resExpr[,5])>=0.5),]
					tab.int.p.ceta=c(resExpr.p[1],1,resExpr.p[5])
					tab.int.p.ceta=t(as.matrix(tab.int.p.ceta))
					}


					
TabA=NULL
TabAdd=NULL
if(length(which(add[,2]>=0.5))>0) TabA=add[which(add[,2]>=0.5),]
if(!is.null(TabA)) TabAdd=cbind(TabA[,1],rep(1,nrow(TabA)),TabA[,2])
#### tp.fp.CETA function
tpfp.ceta=tp.fp.CETA(true.int,tab.add.p=TabAdd,tab.int.p=tab.int.p.ceta,thres=0.15)

TP.add.ceta<-tpfp.ceta$TP.add
TP.int.ceta<-tpfp.ceta$TP.int
FP.add.ceta<-tpfp.ceta$FP.add
FP.int.ceta<-tpfp.ceta$FP.int
FDR.a.ceta<-tpfp.ceta$FDR.a
FDR.i.ceta<-tpfp.ceta$FDR.i
FDR.ceta<-tpfp.ceta$FDR


write.table(add,file=paste(outFile3,"_add.txt"),append = TRUE, quote = FALSE, sep = "\t ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")		
write.table(resExpr,file=paste(outFile3,"_Int.txt"),append = TRUE, quote = FALSE, sep = "\t ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

fourways=NULL
res.int4=NULL
where4way=NULL
if(!is.null(res.int)){
where4way=apply(res.int,1,is.4way)
fourways=res.int[which(where4way),]
}
fourways

TabInt4<-NULL
if(!is.null(fourways)){
if(dim(fourways)[1]>0){
TabInt4<-fourways
TabInt4[,1:3]<-matrix(rep(rep(1000,3),nrow(fourways)),ncol=3)
}
}

if(!is.null(TabInt4)){
Ind=which(as.numeric(TabInt4[,8])>=0.5)
if(length(Ind)>1) TabInt4=TabInt4[which(as.numeric(TabInt4[,8])>=0.5),]
if(length(Ind)==1){
				TabInt4=TabInt4[which(as.numeric(TabInt4[,8])>=0.5),]
				TabInt4=t(as.matrix(TabInt4))
				}
if(length(Ind)==0) TabInt4=NULL
TabInt<-rbind(TabInt,TabInt4)
}


if(!is.null(TabInt)){
if(nrow(TabInt)>1) TabInt<-cbind(TabInt[,1:7],rep(1,nrow(TabInt)),TabInt[,8])
if(nrow(TabInt)==1){
					TabInt<-c(TabInt[,1:7],1,TabInt[,8])
					TabInt<-t(as.matrix(TabInt))
					}
					}
TabInt.p=NULL
if(!is.null(TabInt)){
if(nrow(TabInt)>1) TabInt.p=cbind(intpaste(TabInt[,1:7]),TabInt[,8:9])
if(nrow(TabInt)==1){
					TabInt.p<-c(intpaste(TabInt[,1:7]),TabInt[,8:9])
					TabInt.p<-t(as.matrix(TabInt.p))
					}
					}
					
tpfp.ILA<-tp.fp.ILA(true.int,tab.add.p=TabAdd,tab.int.p=TabInt.p,thres=0.85)
###########################################
TP.add.ila<-tpfp.ILA$TP.add
TP.int.ila<-tpfp.ILA$TP.int
FP.add.ila<-tpfp.ILA$FP.add
FP.int.ila<-tpfp.ILA$FP.int
FDR.a.ila<-tpfp.ILA$FDR.a
FDR.i.ila<-tpfp.ILA$FDR.i
FDR.ila<-tpfp.ILA$FDR


# res.int4=cbind(apply(fourways[,1:7],1,intpaste),fourways[,8])
res.int4=NULL
if(!is.null(fourways)){
				if(nrow(fourways)>1){
				mpom=rbind(matrix(rep(1000,3*nrow(fourways)),ncol=nrow(fourways)),t(as.matrix(abs(fourways[,4:7]))))
				res.int4=cbind(apply(t(mpom),1,intpaste),fourways[,8])
				}
				if(nrow(fourways)==1){
				mpom=rbind(matrix(rep(1000,3*nrow(fourways)),ncol=nrow(fourways)),t(as.matrix(abs(fourways[,4:7]))))
				res.int4=cbind(apply(t(mpom),1,intpaste),fourways[,8])
				}
				if(nrow(fourways)==0){
				res.int4=NULL
				}
				}

#########################################
return(list(add=add, int2=res.int2,int3=res.int3,int4=res.int4,res.add.T=res.add,sig.int.T=sig.int,res.expr=resExpr,TP.add.ila=TP.add.ila,TP.int.ila=TP.int.ila,FP.add.ila=FP.add.ila,FP.int.ila=FP.int.ila,FDR.a.ila=FDR.a.ila,FDR.i.ila=FDR.i.ila,FDR.ila=FDR.ila,TP.add.eta=TP.add.eta,TP.int.eta=TP.int.eta,FP.add.eta=FP.add.eta,FP.int.eta=FP.int.eta,FDR.a.eta=FDR.a.eta,FDR.i.eta=FDR.i.eta,FDR.eta=FDR.eta,TP.add.ceta=TP.add.ceta,TP.int.ceta=TP.int.ceta,FP.add.ceta=FP.add.ceta,FP.int.ceta=FP.int.ceta,FDR.a.ceta=FDR.a.ceta,FDR.i.ceta=FDR.i.ceta,FDR.ceta=FDR.ceta))
}


is.4way=function(vec){
return(length(which(vec[1:3]>=1000))==3)
}


is.3way=function(vec){
return(length(which(vec[1:3]>=1000))==2)
}

is.2way=function(vec){
return(length(which(vec[1:3]>=1000))==1)
}

is.1way=function(vec){
return(length(which(vec[1:3]>=1000))==0)
}


#######################################################
fitmodel = function(nsnp,nrows,nitr,ntrs,nlvs,p,a=0.5,n.tree.mod,tree.lvs,eff.sizes,slope,tree.ind,op)
# single data simulation and model fitting by MCLR
{
	X=snpmatrix(nsnp,nrows,p)
	mod=make.model(X,n.tree.mod,tree.lvs,eff.sizes,slope,tree.ind,op)
	Y=mod$resp(X)
	
	mymccontrol <-logreg.mc.control(nburn=1000, niter=nitr, hyperpars=log(1/a), output=4)
	
	lmod <- logreg(
			resp=binMean(Y),bin=X, 
			type=3,select = 7, 
			ntrees=ntrs,
			nleaves=nlvs*ntrs,
			mc.control=mymccontrol,
			tree.control=logreg.tree.control(treesize=nlvs
			,opers=2)
	)
	return(lmod)
}
 
##########################################################################
binMean<-function (X) 
{
    t = mean(X)
    for (i in seq(1, length(X))) {
        if (X[i] < t) 
            X[i] = 0
        else X[i] = 1
    }
    X
}
###############################################
#############
simulate = function(nsym,nsnp,nrows,nitr,ntrs,nlvs,p,a,lambda,outFile="mypost_results",n.tree.mod,tree.lvs,eff.sizes,slope,tree.ind,op)
## simulation procedure 
# nsym...  no. of different simulated datasets
# nsnp,nrows,nitr,ntrs,nlvs,p,a,n.tree.mod,tree.lvs,eff.sizes,slope,tree.ind,op ... parameters for fitmodel function 
# lambda... Poisson distribution parameter
# outFile ... name of a file to store  results 
{
	
	TP.add.poiss.ila<-TP.int.poiss.ila<-FP.add.poiss.ila<-FP.int.poiss.ila<-TP.add.poiss.eta<-TP.int.poiss.eta<-FP.add.poiss.eta<-FP.int.poiss.eta<-TP.add.poiss.ceta<-TP.int.poiss.ceta<-FP.add.poiss.ceta<-FP.int.poiss.ceta<-vector("list",nsym)
	FDR.a.poiss.ila<-FDR.i.poiss.ila<-FDR.poiss.ila<-FDR.a.poiss.eta<-FDR.i.poiss.eta<-FDR.poiss.eta<-FDR.a.poiss.ceta<-FDR.i.poiss.ceta<-FDR.poiss.ceta<-NULL
	
	TP.add.MCLR.ila<-TP.int.MCLR.ila<-FP.add.MCLR.ila<-FP.int.MCLR.ila<-TP.add.MCLR.eta<-TP.int.MCLR.eta<-FP.add.MCLR.eta<-FP.int.MCLR.eta<-TP.add.MCLR.ceta<-TP.int.MCLR.ceta<-FP.add.MCLR.ceta<-FP.int.MCLR.ceta<-vector("list",nsym)
	FDR.a.MCLR.ila<-FDR.i.MCLR.ila<-FDR.MCLR.ila<-FDR.a.MCLR.eta<-FDR.i.MCLR.eta<-FDR.MCLR.eta<-FDR.a.MCLR.ceta<-FDR.i.MCLR.ceta<-FDR.MCLR.ceta<-NULL
	
	#######################################
	for(i in seq(1,nsym)){
	add.poiss=NULL
	int2.poiss=NULL
	int3.poiss=NULL
	int4.poiss=NULL
	add.poissT=NULL
	int.poissT=NULL
	add.MCLR=NULL
	int2.MCLR=NULL
	int3.MCLR=NULL
	res.expr=NULL
	res=fitmodel(
				nsnp, nrows,
				nitr,
				ntrs,	nlvs,
				p,a,n.tree.mod,tree.lvs,eff.sizes,slope,tree.ind,op)
	#mc.chain.summary(filename="slogiclisting.tmp",OutFile=paste("/",sep=""))
res.poiss=NULL
res.poiss=model.PI.poiss(filename="slogiclisting.tmp",nlvs,ntrs, lambda=lambda,nrows=nrows,outFile=paste(outFile,"_poiss",sep=""),nsnp=nsnp)
		TP.add.poiss.ila[[i]]<-res.poiss$TP.add.ila
		TP.int.poiss.ila[[i]]<-res.poiss$TP.int.ila
		FP.add.poiss.ila[[i]]<-res.poiss$FP.add.ila
		FP.int.poiss.ila[[i]]<-res.poiss$FP.int.ila
		TP.add.poiss.eta[[i]]<-res.poiss$TP.add.eta
		TP.int.poiss.eta[[i]]<-res.poiss$TP.int.eta
		FP.add.poiss.eta[[i]]<-res.poiss$FP.add.eta
		FP.int.poiss.eta[[i]]<-res.poiss$FP.int.eta
		TP.add.poiss.ceta[[i]]<-res.poiss$TP.add.ceta
		TP.int.poiss.ceta[[i]]<-res.poiss$TP.int.ceta
		FP.add.poiss.ceta[[i]]<-res.poiss$FP.add.ceta
		FP.int.poiss.ceta[[i]]<-res.poiss$FP.int.ceta
		FDR.a.poiss.ila<-c(FDR.a.poiss.ila,res.poiss$FDR.a.ila)
		FDR.i.poiss.ila<-c(FDR.i.poiss.ila,res.poiss$FDR.i.ila)
		FDR.poiss.ila<-c(FDR.poiss.ila,res.poiss$FDR.ila)
		FDR.a.poiss.eta<-c(FDR.a.poiss.eta,res.poiss$FDR.a.eta)
		FDR.i.poiss.eta<-c(FDR.i.poiss.eta,res.poiss$FDR.i.eta)
		FDR.poiss.eta<-c(FDR.poiss.eta,res.poiss$FDR.eta)
		FDR.a.poiss.ceta<-c(FDR.a.poiss.ceta,res.poiss$FDR.a.ceta)
		FDR.i.poiss.ceta<-c(FDR.i.poiss.ceta,res.poiss$FDR.i.ceta)
		FDR.poiss.ceta<-c(FDR.poiss.ceta,res.poiss$FDR.ceta)
		
res.MCLR=NULL
res.MCLR=model.PIMCLR(res,nitr,outFile=paste(outFile,sep=""),filename="slogiclisting.tmp",nlvs,ntrs)
		TP.add.MCLR.ila[[i]]<-res.MCLR$TP.add.ila
		TP.int.MCLR.ila[[i]]<-res.MCLR$TP.int.ila
		FP.add.MCLR.ila[[i]]<-res.MCLR$FP.add.ila
		FP.int.MCLR.ila[[i]]<-res.MCLR$FP.int.ila
		TP.add.MCLR.eta[[i]]<-res.MCLR$TP.add.eta
		TP.int.MCLR.eta[[i]]<-res.MCLR$TP.int.eta
		FP.add.MCLR.eta[[i]]<-res.MCLR$FP.add.eta
		FP.int.MCLR.eta[[i]]<-res.MCLR$FP.int.eta
		TP.add.MCLR.ceta[[i]]<-res.MCLR$TP.add.ceta
		TP.int.MCLR.ceta[[i]]<-res.MCLR$TP.int.ceta
		FP.add.MCLR.ceta[[i]]<-res.MCLR$FP.add.ceta
		FP.int.MCLR.ceta[[i]]<-res.MCLR$FP.int.ceta
		FDR.a.MCLR.ila<-c(FDR.a.MCLR.ila,res.MCLR$FDR.a.ila)
		FDR.i.MCLR.ila<-c(FDR.i.MCLR.ila,res.MCLR$FDR.i.ila)
		FDR.MCLR.ila<-c(FDR.MCLR.ila,res.MCLR$FDR.ila)
		FDR.a.MCLR.eta<-c(FDR.a.MCLR.eta,res.MCLR$FDR.a.eta)
		FDR.i.MCLR.eta<-c(FDR.i.MCLR.eta,res.MCLR$FDR.i.eta)
		FDR.MCLR.eta<-c(FDR.MCLR.eta,res.MCLR$FDR.eta)
		FDR.a.MCLR.ceta<-c(FDR.a.MCLR.ceta,res.MCLR$FDR.a.ceta)
		FDR.i.MCLR.ceta<-c(FDR.i.MCLR.ceta,res.MCLR$FDR.i.ceta)
		FDR.MCLR.ceta<-c(FDR.MCLR.ceta,res.MCLR$FDR.ceta)
		
		write.table(c(TP.add.MCLR.ila[[i]],TP.int.MCLR.ila[[i]]), file=paste(outFile,"_TP_ILA_MCLR.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
		write.table(c(TP.add.MCLR.eta[[i]],TP.int.MCLR.eta[[i]]), file=paste(outFile,"_TP_ETA_MCLR.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
		write.table(c(TP.add.MCLR.ceta[[i]],TP.int.MCLR.ceta[[i]]), file=paste(outFile,"_TP_CETA_MCLR.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
		# ########################################
		write.table(c(TP.add.poiss.ila[[i]],TP.int.poiss.ila[[i]]), file=paste(outFile,"_TP_ILA_poiss.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
		write.table(c(TP.add.poiss.eta[[i]],TP.int.poiss.eta[[i]]), file=paste(outFile,"_TP_ETA_poiss.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
		write.table(c(TP.add.poiss.ceta[[i]],TP.int.poiss.ceta[[i]]), file=paste(outFile,"_TP_CETA_poiss.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
		# ########################################
		}
	return(list(TP.add.poiss.ila=TP.add.poiss.ila, TP.int.poiss.ila= TP.int.poiss.ila, FP.add.poiss.ila=FP.add.poiss.ila, FP.int.poiss.ila=FP.int.poiss.ila, TP.add.poiss.eta=TP.add.poiss.eta, TP.int.poiss.eta=TP.int.poiss.eta, FP.add.poiss.eta=FP.add.poiss.eta, FP.int.poiss.eta=FP.int.poiss.eta, TP.add.poiss.ceta=TP.add.poiss.ceta, TP.int.poiss.ceta=TP.int.poiss.ceta, FP.add.poiss.ceta=FP.add.poiss.ceta, FP.int.poiss.ceta=FP.int.poiss.ceta, 
	FDR.a.poiss.ila=FDR.a.poiss.ila, FDR.i.poiss.ila=FDR.i.poiss.ila, FDR.poiss.ila=FDR.poiss.ila, FDR.a.poiss.eta=FDR.a.poiss.eta, FDR.i.poiss.eta=FDR.i.poiss.eta, FDR.poiss.eta=FDR.poiss.eta, FDR.a.poiss.ceta=FDR.a.poiss.ceta, FDR.i.poiss.ceta=FDR.i.poiss.ceta, FDR.poiss.ceta=FDR.poiss.ceta,TP.add.MCLR.ila=TP.add.MCLR.ila, TP.int.MCLR.ila=TP.int.MCLR.ila, FP.add.MCLR.ila=FP.add.MCLR.ila, FP.int.MCLR.ila=FP.int.MCLR.ila, TP.add.MCLR.eta=TP.add.MCLR.eta, TP.int.MCLR.eta=TP.int.MCLR.eta, FP.add.MCLR.eta=FP.add.MCLR.eta, FP.int.MCLR.eta=FP.int.MCLR.eta, TP.add.MCLR.ceta=TP.add.MCLR.ceta, TP.int.MCLR.ceta=TP.int.MCLR.ceta, FP.add.MCLR.ceta=FP.add.MCLR.ceta, FP.int.MCLR.ceta=FP.int.MCLR.ceta, 	FDR.a.MCLR.ila=FDR.a.MCLR.ila, FDR.i.MCLR.ila=FDR.i.MCLR.ila, FDR.MCLR.ila=FDR.MCLR.ila, FDR.a.MCLR.eta=FDR.a.MCLR.eta, FDR.i.MCLR.eta=FDR.i.MCLR.eta, FDR.MCLR.eta=FDR.MCLR.eta, FDR.a.MCLR.ceta=FDR.a.MCLR.ceta, FDR.i.MCLR.ceta=FDR.i.MCLR.ceta, FDR.MCLR.ceta=FDR.MCLR.ceta))	
}




##########################################################################
make.model<-function(X=snpmatrix(nsnp=n.snp,nrows=nobs,p=0.5),n.tree.mod=2,tree.lvs=c(3,3),eff.sizes=c(1.0,1.3,0.5),slope=2.3,tree.ind=c(1,-2,3,23,4,-42),op=c("&","||","||","&"))
## function to create a model 
# n.tree.mod....... number of trees in a true model 
# tree.lvs....number of leaves in the following trees. Must be sorted decreasingly.
# eff.sizes..... effect sizes for trees
# slope.... beta 0
# tree.ind .....indices of variables within the  trees 
# op....... logic operators in the following trees
{
tr<-NULL
I<-NULL
model=NULL
logic.op<-NULL
if(n.tree.mod>0){
tr<-vector("list",n.tree.mod)
I=seq(1:tree.lvs[1]) 
tr[[1]]=tree.ind[I]
if(n.tree.mod>1){
for(i in 2:n.tree.mod){
start.i<-sum(tree.lvs[1:(i-1)])+1
end.i<-sum(tree.lvs[1:i])
tr[[i]]<-tree.ind[start.i:end.i]
}
}

ex.ind<-vector("list", n.tree.mod)# indices of variables within the  trees 
ex.ind=tr
# which are longer than one
tree.int=NULL
for( i in seq(1:n.tree.mod)) tree.int[i]<-length(ex.ind[[i]])
tree.int

logic.op<-NULL
logic.length<-tree.lvs-1 # how many logic operators per tree
if(!sum(logic.length)==0){
logic.op<-vector("list",length(which(tree.int>1)))# logic operators 
I=seq(1:logic.length[1]) 
logic.op[[1]]=op[I]
if(length(which(tree.int>1))>1){
for(i in 2:length(which(tree.int>1))){
start.i<-sum(logic.length[1:(i-1)])+1
end.i<-sum(logic.length[1:i])
logic.op[[i]]=op[start.i:end.i]
}
}
}

Var<-vector("list",n.tree.mod)
for(i in 1:n.tree.mod){
var.vec<-NULL
if(length(which(ex.ind[[i]]<0))==0){
var.vec<-paste(" X[, ",ex.ind[[i]][1],"] ",sep="")
if(!logic.length[i]==0){
for(k in 1:length(logic.op[[i]])){
var.vec=c(var.vec,paste(logic.op[[i]][k]," "," X[, ",ex.ind[[i]][k+1],"] ",sep=""))
}
var.vec<-paste(var.vec,collapse="")
}
}


if(length(which(ex.ind[[i]]<0))>0){
var.vec<-ifelse(ex.ind[[i]][1]<0,paste("!X[, ",abs(ex.ind[[i]][1]),"] ",sep=""),paste("X[, ",ex.ind[[i]][1],"] ",sep=""))
if(!logic.length[i]==0){
for(k in 1:length(logic.op[[i]])){
v1<-c(var.vec,paste(logic.op[[i]][k]," "," !X[, ",abs(ex.ind[[i]][k+1]),"] ",sep=""))
v11<-paste(v1,collapse="")
v2<-c(var.vec,paste(logic.op[[i]][k]," "," X[, ",abs(ex.ind[[i]][k+1]),"] ",sep=""))
v21<-paste(v2,collapse="")
var.vec<-ifelse(ex.ind[[i]][k+1]<0,v11,v21)
}
}
}
Var[[i]]<-var.vec
}
Var
Var.vec<-NULL
Var.vec<-unlist(Var)
Var.vec

mod.vec<-NULL
for(i in 1:n.tree.mod)  mod.vec[[i]]<-paste(eff.sizes[i],"*(",Var.vec[i],")", sep="")
mod.vec
Mod.vec<-NULL
Mod.vec<-paste(mod.vec,collapse="+")
model=NULL
model=paste(c(slope,Mod.vec),collapse="+")
}
#######snp model
n.snp=ncol(X)/2
###############################################
namesL<-NULL
ind1<-2*c(1:n.snp)-1
ind2<-2*c(1:n.snp)
i=c(1:n.snp)
namesL[ind1]<-paste("D_{",i,"}",sep="")
namesL[ind2]<-paste("R_{",i,"}",sep="")
################ model int 
Mnam<-NULL
for(i in 1:length(namesL))
  Mnam[i]<-paste(as.name(paste(namesL)[i]))
  
snpmodel=NULL
if(n.tree.mod>0){
SnpVar<-vector("list",n.tree.mod)
for(i in 1:n.tree.mod){
snpvar.vec<-NULL
if(length(which(ex.ind[[i]]<0))==0){
snpvar.vec<-paste(Mnam[ex.ind[[i]][1]])
if(!logic.length[i]==0){
for(k in 1:length(logic.op[[i]])){
snpvar.vec=c(snpvar.vec,paste(logic.op[[i]][k]," ",Mnam[ex.ind[[i]][k+1]],sep=""))
}
snpvar.vec<-paste(snpvar.vec,collapse="")
}
}



if(length(which(ex.ind[[i]]<0))>0){
snpvar.vec<-ifelse(ex.ind[[i]][1]<0,paste(Mnam[abs(ex.ind[[i]][1])],"^C",sep=""),paste(Mnam[ex.ind[[i]][1]],sep=""))
if(!logic.length[i]==0){
for(k in 1:length(logic.op[[i]])){
v1<-c(snpvar.vec,paste(logic.op[[i]][k]," ",Mnam[abs(ex.ind[[i]][k+1])],"^C ",sep=""))
v11<-paste(v1,collapse="")
v2<-c(snpvar.vec,paste(logic.op[[i]][k]," ",Mnam[abs(ex.ind[[i]][k+1])],sep=""))
v21<-paste(v2,collapse="")
snpvar.vec<-ifelse(ex.ind[[i]][k+1]<0,v11,v21)
}
}
}
SnpVar[[i]]<-snpvar.vec
}
SnpVar

SnpVar.vec<-NULL
SnpVar.vec<-unlist(SnpVar)
SnpVar.vec

snpmod.vec<-NULL
for(i in 1:n.tree.mod)  snpmod.vec[[i]]<-paste(eff.sizes[i],"*(",SnpVar.vec[i],")", sep="")
snpmod.vec
SnpMod.vec<-NULL
SnpMod.vec<-paste(snpmod.vec,collapse="+")
snpmodel=NULL
snpmodel=paste(c(slope,SnpMod.vec),collapse="+")
}

resp=function(X){
expr0<-paste("Y=",model," + rnorm(X[,1],mean=0,sd=1)")
cmd<-paste(expr0)
return(eval(parse(text=cmd)))
}
return(list(resp=resp,model=model,tr=tr, logic.op=logic.op,#single=Single,double=Double,triple=Triple,fourway=Fourway,
snpmodel=snpmodel))
}

# ##########################################################################


########################################################################
model_int<-function(X=snpmatrix(nsnp=n.snp,nrows=nobs,p=0.5),nlvs=n.lvs,n.tree.mod=3,tree.lvs=c(3,3,3),eff.sizes=c(1.0,1.3,0.5),slope=2.3,tree.ind=c(1,-2,3,23,4,-42),op=c("&","&","&"))
# what are model interactions 
# returned in a tree-structured notation (1000-2000 coding)
{
   m.int=NULL
   Mm=make.model(X,n.tree.mod,tree.lvs,eff.sizes,slope,tree.ind,op)
  tr=Mm$tr
  logic.op=Mm$logic.op
  logicop=vector('list',length(tr))
  nullvec=vector('list',length(tr))# vector of zeros to build a complete tree		
if(length(tr)>0){
for(i in 1:length(tr)){
	if(length(logic.op)>0){
		for(j in 1:length(logic.op)){
		Lo=NULL
		if(length(logic.op[[j]])>0){
		for(ss in 1:length(logic.op[[j]])){
		if(logic.op[[j]][ss]=="&") Lo<-c(Lo,1000)
		if(logic.op[[j]][ss]=="||") Lo<-c(Lo,2000)
		}
		logicop[[j]]=Lo
		if(length(tr[[j]])+length(logic.op[[j]])<nlvs)
		nullvec[[j]]=rep(0,nlvs-(length(tr[[j]])+length(logic.op[[j]])))
		}
		}		
		v=tr[[i]]
		v=c(logicop[[i]],v,nullvec[[i]])
		if(length(v)>1) v=tree.merge(v)
		m.int[i]=intpaste(v)
		}		
	if(length(tr[[i]])==1){				
		v=c(tr[[i]],rep(0,6))
		m.int[i]=intpaste(v)
		}
		}
		}		
  return(list(m.int=m.int))#,m.int2=m.int2))
}


##########################################################################
snpmatrix=function(nsnp,nrows,p)
## simulates a matrix of snps containing nsnp snps, nrows observations (individuals) and a minor allele frequency of p 
{
	ncols=2*nsnp
	X=matrix(rbinom(ncols*nrows,1,p),nrows,ncols)
	
	for(i in 1:nrows)
		for(j in 1:nsnp)
			if(X[i,2*j-1]==0 && X[i,2*j]==1)
			{
				X[i,2*j]=0
				X[i,2*j-1]=1
			}
	
	X
}

##########################################################################
tree.to.int=function(tree=tree,n.snp=n.snp)
### function to create interaction names from a tree structure
{
###############################################
namesL<-NULL
ind1<-2*c(1:n.snp)-1
ind2<-2*c(1:n.snp)
i=c(1:n.snp)
namesL[ind1]<-paste("D_{",i,"}",sep="")
namesL[ind2]<-paste("R_{",i,"}",sep="")
################ model int 
Mnam<-NULL
for(i in 1:length(namesL))
  Mnam[i]<-paste(as.name(paste(namesL)[i])) # names needed to create interactions
# ################################################
nlvs=length(tree) # number of nodes in a tree 
ntlvs=length(which(sapply(tree,is.leaf)))# number of leaves
nop=length(which(sapply(tree,is.operator)))# number of operators

Nlvs=0
i=1
fname=NULL
while(Nlvs<ntlvs){
	if(is.leaf(tree[i])){
		fname=c(fname,ifelse(tree[i]<0,paste(Mnam[abs(tree[i])],"^C",sep=""),paste(Mnam[tree[i]],sep="")))
		Nlvs=Nlvs+1
						}
	else{
	if(is.operator(tree[i]))
		if((is.leaf(tree[2*i]))&&(is.leaf(tree[2*i+1]))){
		optopaste=ifelse(tree[i]==1000,':',':')
		snp1=ifelse(tree[2*i]<0,paste(Mnam[abs(tree[2*i])],"^C",sep=""),paste(Mnam[tree[2*i]],sep=""))
		snp2=ifelse(tree[2*i+1]<0,paste(Mnam[abs(tree[2*i+1])],"^C",sep=""),paste(Mnam[tree[2*i+1]],sep=""))
		fname=c(fname,paste(snp1,optopaste,snp2,sep=" "))
		Nlvs=Nlvs+2
		i=i+1
		} else {i=i+1}
		}
				}
				
if(length(fname)>1){
k=0
fnam=NULL
while(k<(length(fname)-1)){
		fnam=c(fnam,paste("(",fname[length(fname)-k],") ",ifelse(tree[k+1]==1000,':',':')," (",fname[length(fname)-(k+1)],") ",sep=""))
		k=k+1
						}
						} else {fnam=fname}
						
return(fnam)
}


##########################################################################
tree.to.int_for_plots=function(tree=tree,n.snp=n.snp)
### function to create interaction names from a tree structure
{
###############################################
namesL<-NULL
ind1<-2*c(1:n.snp)-1
ind2<-2*c(1:n.snp)
i=c(1:n.snp)
namesL[ind1]<-paste("D[",i,"]",sep="")
namesL[ind2]<-paste("R[",i,"]",sep="")
################ model int 
Mnam<-NULL
for(i in 1:length(namesL))
  Mnam[i]<-paste(as.name(paste(namesL)[i])) # names needed to create interactions
# ################################################
nlvs=length(tree) # number of nodes in a tree 
ntlvs=length(which(sapply(tree,is.leaf)))# number of leaves
nop=length(which(sapply(tree,is.operator)))# number of operators

Nlvs=0
i=1
fname=NULL
while(Nlvs<ntlvs){
	if(is.leaf(tree[i])){
		fname=c(fname,ifelse(tree[i]<0,paste(Mnam[abs(tree[i])],"^C",sep=""),paste(Mnam[tree[i]],sep="")))
		Nlvs=Nlvs+1
						}
	else{
	if(is.operator(tree[i]))
		if((is.leaf(tree[2*i]))&&(is.leaf(tree[2*i+1]))){
		optopaste=ifelse(tree[i]==1000,':',':')
		snp1=ifelse(tree[2*i]<0,paste(Mnam[abs(tree[2*i])],"^C",sep=""),paste(Mnam[tree[2*i]],sep=""))
		snp2=ifelse(tree[2*i+1]<0,paste(Mnam[abs(tree[2*i+1])],"^C",sep=""),paste(Mnam[tree[2*i+1]],sep=""))
		fname=c(fname,paste(snp1,optopaste,snp2,sep=" "))
		Nlvs=Nlvs+2
		i=i+1
		} else {i=i+1}
		}
				}
				
if(length(fname)>1){
k=0
fnam=NULL
while(k<(length(fname)-1)){
		fnam=c(fnam,paste("(",fname[length(fname)-k],") ",ifelse(tree[k+1]==1000,':',':')," (",fname[length(fname)-(k+1)],") ",sep=""))
		k=k+1
						}
						} else {fnam=fname}
						
return(fnam)
}


tree.to.int_Logic=function(tree=tree,n.snp=n.snp)
### function to create interaction names  with logic operators from a tree structure
{
###############################################
namesL<-NULL
ind1<-2*c(1:n.snp)-1
ind2<-2*c(1:n.snp)
i=c(1:n.snp)
namesL[ind1]<-paste("D_{",i,"}",sep="")
namesL[ind2]<-paste("R_{",i,"}",sep="")
################ model int 
Mnam<-NULL
for(i in 1:length(namesL))
  Mnam[i]<-paste(as.name(paste(namesL)[i])) # names needed to create interactions
# ################################################
nlvs=length(tree) # number of nodes in a tree 
ntlvs=length(which(sapply(tree,is.leaf)))# number of leaves
nop=length(which(sapply(tree,is.operator)))# number of operators

Nlvs=0
i=1
fname=NULL
while(Nlvs<ntlvs){
	if(is.leaf(tree[i])){
		fname=c(fname,ifelse(tree[i]<0,paste(Mnam[abs(tree[i])],"^C",sep=""),paste(Mnam[tree[i]],sep="")))
		Nlvs=Nlvs+1
						}
	else{
	if(is.operator(tree[i]))
		if((is.leaf(tree[2*i]))&&(is.leaf(tree[2*i+1]))){
		optopaste=ifelse(tree[i]==1000,' * ',' | ')
		snp1=ifelse(tree[2*i]<0,paste(Mnam[abs(tree[2*i])],"^C",sep=""),paste(Mnam[tree[2*i]],sep=""))
		snp2=ifelse(tree[2*i+1]<0,paste(Mnam[abs(tree[2*i+1])],"^C",sep=""),paste(Mnam[tree[2*i+1]],sep=""))
		fname=c(fname,paste(snp1,optopaste,snp2,sep=" "))
		Nlvs=Nlvs+2
		i=i+1
		} else {i=i+1}
		}
				}
				
if(length(fname)>1){
k=0
fnam=NULL
while(k<(length(fname)-1)){
		fnam=c(fnam,paste("(",fname[length(fname)-k],") ",ifelse(tree[k+1]==1000,' * ',' | ')," (",fname[length(fname)-(k+1)],") ",sep=""))
		k=k+1
						}
						} else {fnam=fname}
						
return(fnam)
}

##########################################################################
tex.tab=function(vec)
## function to rewrite a vector (vec) onto a row in a tex-formated table
{
res=paste(vec, collapse=" & ")
res=paste(res, '\\ \\hline',sep="")
return(res)
}


##########################################################################
int.lev=function(vec)
## level of interaction in a tree (vec) 
{
Lop=length(which(vec>=1000))
return(7-(length(which(vec==0))+Lop))
}
##########################################################################

summarize.ILA=function(OutFile,A,p.param,n.sym,true.int,object,n.snp,n.obs,n.itr,boxplots.all=FALSE)
## function to summarize results of simulations for ILA approach
## files containing tables in tex-syntax are created
{
fs=file.info(paste(OutFile,".txt",sep=""))$size
if(fs>0){		
		M1.1<-NULL
		M1.1<-read.table(file=paste(OutFile,".txt",sep=""), header = FALSE, sep = "\t", quote = "\"'", dec = ".",as.is = TRUE,na.strings = "NA", colClasses = NA, nrows = -1,
				   skip = 0, check.names = TRUE, strip.white = FALSE, blank.lines.skip = TRUE,
				   comment.char = "#",
				   allowEscapes = FALSE, flush = FALSE,
				   stringsAsFactors = FALSE,
				   fileEncoding = "", encoding = "unknown")
add.eff<-NULL
int2<-NULL
int3<-NULL
str2 <- gsub(' {2,}',' ',M1.1[,1])
int.level=unlist(lapply(strsplit(str2,' '),length))
add.eff<-M1.1[which(int.level==1),]
int=M1.1[which(int.level>1),]
res.post.add=NULL
if(nrow(add.eff)>0){
Lv=count.rows(add.eff[,1])
L.p=Lv[,2]
res.post.add=vector("list",length(L.p))
for(i in 1:length(L.p)){
p=NULL
if(length(which(Lv[,3]==i))>0) p=add.eff[L.p[which(Lv[,3]==i)][[1]],2]
res.post.add[[i]]=list(v=as.numeric(Lv[i,3]),n=Lv[i,1],p=p)
}
}
#####
write.table(int, file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
Alli=read.table(file = "M.int.txt")
fn <- "M.int.txt"
if (file.exists(fn)) file.remove(fn)
#####
Alli1=t(apply(Alli[,1:7],1,tree.merge))
Alli1.vec=apply(Alli1[,1:7],1,intpaste)
int.1=cbind(Alli1.vec,int[,2])
int=int.1
#####
res.post.int=NULL
Lvi=NULL
if(nrow(int)>0){
Lvi=count.rows(int[,1])
res.post.int=vector("list",nrow(Lvi))
res.post.int=apply(Lvi,1,do.list,int=int)
}
###############################################
namesL1<-NULL
ind1<-2*c(1:n.snp)-1
ind2<-2*c(1:n.snp)
i=c(1:n.snp)
namesL1[ind1]<-paste("D[",i,"]",sep="")
namesL1[ind2]<-paste("R[",i,"]",sep="")
################ BOXPLOTS
mint.v<-NULL
Ll=NULL
if(!is.null(true.int)){
write.table(true.int, file= "trint.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
mint.v<-read.table(file = "trint.txt")
fn <- "trint.txt"
if (file.exists(fn)) file.remove(fn)
va=t(apply(mint.v,1,variables.from.model))
vva=sort(as.vector(unlist(va)))
TNa<-ifelse(vva<0,paste(namesL1[abs(vva)],"^C",sep=""),paste(namesL1[as.numeric(abs(vva))],sep=""))
ex.nam<-TNa

if(!is.null(res.post.add)){
Ll=res.post.add[vva]
if(boxplots.all){
for(ss in 1:length(vva)){
jpeg(file=file.path(OutFile,"_boxplot_add",ss,".jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
allps=Ll[[ss]]$p
if(length(Ll[[ss]]$p)<n.sym) allps=c(Ll[[ss]]$p,rep(0,(n.sym-length(Ll[[ss]]$p))))
boxplot(allps,
        col = "yellow",
        main = paste("Boxplot of posterior probability for",ex.nam[ss]),
        xlab = " ",
        ylab = "Posterior probability",
        ylim = c(0, 2), yaxs = "i")
		axis(1, c(1), paste(ex.nam[ss],": Number of appearances :",Ll[[ss]]$n,", Number of post.pr >=0.5 :",length(which(allps>=0.5)),sep=""), cex.lab=1.2, cex.axis=0.4, cex.main=0.5, cex.sub=0.5, col.axis = "black",las=1)
dev.off()
}
}
}
}
####### boxplots-altogether
ps=NULL
res=NULL
for(ss in 1:length(vva)){
ps=Ll[[ss]]$p
if(length(Ll[[ss]]$p)<n.sym) ps=c(Ll[[ss]]$p,rep(0,(n.sym-length(Ll[[ss]]$p))))
res=cbind(res,ps)
}
colnames(res)=NULL
if(nrow(res)>1){
jpeg(file=file.path(OutFile,"_boxplot_add_all.jpeg",fsep=""),height=12,width=24,unit="cm",res=300)
boxplot(res[,1:ncol(res)], 
        names=paste(ex.nam[1:ncol(res)],sep=""),
        col = "yellow",
        main = paste("Boxplots of posterior probabilities"),
        xlab = " ",
        ylab = "Posterior probability",
        ylim = c(0, 1.2), yaxs = "i",las=2)
dev.off()
}
Ll=NULL


###############################################################
##### indices of models containing true interactions
# ################################
if(!is.null(true.int)){
 write.table(true.int, file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  true.int1=read.table(file = "M.int.txt")
   fn <- "M.int.txt"
 if (file.exists(fn)) file.remove(fn)
 }
  Alli=NULL
  ind.tr.i=NULL
   if(!is.null(res.post.int)){
   write.table(Lvi[,3], file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  Alli=read.table(file = "M.int.txt")
   fn <- "M.int.txt"
 if (file.exists(fn)) file.remove(fn)
 }
 ind.tr.i=NULL
  if(!is.null(true.int1)){
  if(length(unique(unlist(lapply(apply(true.int1,1,variables.from.model),length))))==1){
  v=t(apply(true.int1,1,variables.from.model))
  if(nrow(v)>1){
  NO.Int=NULL
  NO.Int=nrow(true.int1)
  ind.tr.i=vector("list",nrow(true.int1))
  for(ij in 1:nrow(true.int1)){
  ind.tr.i[[ij]]=0
  sumi=which(rowSums(Alli)==sum(true.int1[ij,]))
  m1=Alli[which(rowSums(Alli)==sum(true.int1[ij,])),]
  kt=NULL
  for(rr in 1:nrow(m1)) kt=c(kt,(sum(v[ij,]%in%m1[rr,])==length(v[ij,])))
  if(length(which(kt))>0) ind.tr.i[[ij]]=sumi[which(kt)]
  }}
  }else{v=apply(true.int1,1,variables.from.model)# a list
		NO.Int=NULL
		NO.Int=length(which(unlist(lapply(v,length))>1))
		ind.tr.i=vector("list",NO.Int)
		for(ij in 1:NO.Int){
		  ind.tr.i[ij]=0
		  sumi=which(as.vector(rowSums(Alli)==sum(true.int1[ij,])))
		  m1=Alli[which(rowSums(Alli)==sum(true.int1[ij,])),]
		  kt=NULL
		  for(rr in 1:nrow(m1)) kt=c(kt,(sum(v[[ij]]%in%m1[rr,])==length(v[[ij]])))
		  if(length(which(kt))>0) ind.tr.i[[ij]]=sumi[which(kt)]
		  }
	  }
  }
##### boxplots for interactions : individually
  m1.int.v<-NULL
  ex.nam.int<-NULL
  if(!is.null(true.int)){
  write.table(true.int, file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  m1.int.v<-read.table(file = "M.int.txt")
  fn <- "M.int.txt"
  if (file.exists(fn)) file.remove(fn)
  ex.nam.int<-apply(m1.int.v,1,tree.to.int_for_plots,n.snp=n.snp)
  }
if(boxplots.all){ 
if(!is.null(ind.tr.i)){ 
for(si in 1:NO.Int){ 
if(ind.tr.i[[si]]>0){
jpeg(file=file.path(OutFile,"_boxplot_int",si,".jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
allpints=res.post.int[[ind.tr.i[[si]]]]$p
if(length(res.post.int[[ind.tr.i[[si]]]]$p)<n.sym) allpints=c(res.post.int[[ind.tr.i[[si]]]]$p,rep(0,(n.sym-length(res.post.int[[ind.tr.i[[si]]]]$p))))
boxplot(allpints,
		col = "yellow",
        ylab = "Posterior probability",
        ylim = c(0, 1.5), yaxs = "i")
		axis(1, c(1),ex.nam.int[si], cex.lab=1.2, cex.axis=0.4, cex.main=0.5, cex.sub=0.5, col.axis = "black",las=2)
dev.off()
}
}
}
}
### boxplots -all int together
  resint=NULL
  for(si in 1:NO.Int){ 
  pints=NULL
  if(ind.tr.i[[si]]>0){
  pints=res.post.int[[ind.tr.i[[si]]]]$p
  if(length(res.post.int[[ind.tr.i[[si]]]]$p)<n.sym) pints=c(pmin(1,res.post.int[[ind.tr.i[[si]]]]$p),rep(0,(n.sym-length(res.post.int[[ind.tr.i[[si]]]]$p))))
  resint=cbind(resint,pints)
  }
  colnames(resint)=NULL
  }
  
if(length(which(unlist(ind.tr.i)>0))>0){ 
jpeg(file=file.path(OutFile,"_boxplot_int_all.jpeg",fsep=""),height=24,width=24,unit="cm",res=300)
par(mar=c(12,8,2,2)+0.3,mgp=c(3,1,0))
boxplot(resint[,1:ncol(resint)],
		names=paste(ex.nam.int[which(ind.tr.i>0)],sep=""),
		col = "yellow",
		main = paste("Boxplots of posterior probabilities"),
        ylab = "Posterior probability",
        ylim = c(0, 1.1), yaxs = "i", las=2)
dev.off()
}
# ################### BOXPLOT OF TREE SIZES
i=NULL
s=NULL
ex.nam=NULL
ex.nam.int=NULL
m1.int.v=NULL
##################
tab.add<-NULL
if(length(res.post.add)>0){
tab.add<-matrix(nrow=length(res.post.add),ncol=3)
for(i in 1:length(res.post.add)){ 
tab.add[i,]<-c(res.post.add[[i]]$v,res.post.add[[i]]$n,sum(as.numeric(res.post.add[[i]]$p)))
}
}

tab.int<-NULL
if(!is.null(res.post.int)){
tab.int<-matrix(nrow=length(res.post.int),ncol=3)
for(i in 1:length(res.post.int)){
tab.int[i,]<-c(res.post.int[[i]]$v$x,res.post.int[[i]]$n$counts,sum(as.numeric(res.post.int[[i]]$p)))
}
}

 tab.add[,3]<-as.numeric(tab.add[,3])/n.sym
 tab.int[,3]<-as.numeric(tab.int[,3])/n.sym 
 if(!is.null(tab.add)) 
 if(length(which(as.numeric(tab.add[,3])>=0))>0)
	{tab.add=tab.add[which(as.numeric(tab.add[,3])>=0),]
	}else{tab.add=NULL}
 if(!is.null(tab.int)) 
 if(length(which(as.numeric(tab.int[,3])>=0))>0){
 tab.int=tab.int[which(as.numeric(tab.int[,3])>=0),]}else{tab.int=NULL}
 
 ind=NULL
 No.a=NULL
for(jj in 1:length(res.post.add)){ 
  if(length(which(res.post.add[[jj]]$p>=0.5))>0){
  ind=c(ind,jj)
  No.a=c(No.a,length(which(res.post.add[[jj]]$p>=0.5)))
  }
 }
 tab.add.p=cbind(tab.add[ind,1],No.a,tab.add[ind,3])
 N.A=NULL
 N.A=sum(No.a)
 
 ind=NULL
 No.i=NULL
if(!is.null(res.post.int)){
for(jj in 1:length(res.post.int)){ 
  if(length(which(res.post.int[[jj]]$p>=0.5))>0){
  ind=c(ind,jj)# which interaction had at lest once a posterior larger than 0.5
  No.i=c(No.i,length(which(res.post.int[[jj]]$p>=0.5)))
  }
 }
 tab.int.p=cbind(tab.int[ind,1],No.i,tab.int[ind,3])
colnames(tab.int.p)=NULL
 tab.add.red=tab.add.p
 tab.int.red=tab.int.p
 Li=NULL
 if(length(tab.int.p)>0){
 write.table(tab.int.p[,1], file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
 Li=read.table(file = paste("M.int.txt",sep=""))
 fn <- "M.int.txt"
 if (file.exists(fn)) file.remove(fn)
 }
 i.lev=NULL
if(!is.null(Li)) i.lev=apply(Li,1,int.lev)
N.2=0
N.3=0
N.4=0
 if(length(which(i.lev==2))>0) N.2=sum(No.i[which(i.lev==2)])
  if(length(which(i.lev==3))>0) N.3=sum(No.i[which(i.lev==3)])
  if(length(which(i.lev==4))>0) N.4=sum(No.i[which(i.lev==4)])
}  
dist.lvs=c(rep(1,N.A),rep(2,N.2),rep(3,N.3),rep(4,N.4))

jpeg(file=file.path(OutFile,"_BoxplotTreeSizes.jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
boxplot(dist.lvs,
        col = "yellow")
dev.off()
########### HISTOGRAM OF TREE SIZES
jpeg(file=file.path(OutFile,"_HistTreeSizes.jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
hist(dist.lvs,col = "lightblue", border = "pink",main="Histogram of selected trees size",xlab="Tree size")
dev.off()

# ##################################################
# ##########   FP and TP 
#

# ##########   FP and TP 
res.tpfp<-tp.fp.ILA(true.int,tab.add.p,tab.int.p,thres=0.85)
TP.add=res.tpfp$TP.add
TP.int=res.tpfp$TP.int
FP.add=res.tpfp$FP.add
FP.int=res.tpfp$FP.int
tp.fp.ILA.plot(true.int,tab.add.p,tab.int.p,max.thres=0.85,OutFile=OutFile)

# # #### FALSE DISCOVERY RATE
FDR=res.tpfp$FDR
FDR.a=res.tpfp$FDR.a
FDR.i=res.tpfp$FDR.i


###############################################
#######		PLOTS AND TABLES WITH SUMMARIES
###############################################
namesL<-NULL
ind1<-2*c(1:n.snp)-1
ind2<-2*c(1:n.snp)
i=c(1:n.snp)
namesL[ind1]<-paste("D_{",i,"}",sep="")
namesL[ind2]<-paste("R_{",i,"}",sep="")
#####################################################
######## table with settings to the file 
#####################################################
set.tab=NULL
set.tab<-rbind(c(paste('True model:'),paste("$",object$snpmodel,"$")),
c(paste("No. of snps"), n.snp),c(paste("No. of MC iterations"),n.itr),c(paste("'a' parameter"),A),c(paste("'p' parameter"), p.param),c(paste("No of simulations"),n.sym),c(paste("No of individuals"),n.obs),c(paste("Max no. of leaves per tree"), (n.lvs+1)/2),c(paste("Max no. of trees"), n.trs))
newset.tab=as.matrix(apply(set.tab,1,tex.tab))
write.table(paste("\\begin{table}[!h] \\caption{Settings }\\label{} \\centering \\begin{tabular}{p{4cm}|p{8cm}}\\hline"), file = paste(OutFile,"_tab.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(newset.tab, file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"),file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\begin{table}[!h]\\caption{Selected expressions}\\label{}\\centering\\begin{tabular}{|c|c|c|c|}\\hline"), file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
tp.add=TP.add[,1]
fp.add=FP.add[,1]
library(plotrix)
Intnames<- c(tp.add,fp.add)
PostP<-NULL
if(!is.null(tab.add.p)){
if(is.matrix(tab.add.p)) PostP<-c(as.numeric(tab.add.p[,3]))
if(!is.matrix(tab.add.p)) PostP<-c(as.numeric(tab.add.p[3]))
}
T1<- NULL
T1<-c(tp.add,fp.add)
T1<-as.numeric(T1)
TNa<-ifelse(T1<0,paste(namesL[abs(T1)],"^C",sep=""),paste(namesL[as.numeric(abs(T1))],sep=""))
ex.nam<-TNa
tabnew_add=NULL
##########################################################
write.table(c(paste('Expression & No.  & Frequency & Posterior probability\\ \\hline')), file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
if(length(ex.nam)>0){
 if(is.matrix(tab.add.p)) tabnew_add=cbind(paste("$",ex.nam,"$"),c(TP.add[,2],FP.add[,2]),round(as.numeric(c(TP.add[,2],FP.add[,2]))/n.sym,digits=5),c(TP.add[,3],FP.add[,3]))
 if(!is.matrix(tab.add.p)) tabnew_add=cbind(paste("$",ex.nam,"$"),c(TP.add[2],FP.add[2]),round(as.numeric(c(TP.add[2],FP.add[2]))/n.sym,digits=5),c(TP.add[3],FP.add[3]))
colnames(tabnew_add)=NULL
tabnew_add1=NULL
if(is.matrix(tabnew_add)){
if(nrow(tabnew_add)>0)
tabnew_add1=as.matrix(apply(tabnew_add,1,tex.tab))
}
if(!is.matrix(tabnew_add)){
if(length(tabnew_add)>0)
tabnew_add1=tex.tab(tabnew_add)
}
write.table(tabnew_add1, file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}
############# tex tables with fp and tp for add effects
l.tp.add=0
l.fp.add=0
if(!is.null(TP.add)) l.tp.add=nrow(TP.add)
if(!is.null(FP.add)) l.fp.add=nrow(FP.add)
write.table(paste("\\begin{table}[!h]\\caption{True positives}\\label{}\\centering\\begin{tabular}{|c|c|c|c|}\\hline"), file = paste(OutFile,"_tp_tab.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(c(paste("Expression & No.  & Frequency & Posterior probability\\ \\hline")), file = paste(OutFile,"_tp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
if(l.tp.add>0)
write.table(tabnew_add1[1:l.tp.add,], file = paste(OutFile,"_tp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

write.table(paste("\\begin{table}[!h]\\caption{False positives}\\label{}\\centering\\begin{tabular}{|c|c|c|c|}\\hline"), file = paste(OutFile,"_fp_tab.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(c(paste('Expression & No.  & Frequency & Posterior probability\\ \\hline')), file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
if(l.tp.add==0){
if(l.fp.add>0)
write.table(tabnew_add1[1:l.fp.add,], file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))}else{
if(l.fp.add>0)
write.table(tabnew_add1[(l.tp.add+1):l.fp.add,], file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}
##########################################################
m.tp.int=NULL
m.fp.int=NULL
if(!is.null(TP.int)) if(nrow(TP.int)>0) m.tp.int=TP.int[,1:7]
if(!is.null(FP.int)) if(nrow(FP.int)>0) m.fp.int=FP.int[,1:7]
# m.int<-rbind(m.tp.int,m.fp.int)
m.int<-NULL
if(!is.null(m.tp.int)&&!is.null(m.fp.int)){m.int=rbind(as.matrix(m.tp.int),as.matrix(m.fp.int))}else{m.int=mat.bind(m.tp.int,m.fp.int)}


if(!is.null(m.int)){ 
  PP<-NULL
  if(!is.null(nrow(tab.int.p))){ 
  for(i in 1:nrow(tab.int.p)) PP[i]<-as.double(as.numeric(tab.int.p[i,3]))
  PostP<-c(PostP,PP)
  m1.int.v<-NULL
  write.table(tab.int.p[,1], file = "al.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  }
  if(is.null(nrow(tab.int.p))){
  PP<-as.double(as.numeric(tab.int.p[3]))
  PostP<-c(PostP,PP)
  m1.int.v<-NULL
  write.table(tab.int.p[1], file = "al.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  }
  m1.int.v<-read.table(file = "al.int.txt")
  fn <- "al.int.txt"
  if (file.exists(fn)) file.remove(fn)
  ex.nam<-c(ex.nam,apply(m1.int.v,1,tree.to.int,n.snp=n.snp))
}

b.s<-NULL
b.s<-cbind(ex.nam,PostP)
if(nrow(b.s)>0){
jpeg(file=file.path(OutFile,"_Posterior_prob.jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
xb<-c(1:nrow(b.s))
plot(x=xb,y=PostP,lty=1,lwd=2,type="h",axes=FALSE,xlab=NA,ylab="Posterior probability",
main=paste("Posterior probability for logic expressions.\nTrue model:",object$snpmodel),
ylim=c(0.0,1.5),cex.lab=1.0, cex.axis=0.5, cex.main=0.4, cex.sub=0.5,col="blue")
    axis(1, xb, ex.nam, cex.lab=1.2, cex.axis=0.4, cex.main=0.5, cex.sub=0.5, col.axis = "black",las=2)
	  axis(2, at=seq(0,1,0.1), cex.lab=1.2, cex.axis=1.0, cex.main=0.5, cex.sub=0.5,col.axis = "black")
    abline(h=0.5,lty=1,col="red")
  legend(x=mean(xb),y=1.5,c(paste("No. of snps=", n.snp),paste("No. of MC iterations=",n.itr),paste("'a' parameter=",A),paste("'p' parameter=", p.param),paste("No of simulations=",n.sym),paste("No of individuals=",n.obs),paste("Max no. of leaves per tree=", (n.lvs+1)/2),paste("Max no. of trees=", n.trs)),pch=c(1,1,1,1,1),col=c("black","black","black","black","black"),inset = .01,cex=0.4,bty="n", title.adj=0.15,title="Simulation settings")
dev.off()
}
###########################################################
Intnames<- NULL
PostP<-NULL
ex.nam<-NULL
b.s<-NULL
m.int<-NULL
if(!is.null(m.tp.int)&&!is.null(m.fp.int)){m.int=rbind(as.matrix(m.tp.int),as.matrix(m.fp.int))}else{m.int=mat.bind(m.tp.int,m.fp.int)}


#########################################################
tabnew=NULL
tabnew1=NULL
if(!is.null(TP.int)) colnames(TP.int)=c(1:ncol(TP.int))
if(!is.null(FP.int)) colnames(FP.int)=c(1:ncol(FP.int))
tab.int.P=NULL 
tab.int.P<-mat.bind(TP.int,FP.int)
tab.int.pp=NULL
if(!is.null(tab.int.P))
if(nrow(tab.int.P)>0){
tab.int.P=as.matrix(tab.int.P)
tab.int.pp=t(matrix(as.numeric(t(tab.int.P[,1:7])),nrow=7))
}
if(!is.null(tab.int.pp)){
if(is.matrix(tab.int.p)){
if(nrow(tab.int.pp)>0){
ex.nam<-apply(tab.int.pp,1,tree.to.int,n.snp=n.snp)
tabnew=cbind(paste("$",ex.nam,"$"),tab.int.P[,8],round(as.numeric(tab.int.P[,8])/n.sym,digits=5),tab.int.P[,9])
colnames(tabnew)=NULL
tabnew1=NULL
if(is.matrix(tabnew)){
if(nrow(tabnew)>0)
tabnew1=as.matrix(apply(tabnew,1,tex.tab))
}
if(!is.matrix(tabnew)){
if(length(tabnew)>0)
tabnew1=tex.tab(tabnew)
}
write.table(tabnew1, file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}
}

if(!is.matrix(tab.int.p)){
if(length(tab.int.p)>0){
write.table(tab.int.p[1], file = "al.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  m1.int.v<-read.table(file = "al.int.txt")
  fn <- "al.int.txt"
  if (file.exists(fn)) file.remove(fn)
  ex.nam<-apply(m1.int.v,1,tree.to.int,n.snp=n.snp)
tabnew=cbind(paste("$",ex.nam,"$"),tab.int.p[2],round(as.numeric(tab.int.p[2])/n.sym,digits=5),tab.int.p[3])
colnames(tabnew)=NULL
tabnew1=NULL
if(is.matrix(tabnew)){
if(nrow(tabnew)>0)
tabnew1=as.matrix(apply(tabnew,1,tex.tab))
}
if(!is.matrix(tabnew)){
if(length(tabnew)>0)
tabnew1=tex.tab(tabnew)
}

write.table(tabnew1, file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}
}
}

write.table(paste("\\end{tabular}\\end{table}"),file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
############# tex tables with fp and tp for int effects
# l.tp.int=0
# l.fp.int=0
# if(!is.null(TP.int)) l.tp.int=nrow(TP.int)
# if(!is.null(FP.int)) l.fp.int=nrow(FP.int)
# if(l.tp.int>0)
# write.table(tabnew1[1:l.tp.int,], file = paste(OutFile,"_tp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
ila.tp.i=NULL
ila.tp.i=res.tpfp$ila.tp.i
true.names=apply(true.int1,1,tree.to.int,n.snp=n.snp)
if(length(ila.tp.i)>0){
for(rr in 1:length(ila.tp.i)){
if(!is.null(ila.tp.i[[rr]])){
if(is.matrix(ila.tp.i[[rr]])){
if(nrow(ila.tp.i[[rr]])>0){
write.table(paste("\\begin{table}[!h]\\caption{True positives for $",true.names[rr],"$}\\label{}\\centering\\begin{tabular}{|c|c|c|c|}\\hline"), file = paste(OutFile,"_tp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(c(paste("Expression & No.  & Frequency & Posterior probability\\ \\hline")), file = paste(OutFile,"_tp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
foo=ila.tp.i[[rr]]
ila.tp.i[[rr]]=foo[order(as.numeric(foo[,9]),foo[,1],decreasing=TRUE),]
ex.nam<-apply(ila.tp.i[[rr]][,1:7],1,tree.to.int,n.snp=n.snp)
tabnew=cbind(paste("$",ex.nam,"$"),ila.tp.i[[rr]][,8],round(as.numeric(ila.tp.i[[rr]][,8])/n.sym,digits=5),ila.tp.i[[rr]][,9])
colnames(tabnew)=NULL
tabnew1=NULL
if(is.matrix(tabnew)){
if(nrow(tabnew)>0)
tabnew1=as.matrix(apply(tabnew,1,tex.tab))
}
if(!is.matrix(tabnew)){
if(length(tabnew)>0)
tabnew1=tex.tab(tabnew)
}
write.table(tabnew1, file = paste(OutFile,"_tp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"),file = paste(OutFile,"_tp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}
}

}
}
}



# if(!is.matrix(tab.int.p)){
# if(length(tab.int.p)>0){
# write.table(tab.int.p[1], file = "al.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  # m1.int.v<-read.table(file = "al.int.txt")
  # fn <- "al.int.txt"
  # if (file.exists(fn)) file.remove(fn)
  # ex.nam<-apply(m1.int.v,1,tree.to.int,n.snp=n.snp)
# tabnew=cbind(paste("$",ex.nam,"$"),tab.int.p[2],round(as.numeric(tab.int.p[2])/n.sym,digits=5),tab.int.p[3])
# colnames(tabnew)=NULL
# tabnew1=NULL
# if(is.matrix(tabnew)){
# if(nrow(tabnew)>0)
# tabnew1=as.matrix(apply(tabnew,1,tex.tab))
# }
# if(!is.matrix(tabnew)){
# if(length(tabnew)>0)
# tabnew1=tex.tab(tabnew)
# }

# write.table(tabnew1, file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
# }
# }
# }

# ila.tp.i 
# write.table(paste("\\end{tabular}\\end{table}"),file = paste(OutFile,"_tp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
####################################




###########################################
# if(l.tp.int==0){
# if(l.fp.int>0)
# write.table(tabnew1[1:l.fp.int,], file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))}else{
# if(l.fp.int>0)
# write.table(tabnew1[(l.tp.int+1):l.fp.int,], file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
# }

foo=FP.int
FP.int=foo[order(as.numeric(foo[,9]),foo[,1],decreasing=TRUE),]
ex.nam<-apply(FP.int[,1:7],1,tree.to.int,n.snp=n.snp)
tabnew=cbind(paste("$",ex.nam,"$"),FP.int[,8],round(as.numeric(FP.int[,8])/n.sym,digits=5),FP.int[,9])
colnames(tabnew)=NULL
tabnew1=NULL
if(is.matrix(tabnew)){
if(nrow(tabnew)>0)
tabnew1=as.matrix(apply(tabnew,1,tex.tab))
}
if(!is.matrix(tabnew)){
if(length(tabnew)>0)
tabnew1=tex.tab(tabnew)
}
write.table(tabnew1, file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"),file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

write.table(paste("\\begin{table}[!h]\\caption{FDR}\\label{}\\centering\\begin{tabular}{|c|c|c|}\\hline"), file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(c(paste('FDR  & additive FDR & interaction FDR\\ \\hline')), file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
FDR.vec=c(FDR,FDR.a,FDR.i)
FDR.ttab=tex.tab(FDR.vec)
write.table(FDR.ttab, file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"),file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

sorted.add=NULL
if(!is.null(tabnew_add)) sorted.add=tabnew_add[order(as.numeric(tabnew_add[,2]),tabnew_add[,1],decreasing=TRUE),]
sorted.int=NULL
if(!is.null(tabnew)) sorted.int=tabnew[order(as.numeric(tabnew[,2]),tabnew[,1],decreasing=TRUE),]
tabnew1=NULL
if(is.matrix(sorted.add)){
if(nrow(sorted.add)>0)
tabnew1=as.matrix(apply(sorted.add,1,tex.tab))
}
if(!is.matrix(sorted.add)){
if(length(sorted.add)>0)
tabnew1=tex.tab(sorted.add)
}
write.table(paste("\\begin{tabular}{|c|c|c|c|}\\hline"), file = paste(OutFile,"_sorted_tab.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(c(paste('Expression & No.  & Frequency & Posterior probability\\ \\hline')), file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tabnew1, file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}"),file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

tabnew1=NULL
if(is.matrix(sorted.int)){
if(nrow(sorted.int)>0)
tabnew1=as.matrix(apply(sorted.int,1,tex.tab))
}
if(!is.matrix(sorted.int)){
if(length(sorted.int)>0)
tabnew1=tex.tab(sorted.int)
}
write.table(paste("\\begin{tabular}{|c|c|c|c|}\\hline"), file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(c(paste('Expression & No.  & Frequency & Posterior probability\\ \\hline')), file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tabnew1, file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}"),file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

########## power for ILA
pow.add=TP.add
pow.int=TP.int
mean.pow.add=0
mean.pow.int=0
mean.pow.total=0
if(!is.null(pow.add)) mean.pow.add=mean(pow.add[,3])
if(!is.null(pow.int)) mean.pow.int=mean(pow.int[,9])
if(!(is.null(pow.int)&&is.null(pow.int))) mean.pow.total=mean(c(pow.add[,3],pow.int[,9]))

res.tab=NULL
res.tab=rbind(c(FDR,FDR.a,FDR.i),c(mean.pow.total,mean.pow.add,mean.pow.int))

# # # exclude those, which were detected in less than 10% of datasets
# ######### power for ILA
# pow.add.red=TP.add[-which(TP.add[,2]/n.sym<=0.1),]
# pow.int.red=TP.int[-which(TP.int[,8]/n.sym<=0.1),]

# mean.pow.add=mean(pow.add.red[,3])
# mean.pow.int=mean(pow.int.red[,9])
# mean.pow.total=mean(c(pow.add.red[,3],pow.int.red[,9]))

# res.tab.red=NULL
# res.tab.red=rbind(c(FDR,FDR.a,FDR.i),c(mean.pow.total,mean.pow.add,mean.pow.int))


# exclude those, which were detected in less than 15% of datasets
# ######### power for ILA
# pow.add.red=TP.add[-which(TP.add[,2]/n.sym<=0.15),]
# pow.int.red=TP.int[-which(TP.int[,8]/n.sym<=0.15),]

# mean.pow.add=mean(pow.add.red[,3])
# mean.pow.int=mean(pow.int.red[,9])
# mean.pow.total=mean(c(pow.add.red[,3],pow.int.red[,9]))

# res.tab.red=NULL
# res.tab.red=rbind(c(FDR,FDR.a,FDR.i),c(mean.pow.total,mean.pow.add,mean.pow.int))

#FDR.POWER.ILA(true.int,tab.add.p,tab.int.p,max.thres=0.5,OutFile=OutFile)


### power for ILA for true expressions

## adds
adds=true.int1[which(apply(true.int1,1,is.1way)),1]


if(length(adds)>0){
i.pow=NULL
pow.tp.a=NULL
for(s in 1:length(adds))
{
i.pow=c(i.pow,which(tab.add.p[,1]==adds[s]))
if(length(i.pow)>0) pow.tp.a=c(pow.tp.a,tab.add.p[i.pow[s],2]/n.sym)
}

P.tp.a=NULL
P.tp.a.tex=NULL
if(!is.null(pow.tp.a)){
P.tp.a=cbind(namesL[adds],pow.tp.a)
colnames(P.tp.a)=c("True expression","Power")
P.tp.a.tex=apply(P.tp.a,1,tex.tab)


write.table(paste("\\begin{table}[!h]\\caption{Power for ILA}\\label{}\\centering\\begin{tabular}{|c|c|}\\hline"), file = paste(OutFile,"_Power_ILA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tex.tab(c("Expression","Power")),file = paste(OutFile,"_Power_ILA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(P.tp.a.tex,file = paste(OutFile,"_Power_ILA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"), file = paste(OutFile,"_Power_ILA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}
}
# P.tp.a.tex
# P.tp.a

### power for true expressions
IA=which(apply(true.int1,1,is.1way))
if(length(IA)>0){tints=true.int1[-IA,]}else{tints=true.int1}
tints
t.ints=apply(tints,1,intpaste)
tp.ints=tab.int.p[,1]
II=rep(0,length(t.ints))
for(d in 1:length(t.ints))
II[d]=ifelse(!is.na(match(t.ints[d],tp.ints)),match(t.ints[d],tp.ints),0)
 pow.II=NULL
 pow.II=as.numeric(tab.int.p[II,2])/n.sym
P.tp.i=NULL 
P.tp.i.tex=NULL
if(length(pow.II)>0){
if(length(which(II>0))>1){
P.tp.i=cbind(
 as.vector(t(apply(tints[which(II>0),],1,tree.to.int,n.snp=n.snp))),pow.II)
P.tp.i.tex=apply(P.tp.i,1,tex.tab)
}else{P.tp.i=c(tree.to.int(tints,n.snp=n.snp),pow.II)
P.tp.i.tex=tex.tab(P.tp.i)}
}
 
 
write.table(paste("\\begin{table}[!h]\\caption{Power for ILA}\\label{}\\centering\\begin{tabular}{|c|c|}\\hline"), file = paste(OutFile,"_Power_ILA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tex.tab(c("Expression","Power")),file = paste(OutFile,"_Power_ILA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(P.tp.i.tex,file = paste(OutFile,"_Power_ILA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"), file = paste(OutFile,"_Power_ILA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))


E.fp=sum(c(FP.add[,2],FP.int[,8]))/n.sym
FDR1=NULL
FDR1=sum(FP.add[,2])/(sum(FP.add[,2])+sum(TP.add[,2]))

twowaysFP=NULL
where2way=NULL
where2way=apply(FP.int,1,is.2way)
twowaysFP=FP.int[which(where2way),]
twowaysFP
twowaysTP=NULL
int2=NULL
where2way=NULL
where2way=apply(TP.int,1,is.2way)
twowaysTP=TP.int[which(where2way),]
twowaysTP
if(!is.matrix(twowaysFP)) twowaysFP=t(as.matrix(twowaysFP))
if(!is.matrix(twowaysTP)) twowaysTP=t(as.matrix(twowaysTP))

FDR2=NULL
FDR2= sum(twowaysFP[,8])/(sum(twowaysFP[,8])+sum(twowaysTP[,8]))



threewaysFP=NULL
where3way=NULL
where3way=apply(FP.int,1,is.3way)
threewaysFP=FP.int[which(where3way),]
threewaysFP
threewaysTP=NULL
int3=NULL
where3way=NULL
where3way=apply(TP.int,1,is.3way)
threewaysTP=TP.int[which(where3way),]

if(!is.matrix(threewaysFP)) threewaysFP=t(as.matrix(threewaysFP))
if(!is.matrix(threewaysTP)) threewaysTP=t(as.matrix(threewaysTP))
FDR3=NULL
FDR3= sum(threewaysFP[,8])/(sum(threewaysFP[,8])+sum(threewaysTP[,8]))



fourwaysFP=NULL
where4way=NULL
where4way=apply(FP.int,1,is.4way)
fourwaysFP=FP.int[which(where4way),]

fourwaysTP=NULL
int4=NULL
where4way=NULL
where4way=apply(TP.int,1,is.4way)
fourwaysTP=TP.int[which(where4way),]
if(!is.matrix(fourwaysFP)) fourwaysFP=t(as.matrix(fourwaysFP))
if(!is.matrix(fourwaysTP)) fourwaysTP=t(as.matrix(fourwaysTP))

FDR4=NULL
FDR4= sum(fourwaysFP[,8])/(sum(fourwaysFP[,8])+sum(fourwaysTP[,8]))
FDR4



addsFP=NULL
where1way=NULL
where1way=apply(FP.int,1,is.1way)
addsFP=FP.int[which(where1way),]
addsFP
addsTP=NULL
int1=NULL
where1way=NULL
where1way=apply(TP.int,1,is.1way)
addsTP=TP.int[which(where1way),]
addsTP
# FDR1=sum(addsFP[,8])/(sum(addsFP[,8])+sum(addsTP[,8]))
# FDR1
#E.fp=sum(FP.int[,8])/n.sym
Efp1=sum(FP.add[,2])/n.sym
Efp1
Efp2=sum(twowaysFP[,8])/n.sym
Efp2
Efp3=sum(threewaysFP[,8])/n.sym
Efp3
Efp4=sum(fourwaysFP[,8])/n.sym
Efp4
EFp=c(Efp1,Efp2,Efp3,Efp4,E.fp)
efp=cbind(paste(c("$E(FP)_1$","$E(FP)_2$","$E(FP)_3$","$E(FP)_4$","$E(FP)$ Total"),sep=""),EFp)
efp=apply(efp,1,tex.tab)

FDr=cbind(paste(c("FDR1","FDR2","FDR3","FDR4","FDR Total"),sep=""),c(FDR1,FDR2,FDR3,FDR4,FDR))
FDr=apply(FDr,1,tex.tab)

write.table(paste("\\begin{table}[!h]\\caption{False discovery rate :ILA }\\label{}\\centering\\begin{tabular}{|c|c|}\\hline"), file = paste(OutFile,"_Power_ILA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(FDr,file = paste(OutFile,"_Power_ILA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"), file = paste(OutFile,"_Power_ILA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

write.table(paste("\\begin{table}[!h]\\caption{Expected number of false positives :ILA }\\label{}\\centering\\begin{tabular}{|c|c|}\\hline"), file = paste(OutFile,"_Power_ILA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(efp,file = paste(OutFile,"_Power_ILA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"), file = paste(OutFile,"_Power_ILA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}
}




#########################################################################
mat.bind=function(v1,v2)
# function to bind two matrices by rows 
{
if(is.null(v1)&&(!is.null(v2))) return(v2)
if((!is.null(v1))&&(is.null(v2))) return(v1)
if((!is.null(v1))&&(!is.null(v2))) return(rbind(v1,v2))
if((is.null(v1))&&(is.null(v2))) return(NULL)
}
##########################################################################
int.lev.tree=function(vec)
{
Lop=length(which(vec>=1000))
return(7-(length(which(vec==0))+Lop+1))
}
##########################################################################


FDR.POWER.ILA=function(true.int,tab.add.p,tab.int.p,max.thres=max.thres,OutFile=OutFile)
{


write.table("", file = paste(OutFile,"_Power_tab.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

########## ILA
mean.add.pow=NULL
mean.int.pow=NULL
mean.pow=NULL
FDR.a.t.ila=NULL
FDR.i.t.ila=NULL
FDR.t.ila=NULL
Pow.a.ila=NULL
Pow.int.ila=NULL
mat.a=NULL
step.l=max.thres/10
pow.t.a.ila=vector("list",11)
pow.t.i.ila=vector("list",11)
for(hh in 0:10){
res=tp.fp.ILA(true.int,tab.add.p,tab.int.p,thres=hh*step.l)
FDR.a.t.ila=c(FDR.a.t.ila,res$FDR.a.t)
FDR.i.t.ila=c(FDR.i.t.ila,res$FDR.i.t)
FDR.t.ila=c(FDR.t.ila,res$FDR.t)

ia=which(res$TP.add[,2]/nsym>hh*step.l)
if(length(ia)>0){
Pow.a.ila=cbind(res$TP.add.t[ia,1],res$TP.add.t[ia,2]/n.sym,res$TP.add.t[ia,3])
colnames(Pow.a.ila)=c("Index","Power","Posterior")
pow.t.a.ila[[hh+1]]=Pow.a.ila
foo=Pow.a.ila
Pow.a.ila=foo[order(as.numeric(foo[,2]),foo[,1],decreasing=TRUE),]
mat=cbind(Pow.a.ila[,1],matrix(rep(0,6*nrow(Pow.a.ila)),ncol=6))
mat=apply(mat,1,tree.to.int,n.snp=n.snp)
mat.a=cbind(mat,Pow.a.ila[,2:3])
mat.a[,3]=round(as.numeric(mat.a[,3]),digits=5)
p.a.ila=apply(mat.a,1,tex.tab)
}

ii=which(res$TP.int.t[,8]/nsym>hh*step.l)
if(length(ii)>0){
Pow.int.ila=cbind(res$TP.int.t[ii,1:7],res$TP.int.t[ii,8]/n.sym,res$TP.int.t[ii,9])
colnames(Pow.int.ila)=c(c(1:7),"Power","Posterior")
pow.t.i.ila[[hh+1]]=Pow.int.ila
p.i.ila=cbind(apply(Pow.int.ila[,1:7],1,intpaste),Pow.int.ila[,8:9])
foo=p.i.ila
p.i.ila=foo[order(as.numeric(foo[,2]),foo[,1],decreasing=TRUE),]
#p.i.ila[which(as.numeric(p.i.ila[,2])>0.1),]
	write.table(p.i.ila[which(as.numeric(p.i.ila[,2])>0.1),], file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
Pow.int.ila=read.table(file = "M.int.txt")
	fn <- "M.int.txt"
	if (file.exists(fn)) file.remove(fn)
p.i.ila=cbind(apply(Pow.int.ila[,1:7],1,tree.to.int_Logic,n.snp=n.snp),Pow.int.ila[,8:9])
colnames(p.i.ila)=c("Expression","Power","Posterior")

p.i.ila=as.matrix(p.i.ila)
p.i.ila[,3]=round(as.numeric(p.i.ila[,3]),digits=5)

foo=rbind(mat.a,as.matrix(p.i.ila))
p.ila=foo[order(as.numeric(foo[,2]),foo[,1],decreasing=TRUE),]
p.ila=apply(p.ila,1,tex.tab)
}

mean.add.pow=c(mean.add.pow,mean(Pow.a.ila[,2]))
mean.int.pow=c(mean.int.pow,mean(Pow.int.ila[,8]))
mean.pow=c(mean.pow,mean(c(Pow.a.ila[,2],Pow.int.ila[,8])))

## power tables
write.table(paste("Threshold: ",hh*step.l), file = paste(OutFile,"_Power_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))


write.table(paste("\\begin{table}[!h]\\caption{Power for ILA}\\label{}\\centering\\begin{tabular}{|c|c|c|}\\hline"), file = paste(OutFile,"_Power_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tex.tab(c("Expression","Power","Posterior")),file = paste(OutFile,"_Power_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(p.ila,file = paste(OutFile,"_Power_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"), file = paste(OutFile,"_Power_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))


}
################### FDR plot

F.ila=cbind(FDR.a.t.ila,FDR.i.t.ila,FDR.t.ila)
Power.ila=cbind(mean.add.pow,mean.int.pow,mean.pow)

jpeg(file=file.path(OutFile,"_FDR-t-ILA.jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
plot(c(0:10)*step.l,FDR.t.ila,ylim=c(0,0.5),xlab="Threshold",ylab="FDR", main=NULL,type="l",col="red",lty=2,lwd=2)
# points(c(0:10)*step.l,FDR.i.t.eta,type="l",col="red",lty=1,lwd=2)
# points(c(0:10)*step.l,FDR.a.t.ceta,type="l",col="blue",lty=3,lwd=2)
# points(c(0:10)*step.l,FDR.i.t.ceta,type="l",col="blue",lty=5,lwd=2)
points(c(0:10)*step.l,FDR.a.t.ila,type="l",col="blue",lty=3,lwd=2)
points(c(0:10)*step.l,FDR.i.t.ila,type="l",col="black",lty=5,lwd=2)
legend(x=5*step.l,y=0.5,c("ILA FDR total","ILA FDR add", "ILA FDR int"),lty=c(2,1,3,5,5),col=c("red","blue","black"),inset = .01,cex=0.6,bty="n", title.adj=.55)
dev.off()


jpeg(file=file.path(OutFile,"_power-t-ILA.jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
plot(c(0:10)*step.l,mean.add.pow,ylim=c(0,1.5),xlab="Threshold",ylab="Power", main=NULL,type="l",col="red",lty=2,lwd=2)
points(c(0:10)*step.l,mean.int.pow,type="l",col="red",lty=1,lwd=2)
points(c(0:10)*step.l,mean.pow,type="l",col="blue",lty=3,lwd=2)
legend(x=5*step.l,y=1.5,c("power ILA add", "power ILA int","power ILA total"),lty=c(2,1,3),col=c("red","red","blue"),inset = .01,cex=0.6,bty="n", title.adj=.55)
dev.off()



# FDR & power table 
write.table(paste("\\begin{table}[!h]\\caption{False discovery rate}\\label{}\\centering\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}\\hline"), file = paste(OutFile,"_FDR-ila_tab.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tex.tab(c("&Threshold","ILA add","ILA int","ILA Total")),file = paste(OutFile,"_FDR-ila_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(cbind(c(0:10)*step.l,apply(round(F.ila,digits=5),1,tex.tab)),file = paste(OutFile,"_FDR-ila_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"), file = paste(OutFile,"_FDR-ila_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

# 
write.table(paste("\\begin{table}[!h]\\caption{Mean power}\\label{}\\centering\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}\\hline"), file = paste(OutFile,"_power-ila_tab.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tex.tab(c("&Threshold","ILA add","ILA int","ILA Total")),file = paste(OutFile,"_power-ila_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(apply(round(cbind(c(0:10)*step.l,mean.add.pow,mean.int.pow,mean.pow),digits=5),1,tex.tab),file = paste(OutFile,"_power-ila_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"), file = paste(OutFile,"_power-ila_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))



}



FDR.POWER.ETA=function(true.int,tab.add.p,tab.int.p,max.thres=max.thres,OutFile=OutFile){
########## ETA
write.table("", file = paste(OutFile,"_Power_ETA_tab.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

########## ETA
mean.add.pow=NULL
mean.int.pow=NULL
mean.pow=NULL
pow.t.a.eta=vector("list",11)
pow.t.i.eta=vector("list",11)
Pow.a.eta=NULL
Pow.int.eta=NULL
FDR.a.t.eta=NULL
FDR.i.t.eta=NULL
TP.add.t.eta=NULL
TP.int.t.eta=NULL
FP.add.t.eta=NULL
FP.int.t.eta=NULL
FDR.t.eta=NULL
step.l=max.thres/10
for(hh in 0:10){
res=tp.fp.ETA(true.int,tab.add.p,tab.int.p,thres=hh*step.l)
FDR.a.t.eta=c(FDR.a.t.eta,res$FDR.a.t)
FDR.i.t.eta=c(FDR.i.t.eta,res$FDR.i.t)
FDR.t.eta=c(FDR.t.eta,res$FDR.t)

ia=which(res$TP.add[,2]/nsym>hh*step.l)
Pow.a.eta=NULL
if(length(ia)>0){
Pow.a.eta=cbind(res$TP.add.t[ia,1],res$TP.add.t[ia,2]/n.sym,res$TP.add.t[ia,3])
colnames(Pow.a.eta)=c("Index","Power","Posterior")
pow.t.a.eta[[hh+1]]=Pow.a.eta
}

ii=which(res$TP.int.t[,8]/nsym>hh*step.l)
Pow.int.eta=NULL
if(length(ii)>0){
Pow.int.eta=cbind(res$TP.int.t[ii,1:7],res$TP.int.t[ii,8][,1]/n.sym,res$TP.int.t[ii,9])
colnames(Pow.int.eta)=c(c(1:7),"Power","Posterior")
pow.t.i.eta[[hh+1]]=Pow.int.eta
}


mean.add.pow=c(mean.add.pow,ifelse(!is.null(pow.t.a.eta[[hh+1]]),mean(pow.t.a.eta[[hh+1]][,2]),0))


mean.int.pow=c(mean.int.pow,ifelse(!is.null(pow.t.i.eta[[hh+1]]),mean(pow.t.i.eta[[hh+1]][,8]),0))

# mean.pow=c(mean.pow,mean(as.numeric(p.eta[,2])))

# ifelse(!is.null(rbind(pow.t.a.eta[[hh+1]],pow.t.i.eta[[hh+1]])),mean(),0)



mat.a=NULL
p.a.eta=NULL
if(!is.null(Pow.a.eta)){
foo=Pow.a.eta
Pow.a.eta=foo[order(as.numeric(foo[,2]),foo[,1],decreasing=TRUE),]
mat=cbind(Pow.a.eta[,1],matrix(rep(0,6*nrow(Pow.a.eta)),ncol=6))
mat=apply(mat,1,tree.to.int,n.snp=n.snp)
mat.a=cbind(mat,Pow.a.eta[,2:3])
mat.a[,3]=round(as.numeric(mat.a[,3]),digits=5)
p.a.eta=apply(cbind(mat,Pow.a.eta[,2:3]),1,tex.tab)
# mean.add.pow=c(mean.add.pow,mean(Pow.a.eta[,2]))
}


### pow int
p.i.eta=NULL
if(!is.null(Pow.int.eta)){
p.i.eta=cbind(apply(Pow.int.eta[,1:7],1,intpaste),Pow.int.eta[,8:9])
foo=p.i.eta
p.i.eta=foo[order(as.numeric(foo[,2]),foo[,1],decreasing=TRUE),]

	write.table(p.i.eta, file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
Pow.int.eta=read.table(file = "M.int.txt")
	fn <- "M.int.txt"
	if (file.exists(fn)) file.remove(fn)
p.i.eta=cbind(apply(Pow.int.eta[,1:7],1,tree.to.int_Logic,n.snp=n.snp),Pow.int.eta[,8:9])
colnames(p.i.eta)=c("Expression","Power","Posterior")
# p.i.eta=apply(p.i.eta,1,tex.tab)
p.i.eta=as.matrix(p.i.eta)
p.i.eta[,3]=round(as.numeric(p.i.eta[,3]),digits=5)
#mean.int.pow=c(mean.int.pow,mean(Pow.int.eta[,8]))
p.i.eta=as.matrix(p.i.eta)
}


foo=rbind(mat.a,p.i.eta)
p.eta=NULL
if(!is.null(foo)){
p.eta=foo[order(as.numeric(foo[,2]),foo[,1],decreasing=TRUE),]
}

mean.pow[hh+1]=0
if(!is.null(p.eta)){
if(is.matrix(p.eta)){
mean.pow[hh+1]=mean(as.numeric(p.eta[,2]))
}else{
mean.pow[hh+1]=mean(as.numeric(t(as.matrix(p.eta))[,2]))}

p.eta=apply(t(as.matrix(p.eta)),1,tex.tab)


## power tables
write.table(paste("Threshold: ",hh*step.l), file = paste(OutFile,"_Power-ETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))


write.table(paste("\\begin{table}[!h]\\caption{Power for ETA}\\label{}\\centering\\begin{tabular}{|c|c|c|}\\hline"), file = paste(OutFile,"_Power_ETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tex.tab(c("Expression","Power","Posterior")),file = paste(OutFile,"_Power_ETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(p.eta,file = paste(OutFile,"_Power_ETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"), file = paste(OutFile,"_Power_ETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}

}
################### FDR plot

#F.eta=cbind(FDR.a.t.eta,FDR.i.t.eta,FDR.t.eta)
F.eta=cbind(c(0:10)*step.l,FDR.a.t.eta,FDR.i.t.eta,FDR.t.eta)
colnames(F.eta)=c("Frequency threshold", "FDR add", "FDR int", "FDR total")

Power.eta=cbind(mean.add.pow,mean.int.pow,mean.pow)

jpeg(file=file.path(OutFile,"_FDR-t-ETA.jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
plot(c(0:10)*step.l,FDR.t.eta,ylim=c(0,0.5),xlab="Threshold",ylab="FDR", main=NULL,type="l",col="red",lty=2,lwd=2)
# points(c(0:10)*step.l,FDR.i.t.eta,type="l",col="red",lty=1,lwd=2)
# points(c(0:10)*step.l,FDR.a.t.ceta,type="l",col="blue",lty=3,lwd=2)
# points(c(0:10)*step.l,FDR.i.t.ceta,type="l",col="blue",lty=5,lwd=2)
points(c(0:10)*step.l,FDR.a.t.eta,type="l",col="blue",lty=3,lwd=2)
points(c(0:10)*step.l,FDR.i.t.eta,type="l",col="black",lty=5,lwd=2)
legend(x=5*step.l,y=0.5,c("ETA FDR total","ETA FDR add", "ETA FDR int"),lty=c(2,1,3,5,5),col=c("red","blue","black"),inset = .01,cex=0.6,bty="n", title.adj=.55)
dev.off()


jpeg(file=file.path(OutFile,"_power-t-ETA.jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
plot(c(0:10)*step.l,mean.add.pow,ylim=c(0,1.5),xlab="Threshold",ylab="Power", main=NULL,type="l",col="red",lty=2,lwd=2)
points(c(0:10)*step.l,mean.int.pow,type="l",col="red",lty=1,lwd=2)
points(c(0:10)*step.l,mean.pow,type="l",col="blue",lty=3,lwd=2)
legend(x=5*step.l,y=1.5,c("power ETA add", "power ETA int","power ETA total"),lty=c(2,1,3),col=c("red","red","blue"),inset = .01,cex=0.6,bty="n", title.adj=.55)
dev.off()



# FDR & power table 
write.table(paste("\\begin{table}[!h]\\caption{False discovery rate}\\label{}\\centering\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}\\hline"), file = paste(OutFile,"_FDR-eta_tab.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tex.tab(c("Threshold","ETA add","ETA int","ETA Total")),file = paste(OutFile,"_FDR-eta_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(apply(cbind(c(0:10)*step.l,round(F.eta,digits=5)),1,tex.tab),file = paste(OutFile,"_FDR-eta_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"), file = paste(OutFile,"_FDR-eta_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

# 
write.table(paste("\\begin{table}[!h]\\caption{Mean power}\\label{}\\centering\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}\\hline"), file = paste(OutFile,"_power-eta_tab.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tex.tab(c("&Threshold","ETA add","ETA int","ETA Total")),file = paste(OutFile,"_power-eta_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(apply(round(cbind(c(0:10)*step.l,mean.add.pow,mean.int.pow,mean.pow),digits=5),1,tex.tab),file = paste(OutFile,"_power-eta_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"), file = paste(OutFile,"_power-eta_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}



summarize.ETA=function(OutFile,A,p.param,n.sym,true.int,object,n.snp,n.obs,n.itr,boxplots.all=FALSE)
## function to summarize results of simulations for ETA approach
## calculates FDR and power, makes plots if applicable
## files containing tables in tex-syntax are created
{
fs=file.info(paste(OutFile,".txt",sep=""))$size
if(fs>0){			
		M1.1<-NULL
		M1.1<-read.table(file=paste(OutFile,".txt",sep=""), header = FALSE, sep = "\t", quote = "\"'", dec = ".",as.is = TRUE,na.strings = "NA", colClasses = NA, nrows = -1,
				   skip = 0, check.names = TRUE, strip.white = FALSE, blank.lines.skip = TRUE,
				   comment.char = "#",
				   allowEscapes = FALSE, flush = FALSE,
				   stringsAsFactors = FALSE,
				   fileEncoding = "", encoding = "unknown")
add.eff<-NULL
int2<-NULL
int3<-NULL
str2 <- gsub(' {2,}',' ',M1.1[,1])
int.level=unlist(lapply(strsplit(str2,' '),length))
add.eff<-M1.1[which(int.level==1),]
int=M1.1[which(int.level>1),]
res.post.add=NULL
if(nrow(add.eff)>0){
Lv=count.rows(add.eff[,1])
L.p=Lv[,2]
res.post.add=vector("list",length(L.p))
for(i in 1:length(L.p)){
p=NULL
if(length(which(Lv[,3]==i))>0) p=add.eff[L.p[which(Lv[,3]==i)][[1]],2]
res.post.add[[i]]=list(v=as.numeric(Lv[i,3]),n=Lv[i,1],p=p)
}
}
#####
write.table(int, file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
Alli=read.table(file = "M.int.txt")
fn <- "M.int.txt"
if (file.exists(fn)) file.remove(fn)
#####
Alli1=t(apply(Alli[,1:7],1,tree.merge))
Alli1.vec=apply(Alli1[,1:7],1,intpaste)
int.1=cbind(Alli1.vec,int[,2])
int=int.1
#####
res.post.int=NULL
Lvi=NULL
if(nrow(int)>0){
Lvi=count.rows(int[,1])
res.post.int=vector("list",nrow(Lvi))
res.post.int=apply(Lvi,1,do.list,int=int)
}
###############################################
namesL1<-NULL
ind1<-2*c(1:n.snp)-1
ind2<-2*c(1:n.snp)
i=c(1:n.snp)
namesL1[ind1]<-paste("D[",i,"]",sep="")
namesL1[ind2]<-paste("R[",i,"]",sep="")
################ BOXPLOTS
mint.v<-NULL
Ll=NULL
if(!is.null(true.int)){
write.table(true.int, file= "trint.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
mint.v<-read.table(file = "trint.txt")
fn <- "trint.txt"
if (file.exists(fn)) file.remove(fn)
va=t(apply(mint.v,1,variables.from.model))
vva=sort(as.vector(unlist(va)))
TNa<-ifelse(vva<0,paste(namesL1[abs(vva)],"^C",sep=""),paste(namesL1[as.numeric(abs(vva))],sep=""))
ex.nam<-TNa

if(!is.null(res.post.add)){
Ll=res.post.add[vva]
if(boxplots.all){
for(ss in 1:length(vva)){
jpeg(file=file.path(OutFile,"_boxplot_add",ss,".jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
allps=Ll[[ss]]$p
if(length(Ll[[ss]]$p)<n.sym) allps=c(Ll[[ss]]$p,rep(0,(n.sym-length(Ll[[ss]]$p))))
boxplot(allps,
        col = "yellow",
        main = paste("Boxplot of posterior probability for",ex.nam[ss]),
        xlab = " ",
        ylab = "Posterior probability",
        ylim = c(0, 2), yaxs = "i")
		axis(1, c(1), paste(ex.nam[ss],": Number of appearances :",Ll[[ss]]$n,", Number of post.pr >=0.5 :",length(which(allps>=0.5)),sep=""), cex.lab=1.2, cex.axis=0.4, cex.main=0.5, cex.sub=0.5, col.axis = "black",las=1)
dev.off()
}
}
}
}
####### boxplots-altogether
ps=NULL
res=NULL
for(ss in 1:length(vva)){
ps=Ll[[ss]]$p
if(length(Ll[[ss]]$p)<n.sym) ps=c(Ll[[ss]]$p,rep(0,(n.sym-length(Ll[[ss]]$p))))
res=cbind(res,ps)
}
colnames(res)=NULL
jpeg(file=file.path(OutFile,"_boxplot_add_all.jpeg",fsep=""),height=12,width=24,unit="cm",res=300)
boxplot(res[,1:ncol(res)], 
        names=paste(ex.nam[1:ncol(res)],sep=""),
        col = "yellow",
        main = paste("Boxplots of posterior probabilities"),
        xlab = " ",
        ylab = "Posterior probability",
        ylim = c(0, 1.2), yaxs = "i",las=2)
dev.off()
Ll=NULL
true.int1=NULL

###############################################################
##### indices of models containing true interactions
# ################################
if(!is.null(true.int)){
 write.table(true.int, file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  true.int1=read.table(file = "M.int.txt")
   fn <- "M.int.txt"
 if (file.exists(fn)) file.remove(fn)
 }
  Alli=NULL
  ind.tr.i=NULL
   if(!is.null(res.post.int)){
   write.table(Lvi[,3], file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  Alli=read.table(file = "M.int.txt")
   fn <- "M.int.txt"
 if (file.exists(fn)) file.remove(fn)
 }
 ind.tr.i=NULL
  if(!is.null(true.int1)){
  if(length(unique(unlist(lapply(apply(true.int1,1,variables.from.model),length))))==1){
  v=t(apply(true.int1,1,variables.from.model))
  if(nrow(v)>1){
  NO.Int=NULL
  NO.Int=nrow(true.int1)
  ind.tr.i=vector("list",nrow(true.int1))
  for(ij in 1:nrow(true.int1)){
  ind.tr.i[[ij]]=0
  sumi=which(rowSums(Alli)==sum(true.int1[ij,]))
  m1=Alli[which(rowSums(Alli)==sum(true.int1[ij,])),]
  kt=NULL
  for(rr in 1:nrow(m1)) kt=c(kt,(sum(v[ij,]%in%m1[rr,])==length(v[ij,])))
  if(length(which(kt))>0) ind.tr.i[[ij]]=sumi[which(kt)]
  }}
  }else{v=apply(true.int1,1,variables.from.model)# a list
		NO.Int=NULL
		NO.Int=length(which(unlist(lapply(v,length))>1))
		ind.tr.i=vector("list",NO.Int)
		for(ij in 1:NO.Int){
		  ind.tr.i[ij]=0
		  sumi=which(as.vector(rowSums(Alli)==sum(true.int1[ij,])))
		  m1=Alli[which(rowSums(Alli)==sum(true.int1[ij,])),]
		  kt=NULL
		  for(rr in 1:nrow(m1)) kt=c(kt,(sum(v[[ij]]%in%m1[rr,])==length(v[[ij]])))
		  if(length(which(kt))>0) ind.tr.i[[ij]]=sumi[which(kt)]
		  }
	  }
  }
##### boxplots for interactions : individually
  m1.int.v<-NULL
  ex.nam.int<-NULL
  if(!is.null(true.int)){
  write.table(true.int, file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  m1.int.v<-read.table(file = "M.int.txt")
  fn <- "M.int.txt"
  if (file.exists(fn)) file.remove(fn)
  ex.nam.int<-apply(m1.int.v,1,tree.to.int_Logic,n.snp=n.snp)
  }
if(boxplots.all){ 
if(!is.null(ind.tr.i)){ 
for(si in 1:NO.Int){ 
if(ind.tr.i[[si]]>0){
jpeg(file=file.path(OutFile,"_boxplot_int",si,".jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
allpints=res.post.int[[ind.tr.i[[si]]]]$p
if(length(res.post.int[[ind.tr.i[[si]]]]$p)<n.sym) allpints=c(res.post.int[[ind.tr.i[[si]]]]$p,rep(0,(n.sym-length(res.post.int[[ind.tr.i[[si]]]]$p))))
boxplot(allpints,
		col = "yellow",
        ylab = "Posterior probability",
        ylim = c(0, 1.5), yaxs = "i")
		axis(1, c(1),ex.nam.int[si], cex.lab=1.2, cex.axis=0.4, cex.main=0.5, cex.sub=0.5, col.axis = "black",las=2)
dev.off()
}
}
}
}
# ### boxplots -all int together
  # resint=NULL
  # for(si in 1:NO.Int){ 
  # pints=NULL
  # if(ind.tr.i[[si]]>0){
  # pints=res.post.int[[ind.tr.i[[si]]]]$p
  # if(length(res.post.int[[ind.tr.i[[si]]]]$p)<n.sym) pints=c(pmin(1,res.post.int[[ind.tr.i[[si]]]]$p),rep(0,(n.sym-length(res.post.int[[ind.tr.i[[si]]]]$p))))
  # resint=cbind(resint,pints)
  # }
  # colnames(resint)=NULL
  # }
  
# if(length(which(unlist(ind.tr.i)>0))>0){ 
# jpeg(file=file.path(OutFile,"_boxplot_int_all.jpeg",fsep=""),height=24,width=24,unit="cm",res=300)
# par(mar=c(12,8,2,2)+0.3,mgp=c(3,1,0))
# boxplot(resint[,1:ncol(resint)],
		# names=paste(ex.nam.int[which(ind.tr.i>0)],sep=""),
		# col = "yellow",
		# main = paste("Boxplots of posterior probabilities"),
        # ylab = "Posterior probability",
        # ylim = c(0, 1.1), yaxs = "i", las=2)
# dev.off()
# }
################### BOXPLOT OF TREE SIZES
i=NULL
s=NULL
ex.nam=NULL
ex.nam.int=NULL
m1.int.v=NULL
##################
tab.add<-NULL
if(length(res.post.add)>0){
tab.add<-matrix(nrow=length(res.post.add),ncol=3)
for(i in 1:length(res.post.add)){ 
tab.add[i,]<-c(res.post.add[[i]]$v,res.post.add[[i]]$n,sum(as.numeric(res.post.add[[i]]$p)))
}
}

tab.int<-NULL
if(!is.null(res.post.int)){
tab.int<-matrix(nrow=length(res.post.int),ncol=3)
for(i in 1:length(res.post.int)){
tab.int[i,]<-c(res.post.int[[i]]$v$x,res.post.int[[i]]$n$counts,sum(as.numeric(res.post.int[[i]]$p)))
}
}

 tab.add[,3]<-as.numeric(tab.add[,3])/n.sym
 tab.int[,3]<-as.numeric(tab.int[,3])/n.sym 
 if(!is.null(tab.add)) 
 if(length(which(as.numeric(tab.add[,3])>=0))>0)
	{tab.add=tab.add[which(as.numeric(tab.add[,3])>=0),]
	}else{tab.add=NULL}
 if(!is.null(tab.int)) 
 if(length(which(as.numeric(tab.int[,3])>=0))>0){
 tab.int=tab.int[which(as.numeric(tab.int[,3])>=0),]}else{tab.int=NULL}
 
 ind=NULL
 No.a=NULL
for(jj in 1:length(res.post.add)){ 
  if(length(which(res.post.add[[jj]]$p>=0.5))>0){
  ind=c(ind,jj)
  No.a=c(No.a,length(which(res.post.add[[jj]]$p>=0.5)))
  }
 }
 tab.add.p=cbind(tab.add[ind,1],No.a,tab.add[ind,3])
 N.A=NULL
 N.A=sum(No.a)
 
 ind=NULL
 No.i=NULL
if(!is.null(res.post.int)){
for(jj in 1:length(res.post.int)){ 
  if(length(which(res.post.int[[jj]]$p>=0.5))>0){
  ind=c(ind,jj)# which interaction had at lest once a posterior larger than 0.5
  No.i=c(No.i,length(which(res.post.int[[jj]]$p>=0.5)))
  }
 }
 tab.int.p=cbind(tab.int[ind,1],No.i,tab.int[ind,3])
colnames(tab.int.p)=NULL
 tab.add.red=tab.add.p
 tab.int.red=tab.int.p
 Li=NULL
 if(length(tab.int.p)>0){
 write.table(tab.int.p[,1], file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
 Li=read.table(file = paste("M.int.txt",sep=""))
 fn <- "M.int.txt"
 if (file.exists(fn)) file.remove(fn)
 }
 i.lev=NULL
if(!is.null(Li)) i.lev=apply(Li,1,int.lev)
N.2=0
N.3=0
N.4=0
 if(length(which(i.lev==2))>0) N.2=sum(No.i[which(i.lev==2)])
  if(length(which(i.lev==3))>0) N.3=sum(No.i[which(i.lev==3)])
  if(length(which(i.lev==4))>0) N.4=sum(No.i[which(i.lev==4)])
}  
dist.lvs=c(rep(1,N.A),rep(2,N.2),rep(3,N.3),rep(4,N.4))

jpeg(file=file.path(OutFile,"_BoxplotTreeSizes.jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
boxplot(dist.lvs,
        col = "yellow")
dev.off()
########### HISTOGRAM OF TREE SIZES
jpeg(file=file.path(OutFile,"_HistTreeSizes.jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
hist(dist.lvs,col = "lightblue", border = "pink",main="Histogram of selected trees size",xlab="Tree size")
dev.off()

######### FALSE POSITIVES AND TRUE POSITIVES
res.tpfp<-tp.fp.ETA(true.int,tab.add.p,tab.int.p,thres=0.85)
tp.fp.ETA.plot(true.int,tab.add.p,tab.int.p,max.thres=0.85,OutFile=OutFile)
TP.add=res.tpfp$TP.add
TP.int=res.tpfp$TP.int
FP.add=res.tpfp$FP.add
FP.int=res.tpfp$FP.int
# if(!is.null(TP.int)) CRtp=count.rows(as.matrix(TP.int[,1:7]))
# if(!is.null(FP.int)) CRfp=count.rows(as.matrix(FP.int[,1:7]))
# if(!is.null(TP.int)){
	# whichrows=count.rows(TP.int[,1:7])[2]
	# nwr=nrow(whichrows)
	# cmdp<-paste("ptp <- sum(as.numeric(as.matrix(TP.int)[,9])[whichrows[[1]]$'",1:nwr,"'])",sep="")
	# cmdfr<-paste("frtp <- sum(as.numeric(as.matrix(TP.int)[,8])[whichrows[[1]]$'",1:nwr,"'])",sep="")
	# ptp=NULL
	# frtp=NULL
	# for(nn in 1:nwr){
	# ptp=c(ptp,eval(parse(text=cmdp[nn])))
	# frtp=c(frtp,eval(parse(text=cmdfr[nn])))
	# }
	# TP.int<-cbind(CRtp[3:9],frtp,ptp)
	# }
# whichrows=NULL
# whichrows=count.rows(as.matrix(FP.int[,1:7]))[2]
# nwr=nrow(whichrows)
# cmdp<-paste("ptp <- sum(as.numeric(as.matrix(FP.int)[,9])[whichrows[[1]]$'",1:nwr,"'])",sep="")
# cmdfr<-paste("frtp <- sum(as.numeric(as.matrix(FP.int)[,8])[whichrows[[1]]$'",1:nwr,"'])",sep="")
# pfp=NULL
# frfp=NULL
# for(nn in 1:nwr){
# pfp=c(pfp,eval(parse(text=cmdp[nn])))
# frfp=c(frfp,eval(parse(text=cmdfr[nn])))
# }
# FP.int<-cbind(CRfp[3:9],frfp,pfp)

######### FALSE DISCOVERY RATE
FDR=res.tpfp$FDR
FDR.a=res.tpfp$FDR.a
FDR.i=res.tpfp$FDR.i
###############################################
#######		PLOTS AND TABLES WITH SUMMARIES
###############################################
namesL<-NULL
ind1<-2*c(1:n.snp)-1
ind2<-2*c(1:n.snp)
i=c(1:n.snp)
namesL[ind1]<-paste("D_{",i,"}",sep="")
namesL[ind2]<-paste("R_{",i,"}",sep="")
#####################################################
######## table with settings to the file 
#####################################################
set.tab=NULL
set.tab<-rbind(c(paste('True model:'),paste("$",object$snpmodel,"$")),
c(paste("No. of snps"), n.snp),c(paste("No. of MC iterations"),n.itr),c(paste("'a' parameter"),A),c(paste("'p' parameter"), p.param),c(paste("No of simulations"),n.sym),c(paste("No of individuals"),n.obs),c(paste("Max no. of leaves per tree"), (n.lvs+1)/2),c(paste("Max no. of trees"), n.trs))
newset.tab=as.matrix(apply(set.tab,1,tex.tab))
write.table(paste("\\begin{table}[!h] \\caption{Settings }\\label{} \\centering \\begin{tabular}{p{4cm}|p{8cm}}\\hline"), file = paste(OutFile,"_tab.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(newset.tab, file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"),file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

write.table(paste("\\begin{table}[!h]\\caption{Selected expressions}\\label{}\\centering\\begin{tabular}{|c|c|c|c|}\\hline"), file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
tp.add=TP.add[,1]
fp.add=FP.add[,1]
library(plotrix)
Intnames<- c(tp.add,fp.add)
PostP<-NULL
if(!is.null(tab.add.p)){
if(is.matrix(tab.add.p)) PostP<-c(as.numeric(tab.add.p[,3]))
if(!is.matrix(tab.add.p)) PostP<-c(as.numeric(tab.add.p[3]))
}
T1<- NULL
T1<-tab.add.p[,1]#c(tp.add,fp.add)
T1<-as.numeric(T1)
TNa<-ifelse(T1<0,paste(namesL[abs(T1)],"^C",sep=""),paste(namesL[as.numeric(abs(T1))],sep=""))
ex.nam<-TNa
tabnew_add=NULL
##################################################################
write.table(c(paste('Expression & No.  & Frequency & Posterior probability')), file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
if(length(ex.nam)>0){
 if(is.matrix(tab.add.p)) tabnew_add=cbind(paste("$",ex.nam,"$"),tab.add.p[,2],round(as.numeric(tab.add.p[,2])/n.sym,digits=5),tab.add.p[,3])
 if(!is.matrix(tab.add.p)) tabnew_add=cbind(paste("$",ex.nam,"$"),tab.add.p[2],round(as.numeric(tab.add.p[2])/n.sym,digits=5),tab.add.p[3])
colnames(tabnew_add)=NULL
tabnew_add1=NULL
if(is.matrix(tabnew_add)){
if(nrow(tabnew_add)>0)
tabnew_add1=as.matrix(apply(tabnew_add,1,tex.tab))
}
if(!is.matrix(tabnew_add)){
if(length(tabnew_add)>0)
tabnew_add1=tex.tab(tabnew_add)
}
write.table(tabnew_add1, file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}
############# tex tables with fp and tp for add effects
l.tp.add=0
l.fp.add=0
if(!is.null(TP.add)) l.tp.add=nrow(TP.add)
if(!is.null(FP.add)) l.fp.add=nrow(FP.add)
write.table(paste("\\begin{table}[!h]\\caption{True positives}\\label{}\\centering\\begin{tabular}{|c|c|c|c|}\\hline"), file = paste(OutFile,"_tp_tab.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(c(paste("Expression & No.  & Frequency & Posterior probability\\ \\hline")), file = paste(OutFile,"_tp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
if(l.tp.add>0)
write.table(tabnew_add1[1:l.tp.add,], file = paste(OutFile,"_tp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))


write.table(paste("\\begin{table}[!h]\\caption{False positives}\\label{}\\centering\\begin{tabular}{|c|c|c|c|}\\hline"), file = paste(OutFile,"_fp_tab.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(c(paste('Expression & No.  & Frequency & Posterior probability\\ \\hline')), file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
if(l.tp.add==0){
if(l.fp.add>0)
write.table(tabnew_add1[1:l.fp.add,], file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))}else{
if(l.fp.add>0)
write.table(tabnew_add1[(l.tp.add+1):l.fp.add,], file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}

##################################################################
m.tp.int=NULL
m.fp.int=NULL
if(!is.null(TP.int)) if(nrow(TP.int)>0) m.tp.int=TP.int[,1:7]
if(!is.null(FP.int)) if(nrow(FP.int)>0) m.fp.int=FP.int[,1:7]
# m.int<-rbind(m.tp.int,m.fp.int)
m.int<-NULL
if(!is.null(m.tp.int)&&!is.null(m.fp.int)){m.int=rbind(as.matrix(m.tp.int),as.matrix(m.fp.int))}else{m.int=mat.bind(m.tp.int,m.fp.int)}



if(!is.null(m.int)){ 
  PP<-NULL
  colnames(FP.int)=NULL
	colnames(TP.int)=NULL
	tab.int.p1=NULL
	if(!is.null(FP.int)&&!is.null(TP.int)){
					tab.int.p1=rbind(as.matrix(FP.int),as.matrix(TP.int))
										}else{
										if(!is.null(FP.int)) tab.int.p1=as.matrix(FP.int)
										if(!is.null(TP.int)) tab.int.p1=as.matrix(TP.int)
										}
	tab.int.p2=NULL									
	if(!is.null(tab.int.p1)) tab.int.p2=cbind(apply(tab.int.p1[,1:7],1,intpaste),tab.int.p1[,8:9])
	tab.int.p<-tab.int.p2
  if(!is.null(nrow(tab.int.p))){ 
  for(i in 1:nrow(tab.int.p)) PP[i]<-as.double(as.numeric(tab.int.p[i,3]))
  PostP<-c(PostP,PP)
  m1.int.v<-NULL
  write.table(tab.int.p[,1], file = "al.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  }
  if(is.null(nrow(tab.int.p))){
  PP<-as.double(as.numeric(tab.int.p[3]))
  PostP<-c(PostP,PP)
  m1.int.v<-NULL
  write.table(tab.int.p[1], file = "al.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  }
  m1.int.v<-read.table(file = "al.int.txt")
  fn <- "al.int.txt"
  if (file.exists(fn)) file.remove(fn)
  ex.nam<-c(ex.nam,apply(m1.int.v,1,tree.to.int_Logic,n.snp=n.snp))
}

b.s<-NULL
b.s<-cbind(ex.nam,PostP)
if(nrow(b.s)>0){
jpeg(file=file.path(OutFile,"_Posterior_probETA.jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
xb<-c(1:nrow(b.s))
plot(x=xb,y=PostP,lty=1,lwd=2,type="h",axes=FALSE,xlab=NA,ylab="Posterior probability",main=paste("Posterior probability for logic expressions.\nTrue model:",object$snpmodel),ylim=c(0.0,1.5),cex.lab=1.0, cex.axis=0.5, cex.main=0.4, cex.sub=0.5,col="blue")
    axis(1, xb, ex.nam, cex.lab=1.2, cex.axis=0.4, cex.main=0.5, cex.sub=0.5, col.axis = "black",las=2)
	  axis(2, at=seq(0,1,0.1), cex.lab=1.2, cex.axis=1.0, cex.main=0.5, cex.sub=0.5,col.axis = "black")
    abline(h=0.5,lty=1,col="red")
  legend(x=mean(xb),y=1.5,c(paste("No. of snps=", n.snp),paste("No. of MC iterations=",n.itr),paste("'a' parameter=",A),paste("'p' parameter=", p.param),paste("No of simulations=",n.sym),paste("No of individuals=",n.obs),paste("Max no. of leaves per tree=", (n.lvs+1)/2),paste("Max no. of trees=", n.trs)),
pch=c(1,1,1,1,1),col=c("black","black","black","black","black"),inset = .01,cex=0.4,bty="n", title.adj=0.15,title="Simulation settings")
dev.off()
}
#######################
Intnames<- NULL
PostP<-NULL
ex.nam<-NULL
b.s<-NULL
# m.int<-NULL
# m.int=rbind(as.matrix(m.tp.int),as.matrix(m.fp.int))
m.int<-NULL
if(!is.null(m.tp.int)&&!is.null(m.fp.int)){m.int=rbind(as.matrix(m.tp.int),as.matrix(m.fp.int))}else{m.int=mat.bind(m.tp.int,m.fp.int)}



#########################################################
tabnew=NULL
tabnew1=NULL
if(!is.null(TP.int)) colnames(TP.int)=c(1:ncol(TP.int))
if(!is.null(FP.int)) colnames(FP.int)=c(1:ncol(FP.int))
tab.int.P=NULL 
tab.int.P<-mat.bind(TP.int,FP.int)
tab.int.pp=NULL
if(!is.null(tab.int.P))
if(nrow(tab.int.P)>0){
tab.int.P=as.matrix(tab.int.P)
tab.int.pp=t(matrix(as.numeric(t(tab.int.P[,1:7])),nrow=7))
}

if(!is.null(tab.int.pp)){
if(is.matrix(tab.int.p)){
if(nrow(tab.int.pp)>0){
  ex.nam<-apply(tab.int.pp,1,tree.to.int_Logic,n.snp=n.snp)
tabnew=cbind(paste("$",ex.nam,"$"),tab.int.P[,8],round(as.numeric(tab.int.P[,8])/n.sym,digits=5),tab.int.P[,9])
colnames(tabnew)=NULL
tabnew1=NULL
if(is.matrix(tabnew)){
if(nrow(tabnew)>0)
tabnew1=as.matrix(apply(tabnew,1,tex.tab))
}
if(!is.matrix(tabnew)){
if(length(tabnew)>0)
tabnew1=tex.tab(tabnew)
}
write.table(tabnew1, file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}
}


if(!is.matrix(tab.int.p)){
if(length(tab.int.p)>0){
write.table(tab.int.p[1], file = "al.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  m1.int.v<-read.table(file = "al.int.txt")
  fn <- "al.int.txt"
  if (file.exists(fn)) file.remove(fn)
  ex.nam<-apply(m1.int.v,1,tree.to.int,n.snp=n.snp)
tabnew=cbind(paste("$",ex.nam,"$"),tab.int.p[2],round(as.numeric(tab.int.p[2])/n.sym,digits=5),tab.int.p[3])
colnames(tabnew)=NULL
tabnew1=NULL
if(is.matrix(tabnew)){
if(nrow(tabnew)>0)
tabnew1=as.matrix(apply(tabnew,1,tex.tab))
}
if(!is.matrix(tabnew)){
if(length(tabnew)>0)
tabnew1=tex.tab(tabnew)
}

write.table(tabnew1, file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}
}
}

write.table(paste("\\end{tabular}\\end{table}"),file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

############# tex tables with fp and tp for int effects
l.tp.int=0
l.fp.int=0
if(!is.null(TP.int)) l.tp.int=nrow(TP.int)
if(!is.null(FP.int)) l.fp.int=nrow(FP.int)
if(l.tp.int>0)
write.table(tabnew1[1:l.tp.int,], file = paste(OutFile,"_tp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

write.table(paste("\\end{tabular}\\end{table}"),file = paste(OutFile,"_tp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

if(l.tp.int==0){
if(l.fp.int>0)
write.table(tabnew1[1:l.fp.int,], file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))}else{
if(l.fp.int>0)
write.table(tabnew1[(l.tp.int+1):l.fp.int,], file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}
write.table(paste("\\end{tabular}\\end{table}"),file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

write.table(paste("\\begin{table}[!h]\\caption{FDR}\\label{}\\centering\\begin{tabular}{|c|c|c|}\\hline"), file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(c(paste('FDR  & additive FDR & interaction FDR\\ \\hline')), file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
FDR.vec=c(FDR,FDR.a,FDR.i)
FDR.ttab=tex.tab(FDR.vec)
write.table(FDR.ttab, file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"),file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))


sorted.add=NULL
if(!is.null(nrow(tabnew_add)))
sorted.add=tabnew_add[order(as.numeric(tabnew_add[,2]),tabnew_add[,1],decreasing=TRUE),]

sorted.int=tabnew[order(as.numeric(tabnew[,2]),tabnew[,1],decreasing=TRUE),]
tabnew1=NULL
if(is.matrix(sorted.add)){
if(nrow(sorted.add)>0)
tabnew1=as.matrix(apply(sorted.add,1,tex.tab))
}
if(!is.matrix(sorted.add)){
if(length(sorted.add)>0)
tabnew1=tex.tab(sorted.add)
}
write.table(paste("\\begin{tabular}{|c|c|c|c|}\\hline"), file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(c(paste('Expression & No.  & Frequency & Posterior probability\\ \\hline')), file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tabnew1, file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}"),file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

tabnew1=NULL
if(is.matrix(sorted.int)){
if(nrow(sorted.int)>0)
tabnew1=as.matrix(apply(sorted.int,1,tex.tab))
}
if(!is.matrix(sorted.int)){
if(length(sorted.int)>0)
tabnew1=tex.tab(sorted.int)
}
write.table(paste("\\begin{tabular}{|c|c|c|c|}\\hline"), file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(c(paste('Expression & No.  & Frequency & Posterior probability\\ \\hline')), file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tabnew1, file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}"),file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))


### power for true expressions
IA=NULL
tints=true.int1
if(length(IA)>0){tints=true.int1[-IA,]}else{tints=true.int1}


tints
t.ints=apply(tints,1,intpaste)
tp.ints=NULL
if(!is.null(tab.int.p)) tp.ints=tab.int.p[,1]
II=rep(0,length(t.ints))
for(d in 1:length(t.ints))
II[d]=ifelse(!is.na(match(t.ints[d],tp.ints)),match(t.ints[d],tp.ints),0)
nzi=NULL
nzi=which(II>0)
 pow.II=rep(0,nrow(tints))
 pow.II[nzi]=as.numeric(tab.int.p[II[nzi],2])/n.sym
P.tp.i=NULL 
P.tp.i.tex=NULL
if(length(pow.II)>0){
if(length(which(II>0))>1){
P.tp.i=cbind(
 as.vector(t(apply(tints,1,tree.to.int_Logic,n.snp=n.snp))),pow.II)
P.tp.i.tex=apply(P.tp.i,1,tex.tab)
}else{P.tp.i=c(tree.to.int_Logic(tints,n.snp=n.snp),pow.II)
P.tp.i.tex=tex.tab(P.tp.i)}
}
 
write.table(paste("\\begin{table}[!h]\\caption{Power for ETA}\\label{}\\centering\\begin{tabular}{|c|c|}\\hline"), file = paste(OutFile,"_Power_ETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tex.tab(c("Expression","Power")),file = paste(OutFile,"_Power_ETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(P.tp.i.tex,file = paste(OutFile,"_Power_ETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"), file = paste(OutFile,"_Power_ETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))


# E.fp=sum(c(FP.add[,2],FP.int[,8]))/n.sym

FDR1=NULL
# FDR1=sum(FP.add[,2])/(sum(FP.add[,2])+sum(TP.add[,2]))

twowaysFP=NULL
where2way=NULL
where2way=apply(FP.int,1,is.2way)
twowaysFP=FP.int[which(where2way),]
twowaysFP
twowaysTP=NULL
int2=NULL
where2way=NULL
where2way=apply(TP.int,1,is.2way)
twowaysTP=TP.int[which(where2way),]
twowaysTP

FDR2=NULL
FDR2= sum(twowaysFP[,8])/(sum(twowaysFP[,8])+sum(twowaysTP[,8]))
FDR2


threewaysFP=NULL
where3way=NULL
where3way=apply(FP.int,1,is.3way)
threewaysFP=FP.int[which(where3way),]
threewaysFP
threewaysTP=NULL
int3=NULL
where3way=NULL
where3way=apply(TP.int,1,is.3way)
threewaysTP=TP.int[which(where3way),]
threewaysTP

FDR3=NULL
FDR3= sum(threewaysFP[,8])/(sum(threewaysFP[,8])+sum(threewaysTP[,8]))
FDR3


fourwaysFP=NULL
where4way=NULL
where4way=apply(FP.int,1,is.4way)
fourwaysFP=FP.int[which(where4way),]
fourwaysFP
fourwaysTP=NULL
int4=NULL
where4way=NULL
where4way=apply(TP.int,1,is.4way)
fourwaysTP=TP.int[which(where4way),]
fourwaysTP

FDR4=NULL
FDR4= sum(fourwaysFP[,8])/(sum(fourwaysFP[,8])+sum(fourwaysTP[,8]))
FDR4


adds=NULL

addsFP=NULL
where1way=NULL
where1way=apply(FP.int,1,is.1way)
addsFP=FP.int[which(where1way),]
addsFP
addsTP=NULL
int1=NULL
where1way=NULL
where1way=apply(TP.int,1,is.1way)
addsTP=TP.int[which(where1way),]
addsTP
FDR1=sum(addsFP[,8])/(sum(addsFP[,8])+sum(addsTP[,8]))
FDR1
E.fp=sum(FP.int[,8])/n.sym
Efp1=sum(addsFP[,8])/n.sym
Efp1
Efp2=sum(twowaysFP[,8])/n.sym
Efp2
Efp3=sum(threewaysFP[,8])/n.sym
Efp3
Efp4=sum(fourwaysFP[,8])/n.sym
Efp4
EFp=c(Efp1,Efp2,Efp3,Efp4,E.fp)
efp=cbind(paste(c("$E(FP)_1$","$E(FP)_2$","$E(FP)_3$","$E(FP)_4$","$E(FP)$ Total"),sep=""),EFp)
efp=apply(efp,1,tex.tab)

FDr=cbind(paste(c("FDR1","FDR2","FDR3","FDR4","FDR Total"),sep=""),c(FDR1,FDR2,FDR3,FDR4,FDR))
FDr=apply(FDr,1,tex.tab)

write.table(paste("\\begin{table}[!h]\\caption{False discovery rate :ETA }\\label{}\\centering\\begin{tabular}{|c|c|}\\hline"), file = paste(OutFile,"_Power_ETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(FDr,file = paste(OutFile,"_Power_ETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"), file = paste(OutFile,"_Power_ETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

write.table(paste("\\begin{table}[!h]\\caption{Expected number of false positives :ETA }\\label{}\\centering\\begin{tabular}{|c|c|}\\hline"), file = paste(OutFile,"_Power_ETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(efp,file = paste(OutFile,"_Power_ETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"), file = paste(OutFile,"_Power_ETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))


#FDR.POWER.ETA(true.int,tab.add.p,tab.int.p,max.thres=0.5,OutFile=OutFile)

}
}



summarize.CETA=function(OutFile,A,p.param,n.sym,true.int,object,n.snp,n.obs,n.itr,boxplots.all=FALSE)
## function to summarize results of simulations for CETA approach
## calculates FDR and power, makes plots if applicable
## files containing tables in tex-syntax are created
{
fs=file.info(paste(OutFile,".txt",sep=""))$size
if(fs>0){		
		M1.1<-NULL
		M1.1<-read.table(file=paste(OutFile,".txt",sep=""), header = FALSE, sep = "\t", quote = "\"'", dec = ".",as.is = TRUE,na.strings = "NA", colClasses = NA, nrows = -1,
				   skip = 0, check.names = TRUE, strip.white = FALSE, blank.lines.skip = TRUE,
				   comment.char = "#",
				   allowEscapes = FALSE, flush = FALSE,
				   stringsAsFactors = FALSE,
				   fileEncoding = "", encoding = "unknown")
add.eff<-NULL
int2<-NULL
int3<-NULL
str2 <- gsub(' {2,}',' ',M1.1[,1])
int.level=unlist(lapply(strsplit(str2,' '),length))
add.eff<-M1.1[which(int.level==1),]
int=M1.1[which(int.level>1),]

OutFile.a=paste(OutFile,"_add")
M1.2<-NULL
		M1.2<-read.table(file=paste(OutFile.a,".txt",sep="")#,"_add.txt",sep="")
		, header = FALSE, sep = "\t", quote = "\"'", dec = ".",as.is = TRUE,na.strings = "NA", colClasses = NA, nrows = -1,
				   skip = 0, check.names = TRUE, strip.white = FALSE, blank.lines.skip = TRUE,
				   comment.char = "#",
				   allowEscapes = FALSE, flush = FALSE,
				   stringsAsFactors = FALSE,
				   fileEncoding = "", encoding = "unknown")

add.eff=M1.2
res.post.add=NULL
if(nrow(add.eff)>0){
Lv=count.rows(add.eff[,1])
L.p=Lv[,2]
res.post.add=vector("list",length(L.p))
for(i in 1:length(L.p)){
p=NULL
if(length(which(Lv[,3]==i))>0) p=add.eff[L.p[which(Lv[,3]==i)][[1]],2]
res.post.add[[i]]=list(v=as.numeric(Lv[i,3]),n=Lv[i,1],p=p)
}
}


  write.table(int, file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  Alli=read.table(file = "M.int.txt")
   fn <- "M.int.txt"
 if (file.exists(fn)) file.remove(fn)

Alli1=t(apply(Alli[,1:7],1,tree.merge))
Alli1.vec=apply(Alli1[,1:7],1,intpaste)
int.1=cbind(Alli1.vec,int[,2:5])
int=int.1

res.post.int=NULL
Lvi=NULL
if(nrow(int)>0){
Lvi=count.rows(Alli1.vec)
res.post.int=vector("list",nrow(Lvi))
res.post.int=apply(Lvi,1,do.list.cumTP,int=int)
}


	

namesL1<-NULL
ind1<-2*c(1:n.snp)-1
ind2<-2*c(1:n.snp)
i=c(1:n.snp)
namesL1[ind1]<-paste("D[",i,"]",sep="")
namesL1[ind2]<-paste("R[",i,"]",sep="")
################ BOXPLOTS
mint.v<-NULL
Ll=NULL
if(!is.null(true.int)){
write.table(true.int, file= "trint.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
mint.v<-read.table(file = "trint.txt")
fn <- "trint.txt"
if (file.exists(fn)) file.remove(fn)
va=t(apply(mint.v,1,variables.from.model))
vva=sort(as.vector(unlist(va)))
TNa<-ifelse(vva<0,paste(namesL1[abs(vva)],"^C",sep=""),paste(namesL1[as.numeric(abs(vva))],sep=""))
ex.nam<-TNa

if(!is.null(res.post.add)){
Ll=res.post.add[vva]
if(boxplots.all){
for(ss in 1:length(vva)){
jpeg(file=file.path(OutFile,"_boxplot_add",ss,".jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
allps=Ll[[ss]]$p
if(length(Ll[[ss]]$p)<n.sym) allps=c(Ll[[ss]]$p,rep(0,(n.sym-length(Ll[[ss]]$p))))
boxplot(allps,
        col = "yellow",
        main = paste("Boxplot of posterior probability for",ex.nam[ss]),
        xlab = " ",
        ylab = "Posterior probability",
        ylim = c(0, 2), yaxs = "i")
		axis(1, c(1), paste(ex.nam[ss],": Number of appearances :",Ll[[ss]]$n,", Number of post.pr >=0.5 :",length(which(allps>=0.5)),sep=""), cex.lab=1.2, cex.axis=0.4, cex.main=0.5, cex.sub=0.5, col.axis = "black",las=1)
dev.off()
}
}
}
}
####### boxplots-altogether

ps=NULL
res=NULL
for(ss in 1:length(vva)){
ps=Ll[[ss]]$p
if(length(Ll[[ss]]$p)<n.sym) ps=c(Ll[[ss]]$p,rep(0,(n.sym-length(Ll[[ss]]$p))))
res=cbind(res,ps)
}
colnames(res)=NULL
jpeg(file=file.path(OutFile,"_boxplot_add_all.jpeg",fsep=""),height=12,width=24,unit="cm",res=300)
boxplot(res[,1:ncol(res)], 
        names=paste(ex.nam[1:ncol(res)],sep=""),
        col = "yellow",
        main = paste("Boxplots of posterior probabilities"),
        xlab = " ",
        ylab = "Posterior probability",
        ylim = c(0, 1.2), yaxs = "i",las=2)
dev.off()
Ll=NULL
true.int1=NULL

###############################################################

  
# ################################
if(!is.null(true.int)){
 write.table(true.int, file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  true.int1=read.table(file = "M.int.txt")
   fn <- "M.int.txt"
 if (file.exists(fn)) file.remove(fn)
 }
  Alli=NULL
  ind.tr.i=NULL
   if(!is.null(res.post.int)){
   write.table(Lvi[,3], file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  Alli=read.table(file = "M.int.txt")
   fn <- "M.int.txt"
 if (file.exists(fn)) file.remove(fn)
  }
  ind.tr.i=NULL
Pom=res.post.int
  
  if(!is.null(true.int1)){
  if(length(unique(unlist(lapply(apply(true.int1,1,variables.from.model),length))))==1){
  v=t(apply(true.int1,1,variables.from.model))
  if(nrow(v)>1){
  NO.Int=NULL
  NO.Int=nrow(true.int1)
  ind.tr.i=vector("list",nrow(true.int1))
  for(ij in 1:nrow(true.int1)){
  #ind.tr.i[[ij]]=NULL
  sumi=which(rowSums(Alli)==sum(true.int1[ij,]))
  m1=Alli[which(rowSums(Alli)==sum(true.int1[ij,])),]
  kt=NULL
  for(rr in 1:nrow(m1)) kt=c(kt,(sum(v[ij,]%in%m1[rr,])==length(v[ij,])))
  if(length(which(kt))>0) ind.tr.i[[ij]]=sumi[which(kt)]
  }}
  }else{v=apply(true.int1,1,variables.from.model)# a list
		NO.Int=NULL
		NO.Int=length(which(unlist(lapply(v,length))>1))
		ind.tr.i=vector("list",NO.Int)
		for(ij in 1:NO.Int){
		  ind.tr.i[ij]=0
		  sumi=which(as.vector(rowSums(Alli)==sum(true.int1[ij,])))
		  m1=Alli[which(rowSums(Alli)==sum(true.int1[ij,])),]
		  kt=NULL
		  for(rr in 1:nrow(m1)) kt=c(kt,(sum(v[[ij]]%in%m1[rr,])==length(v[[ij]])))
		  if(length(which(kt))>0) ind.tr.i[[ij]]=sumi[which(kt)]
		  }
	  }
  }
  
  
  
  
  
  
  #res.post.int[[sumi[which(kt)]]]$p# w tym wierszu 
  m1.int.v<-NULL
  write.table(true.int, file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  m1.int.v<-read.table(file = "M.int.txt")
   fn <- "M.int.txt"
 if (file.exists(fn)) file.remove(fn)
  ex.nam.int<-apply(m1.int.v,1,tree.to.int_Logic,n.snp=n.snp)
if(!is.null(ind.tr.i)){ 
for(si in 1:length(ind.tr.i)){ 
if(!is.null(ind.tr.i[si][[1]])){
if(ind.tr.i[si][[1]]>0){
jpeg(file=file.path(OutFile,"_boxplot_int",si,".jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
boxplot(Pom[[ind.tr.i[si][[1]]]]$p2,Pom[[ind.tr.i[si][[1]]]]$p3,Pom[[ind.tr.i[si][[1]]]]$p4,
		col = "yellow",
		xlab = paste("Expression ",ex.nam.int[si],#"(",Pom[[ind.tr.i[si]]]$n,")",
		sep=""),
        ylab = "Posterior probability",
        ylim = c(0, 1.5), yaxs = "i")
		axis(1, c(1:3), paste(c("Two-Leaf","Three-Leaf","Four-Leaf"),sep=""), cex.lab=1.2, cex.axis=0.4, cex.main=1.5, cex.sub=1.5, col.axis = "black",las=2)
dev.off()
}
}
}
}


i=NULL
s=NULL
ex.nam=NULL
ex.nam.int=NULL
m1.int.v=NULL

# ############################################################
# ############################################################

################## BOXPLOT OF TREE SIZES
i=NULL
s=NULL
ex.nam=NULL
ex.nam.int=NULL
m1.int.v=NULL
##################
tab.add<-NULL
if(length(res.post.add)>0){
tab.add<-matrix(nrow=length(res.post.add),ncol=3)
for(i in 1:length(res.post.add)){ 
tab.add[i,]<-c(res.post.add[[i]]$v,res.post.add[[i]]$n,sum(as.numeric(res.post.add[[i]]$p)))
}
}

tab.int<-NULL
if(!is.null(res.post.int)){
tab.int<-matrix(nrow=length(res.post.int),ncol=6)
for(i in 1:length(res.post.int)){
NN=res.post.int[[i]]$n$counts
tab.int[i,]<-c(res.post.int[[i]]$v$x,NN,sum(as.numeric(res.post.int[[i]]$p2)),sum(as.numeric(res.post.int[[i]]$p3)),sum(as.numeric(res.post.int[[i]]$p4)),sum(as.numeric(res.post.int[[i]]$p5)))
}
}

 tab.add[,3]<-as.numeric(tab.add[,3])/n.sym
 tab.int[,3:6]<-as.numeric(tab.int[,3:6])/n.sym 
 if(!is.null(tab.add)) 
 if(length(which(as.numeric(tab.add[,3])>=0))>0)
	{tab.add=tab.add[which(as.numeric(tab.add[,3])>=0),]
	}else{tab.add=NULL}
 if(!is.null(tab.int)){
	sum.post=rowSums(cbind(as.numeric(tab.int[,3]),as.numeric(tab.int[,4]),as.numeric(tab.int[,5])))
	pos.i=which(sum.post>0)
 if(length(pos.i)>0){
 tab.int=tab.int[pos.i,]}else{tab.int=NULL}
 }
 ind=NULL
 No.a=NULL
for(jj in 1:length(res.post.add)){ 
  if(length(which(res.post.add[[jj]]$p>=0.5))>0){
  ind=c(ind,jj)
  No.a=c(No.a,length(which(res.post.add[[jj]]$p>=0.5)))
  }
 }
 tab.add.p=cbind(tab.add[ind,1],No.a,tab.add[ind,3])
 N.A=NULL
 N.A=sum(No.a)
 ind2=NULL
 ind3=NULL
 ind4=NULL
 ind5=NULL
 No.i2=NULL
 No.i3=NULL
 No.i4=NULL
 No.i5=NULL
 tab.int.p2=NULL
 tab.int.p3=NULL
 tab.int.p4=NULL
 tab.int.p=NULL
 N.2=0
N.3=0
N.4=0
N.5=0
if(!is.null(res.post.int)){ 
for(jj in 1:length(res.post.int)){ 
  if(length(which(res.post.int[[jj]]$p2>=0.5))>0){
  ind2=c(ind2,jj)
  }
  No.i2=c(No.i2,length(which(res.post.int[[jj]]$p2>=0.5)))
  if(length(which(res.post.int[[jj]]$p3>=0.5))>0){
  ind3=c(ind3,jj)
  }
  No.i3=c(No.i3,length(which(res.post.int[[jj]]$p3>=0.5)))
  if(length(which(res.post.int[[jj]]$p4>=0.5))>0){
  ind4=c(ind4,jj)
  }
  No.i4=c(No.i4,length(which(res.post.int[[jj]]$p4>=0.5)))
  if(length(which(res.post.int[[jj]]$p5>=0.5))>0){
  ind5=c(ind5,jj)
  }
  No.i5=c(No.i5,length(which(res.post.int[[jj]]$p5>=0.5)))
 }
 if(!is.null(ind2)) tab.int.p2=cbind(tab.int[ind2,1],No.i2[which(No.i2>0)],tab.int[ind2,3])
 if(!is.null(ind3)) tab.int.p3=cbind(tab.int[ind3,1],No.i3[which(No.i3>0)],tab.int[ind3,4])
 if(!is.null(ind4)) tab.int.p4=cbind(tab.int[ind4,1],No.i4[which(No.i4>0)],tab.int[ind4,5])
allind=unique(c(ind2,ind3,ind4))
  No.i=No.i2+No.i3+No.i4
 sumn=No.i[which(No.i>0)]
 if(!is.null(allind)) allind=sort(allind)
 sump=as.numeric(tab.int[allind,3])+as.numeric(tab.int[allind,4])+as.numeric(tab.int[allind,5])
 tab.int.p=cbind(tab.int[allind,1],sumn,tab.int[allind,3],tab.int[allind,4],tab.int[allind,5],tab.int[allind,6])
 }
 tab.int.p=tab.int.p[ ,c(1,2,6)]
################################
write.table(tab.int[allind,1], file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  tabint1<-read.table(file = "M.int.txt")
   fn <- "M.int.txt"
 if (file.exists(fn)) file.remove(fn)
tabinT1=t(apply(tabint1,1,tree.merge))
CRt1=count.rows(tabinT1)
whichrows=CRt1[2]
nwr=nrow(whichrows)
cmdp2<-paste("ptp2 <- sum(as.numeric(tab.int[allind,3])[whichrows[[1]]$'",1:nwr,"'])",sep="")
cmdp3<-paste("ptp3 <- sum(as.numeric(tab.int[allind,4])[whichrows[[1]]$'",1:nwr,"'])",sep="")
cmdp4<-paste("ptp4 <- sum(as.numeric(tab.int[allind,5])[whichrows[[1]]$'",1:nwr,"'])",sep="")
cmdNo<-paste("SumNo <- sum(sumn[whichrows[[1]]$'",1:nwr,"'])",sep="")
ptp2=NULL
ptp3=NULL
ptp4=NULL
ptp2=unlist(lapply(c(1:nwr),function(i) eval(parse(text=cmdp2[i]))))
ptp3=unlist(lapply(c(1:nwr),function(i) eval(parse(text=cmdp3[i]))))
ptp4=unlist(lapply(c(1:nwr),function(i) eval(parse(text=cmdp4[i]))))
ptp5=ptp2+ptp3+ptp4
SumNo=unlist(lapply(c(1:nwr),function(i) eval(parse(text=cmdNo[i]))))
TABp=NULL
TABp=cbind(as.matrix(CRt1[,3:9]),SumNo,ptp2,ptp3,ptp4,ptp5)


tab.int.p=cbind(apply(TABp[,1:7],1,intpaste),SumNo,ptp5)

######### FALSE POSITIVES AND TRUE POSITIVES
res.tpfp<-tp.fp.CETA(true.int,tab.add.p,tab.int.p,thres=0.85)
TP.add=res.tpfp$TP.add
TP.int=res.tpfp$TP.int
FP.add=res.tpfp$FP.add
FP.int=res.tpfp$FP.int
######### FALSE DISCOVERY RATE
FDR=res.tpfp$FDR
FDR.a=res.tpfp$FDR.a
FDR.i=res.tpfp$FDR.i 
###############################################
#######		plots and tables with summaries
###############################################
namesL<-NULL
ind1<-2*c(1:n.snp)-1
ind2<-2*c(1:n.snp)
i=c(1:n.snp)
namesL[ind1]<-paste("D_{",i,"}",sep="")
namesL[ind2]<-paste("R_{",i,"}",sep="")
#####################################################
######## table with settings to the file 
#####################################################
set.tab=NULL
set.tab<-rbind(c(paste('True model:'),paste("$",object$snpmodel,"$")),
c(paste("No. of snps"), n.snp),c(paste("No. of MC iterations"),n.itr),c(paste("'a' parameter"),A),c(paste("'p' parameter"), p.param),c(paste("No of simulations"),n.sym),c(paste("No of individuals"),n.obs),c(paste("Max no. of leaves per tree"), (n.lvs+1)/2),c(paste("Max no. of trees"), n.trs))
newset.tab=as.matrix(apply(set.tab,1,tex.tab))
write.table(paste("\\begin{table}[!h] \\caption{Settings }\\label{} \\centering \\begin{tabular}{p{4cm}|p{8cm}}\\hline"), file = paste(OutFile,"_tab.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(newset.tab, file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"),file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

write.table(paste("\\begin{table}[!h]\\caption{Selected expressions}\\label{}\\centering\\begin{tabular}{|c|c|c|c|}\\hline"), file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
tp.add=TP.add[,1]
fp.add=FP.add[,1]
library(plotrix)
Intnames<- c(tp.add,fp.add)
PostP<-NULL
if(!is.null(tab.add.p)){
if(is.matrix(tab.add.p)) PostP<-c(as.numeric(tab.add.p[,3]))
if(!is.matrix(tab.add.p)) PostP<-c(as.numeric(tab.add.p[3]))
}
T1<- NULL
T1<-c(tp.add,fp.add)
T1<-as.numeric(T1)
TNa<-ifelse(T1<0,paste(namesL[abs(T1)],"^C",sep=""),paste(namesL[as.numeric(abs(T1))],sep=""))
ex.nam<-TNa
tabnew_add=NULL
#######################################################
write.table(c(paste('Expression & No.  & Frequency & Posterior probability\\ \\hline')), file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
if(length(ex.nam)>0){
 if(is.matrix(tab.add.p)) tabnew_add=cbind(paste("$",ex.nam,"$"),tab.add.p[,2],round(as.numeric(tab.add.p[,2])/n.sym,digits=5),tab.add.p[,3])
 if(!is.matrix(tab.add.p)) tabnew_add=cbind(paste("$",ex.nam,"$"),tab.add.p[2],round(as.numeric(tab.add.p[2])/n.sym,digits=5),tab.add.p[3])
colnames(tabnew_add)=NULL
tabnew_add1=NULL
if(is.matrix(tabnew_add)){
if(nrow(tabnew_add)>0)
tabnew_add1=as.matrix(apply(tabnew_add,1,tex.tab))
}
if(!is.matrix(tabnew_add)){
if(length(tabnew_add)>0)
tabnew_add1=tex.tab(tabnew_add)
}
write.table(tabnew_add1, file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"),file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}
############# tex tables with fp and tp for add effects
l.tp.add=0
l.fp.add=0
if(!is.null(TP.add)) l.tp.add=nrow(TP.add)
if(!is.null(FP.add)) l.fp.add=nrow(FP.add)
write.table(paste("\\begin{table}[!h]\\caption{True positives}\\label{}\\centering\\begin{tabular}{|c|c|c|c|c|c|c|}\\hline"), file = paste(OutFile,"_tp_tab.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(c(paste("Expression & No.  & Frequency & Two-way & Three-way & Four-way & Total\\ \\hline")), file = paste(OutFile,"_tp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
if(l.tp.add>0)
write.table(tabnew_add1[1:l.tp.add,], file = paste(OutFile,"_tp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

write.table(paste("\\begin{table}[!h]\\caption{False positives}\\label{}\\centering\\begin{tabular}{|c|c|c|c|}\\hline"), file = paste(OutFile,"_fp_tab.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(c(paste('Expression & No.  & Frequency & Posterior probability\\ \\hline')), file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
if(l.tp.add==0){
if(l.fp.add>0)
write.table(tabnew_add1[1:l.fp.add,], file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))}else{
if(l.fp.add>0)
write.table(tabnew_add1[(l.tp.add+1):l.fp.add,], file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}
####################################################
m.tp.int=NULL
m.fp.int=NULL
if(!is.null(TP.int)) if(nrow(TP.int)>0) m.tp.int=TP.int[,1:7]
if(!is.null(FP.int)) if(nrow(FP.int)>0) m.fp.int=FP.int[,1:7]
# m.int<-rbind(m.tp.int,m.fp.int)
m.int<-NULL
if(!is.null(m.tp.int)&&!is.null(m.fp.int)){m.int=rbind(as.matrix(m.tp.int),as.matrix(m.fp.int))}else{m.int=mat.bind(m.tp.int,m.fp.int)}

if(!is.null(m.int)){ 
   PP<-NULL
    # colnames(FP.int)=NULL
	# colnames(TP.int)=NULL
	# tab.int.p1=rbind(as.matrix(FP.int),as.matrix(TP.int))
	# tab.int.p2=cbind(apply(tab.int.p1[,1:7],1,intpaste),tab.int.p1[,8:9])
	# tab.int.p<-tab.int.p2
  if(!is.null(nrow(tab.int.p))){ 
  PP<-as.double(as.numeric(tab.int.p[,3]))
  PostP<-c(PostP,PP)
  m1.int.v<-NULL
  write.table(tab.int.p[,1], file = "al.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  }
  if(is.null(nrow(tab.int.p))){
  PP<-as.double(as.numeric(tab.int.p[3]))
  PostP<-c(PostP,PP)
  m1.int.v<-NULL
  write.table(tab.int.p[1], file = "al.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  }
  m1.int.v<-read.table(file = "al.int.txt")
  fn <- "al.int.txt"
  if (file.exists(fn)) file.remove(fn)
  ex.nam<-c(ex.nam,apply(m1.int.v,1,tree.to.int_Logic,n.snp=n.snp))
}

b.s<-NULL
b.s<-cbind(ex.nam,PostP)
if(nrow(b.s)>0){
jpeg(file=file.path(OutFile,"_Posterior_prob.jpeg",fsep=""),height=12,width=12,unit="cm",res=300)
xb<-c(1:nrow(b.s))
plot(x=xb,y=PostP,lty=1,lwd=2,type="h",axes=FALSE,xlab=NA,ylab="Posterior probability",
main=paste("Posterior probability for logic expressions.\nTrue model:",object$snpmodel),
ylim=c(0.0,1.5),cex.lab=1.0, cex.axis=0.5, cex.main=0.4, cex.sub=0.5,col="blue")
    axis(1, xb, ex.nam, cex.lab=1.2, cex.axis=0.4, cex.main=0.5, cex.sub=0.5, col.axis = "black",las=2)
	  axis(2, at=seq(0,1,0.1), cex.lab=1.2, cex.axis=1.0, cex.main=0.5, cex.sub=0.5,col.axis = "black")
    abline(h=0.5,lty=1,col="red")
  legend(x=mean(xb),y=1.5,c(paste("No. of snps=", n.snp),paste("No. of MC iterations=",n.itr),paste("'a' parameter=",A),paste("'p' parameter=", p.param),paste("No of simulations=",n.sym),paste("No of individuals=",n.obs),paste("Max no. of leaves per tree=", (n.lvs+1)/2),paste("Max no. of trees=", n.trs)),pch=c(1,1,1,1,1),col=c("black","black","black","black","black"),inset = .01,cex=0.4,bty="n", title.adj=0.15,title="Simulation settings")
dev.off()
}
#########################################
Intnames<- NULL
PostP<-NULL
ex.nam<-NULL
b.s<-NULL

m.int<-NULL
if(!is.null(m.tp.int)&&!is.null(m.fp.int)){m.int=rbind(as.matrix(m.tp.int),as.matrix(m.fp.int))}else{m.int=mat.bind(m.tp.int,m.fp.int)}

#########################################
tabnew=NULL
tabnew1=NULL
if(!is.null(TP.int)) colnames(TP.int)=c(1:ncol(TP.int))
if(!is.null(FP.int)) colnames(FP.int)=c(1:ncol(FP.int))
whichii=NULL
for(ww in 1:nrow(TP.int)){
whichii=c(whichii,which(as.numeric(tab.int.p[,3])==TP.int[ww,9]))
}
whichii

tab.int.P=NULL 
tab.int.P<-mat.bind(TP.int,FP.int)
tab.int.pp=NULL
if(!is.null(tab.int.P))
if(nrow(tab.int.P)>0){
tab.int.P=as.matrix(tab.int.P)
tab.int.pp=t(matrix(as.numeric(t(tab.int.P[,1:7])),nrow=7))
}
if(!is.null(tab.int.pp)){
if(is.matrix(tab.int.p)){
if(nrow(tab.int.p)>0){
ex.nam<-apply(tab.int.pp,1,tree.to.int_Logic,n.snp=n.snp)
tab.int.new=rbind(cbind(TABp[whichii,1:11],as.numeric(tab.int.p[whichii,3])),TABp[-whichii,])
tabnew=cbind(paste("$",ex.nam,"$"),tab.int.new[,8],round(as.numeric(tab.int.new[,8])/n.sym,digits=5),tab.int.new[,9],tab.int.new[,10],tab.int.new[,11],tab.int.new[,12])
colnames(tabnew)=NULL
tabnew1=NULL
if(is.matrix(tabnew)){
if(nrow(tabnew)>0)
tabnew1=as.matrix(apply(tabnew,1,tex.tab))
}
if(!is.matrix(tabnew)){
if(length(tabnew)>0)
tabnew1=tex.tab(tabnew)
}
write.table(paste("\\begin{table}[!h]\\caption{Selected expressions}\\label{}\\centering\\begin{tabular}{|c|c|c|c|c|c|c|}\\hline"), file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(c(paste('Expression & No.  & Frequency & Two-way & Three-way & Four-way & Total\\ \\hline')), file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tabnew1, file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}
}

if(!is.matrix(tab.int.p)){
if(length(tab.int.p)>0){
write.table(tab.int.p[1], file = "al.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  m1.int.v<-read.table(file = "al.int.txt")
  fn <- "al.int.txt"
  if (file.exists(fn)) file.remove(fn)
  ex.nam<-apply(m1.int.v,1,tree.to.int_Logic,n.snp=n.snp)
tabnew=cbind(paste("$",ex.nam,"$"),tab.int.p[2],round(as.numeric(tab.int.p[2])/n.sym,digits=5),tab.int.p[3])
colnames(tabnew)=NULL
tabnew1=NULL
if(is.matrix(tabnew)){
if(nrow(tabnew)>0)
tabnew1=as.matrix(apply(tabnew,1,tex.tab))
}
if(!is.matrix(tabnew)){
if(length(tabnew)>0)
tabnew1=tex.tab(tabnew)
}

write.table(tabnew1, file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}
}
}
write.table(paste("\\end{tabular}\\end{table}"),file = paste(OutFile,"_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
############# tex tables with fp and tp for int effects
l.tp.int=0
l.fp.int=0
if(!is.null(TP.int)) l.tp.int=nrow(TP.int)
if(!is.null(FP.int)) l.fp.int=nrow(FP.int)
if(l.tp.int>0)
write.table(tabnew1[1:l.tp.int,], file = paste(OutFile,"_tp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"),file = paste(OutFile,"_tp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

if(l.tp.int==0){
if(l.fp.int>0)
write.table(tabnew1[1:l.fp.int,], file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))}else{
if(l.fp.int>0)
write.table(tabnew1[(l.tp.int+1):l.fp.int,], file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}
write.table(paste("\\end{tabular}\\end{table}"),file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

write.table(paste("\\begin{table}[!h]\\caption{FDR}\\label{}\\centering\\begin{tabular}{|c|c|c|}\\hline"), file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(c(paste('FDR  & additive FDR & interaction FDR\\ \\hline')), file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
FDR.vec=c(FDR,FDR.a,FDR.i)
FDR.ttab=tex.tab(FDR.vec)
write.table(FDR.ttab, file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"),file = paste(OutFile,"_fp_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))



sorted.add=tabnew_add[order(as.numeric(tabnew_add[,2]),tabnew_add[,1],decreasing=TRUE),]

sorted.int=tabnew[order(as.numeric(tabnew[,2]),tabnew[,1],decreasing=TRUE),]
tabnew1=NULL
if(is.matrix(sorted.add)){
if(nrow(sorted.add)>0)
tabnew1=as.matrix(apply(sorted.add,1,tex.tab))
}
if(!is.matrix(sorted.add)){
if(length(sorted.add)>0)
tabnew1=tex.tab(sorted.add)
}
write.table(paste("\\begin{tabular}{|c|c|c|c|}\\hline"), file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(c(paste('Expression & No.  & Frequency & Posterior probability\\ \\hline')), file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tabnew1, file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}"),file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

tabnew1=NULL
if(is.matrix(sorted.int)){
if(nrow(sorted.int)>0)
tabnew1=as.matrix(apply(sorted.int,1,tex.tab))
}
if(!is.matrix(sorted.int)){
if(length(sorted.int)>0)
tabnew1=tex.tab(sorted.int)
}
write.table(paste("\\begin{table}[!h]\\caption{sorted results}\\label{}\\centering\\begin{tabular}{|c|c|c|c|}\\hline"), file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(c(paste('Expression & No.  & Frequency & Posterior probability\\ \\hline')), file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tabnew1, file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"),file = paste(OutFile,"_sorted_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))



### power for CETA for true expressions

## adds
adds=true.int1[which(apply(true.int1,1,is.1way)),1]


if(length(adds)>0){
i.pow=NULL
pow.tp.a=NULL
for(s in 1:length(adds))
{
i.pow=c(i.pow,which(tab.add.p[,1]==adds[s]))
if(length(i.pow)>0) pow.tp.a=c(pow.tp.a,tab.add.p[i.pow[s],2]/n.sym)
}

P.tp.a=NULL
P.tp.a.tex=NULL
if(!is.null(pow.tp.a)){
P.tp.a=cbind(namesL[adds],pow.tp.a)
colnames(P.tp.a)=c("True expression","Power")
P.tp.a.tex=apply(P.tp.a,1,tex.tab)


write.table(paste("\\begin{table}[!h]\\caption{Power for CETA}\\label{}\\centering\\begin{tabular}{|c|c|}\\hline"), file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tex.tab(c("Expression","Power")),file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(P.tp.a.tex,file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"), file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
}
}
# P.tp.a.tex
# P.tp.a

# length(which(Pom[[ind.tr.i[d][[1]]]]$p5>0.5))/nsym
### power for true expressions
IA=which(apply(true.int1,1,is.1way))
if(length(IA)>0){tints=true.int1[-IA,]}else{tints=true.int1}
tints
t.ints=apply(tints,1,intpaste)
tp.ints=tab.int.p[,1]#apply(TP.int[,1:7],1,intpaste)
# II=rep(0,length(t.ints))
 pow.II=NULL
for(d in 1:length(t.ints))
if(!is.null(ind.tr.i[d][[1]]))#||
 #if(ind.tr.i[d][[1]]>0)
 {pow.II=c(pow.II,length(which(Pom[[ind.tr.i[d][[1]]]]$p5>0.5))/n.sym)}else{pow.II=c(pow.II,0)}
P.tp.i=NULL 
P.tp.i.tex=NULL
if(length(pow.II)>0){
# if(length(which(II>0))>1){
P.tp.i=cbind(
 as.vector(t(apply(tints,1,tree.to.int_Logic,n.snp=n.snp))),pow.II)
P.tp.i.tex=apply(P.tp.i,1,tex.tab)
# }else{P.tp.i=c(tree.to.int_Logic(tints,n.snp=n.snp),pow.II) P.tp.i.tex=tex.tab(P.tp.i)}
}
 
 
write.table(paste("\\begin{table}[!h]\\caption{Power for CETA}\\label{}\\centering\\begin{tabular}{|c|c|}\\hline"), file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tex.tab(c("Expression","Power")),file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(P.tp.i.tex,file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"), file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))


 
 
write.table(paste("\\begin{table}[!h]\\caption{Power for CETA}\\label{}\\centering\\begin{tabular}{|c|c|}\\hline"), file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(tex.tab(c("Expression","Power")),file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(P.tp.i.tex,file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"), file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))



E.fp=sum(c(FP.add[,2],FP.int[,8]))/n.sym

FDR1=NULL
FDR1=sum(FP.add[,2])/(sum(FP.add[,2])+sum(TP.add[,2]))

twowaysFP=NULL
where2way=NULL
where2way=apply(FP.int,1,is.2way)
twowaysFP=FP.int[which(where2way),]
twowaysFP
twowaysTP=NULL
int2=NULL
where2way=NULL
where2way=apply(TP.int,1,is.2way)
twowaysTP=TP.int[which(where2way),]
twowaysTP
# if(!is.matrix(twowaysFP)) twowaysFP=t(as.matrix(twowaysFP))
# if(!is.matrix(twowaysTP)) twowaysTP=t(as.matrix(twowaysTP))


FDR2=NULL
FDR2= sum(twowaysFP[,8])/(sum(twowaysFP[,8])+sum(twowaysTP[,8]))
FDR2


threewaysFP=NULL
where3way=NULL
where3way=apply(FP.int,1,is.3way)
threewaysFP=FP.int[which(where3way),]
threewaysFP
threewaysTP=NULL
int3=NULL
where3way=NULL
where3way=apply(TP.int,1,is.3way)
threewaysTP=TP.int[which(where3way),]
threewaysTP
if(!is.matrix(threewaysFP)) threewaysFP=t(as.matrix(threewaysFP))
if(!is.matrix(threewaysTP)) threewaysTP=t(as.matrix(threewaysTP))

FDR3=NULL
FDR3= sum(threewaysFP[,8])/(sum(threewaysFP[,8])+sum(threewaysTP[,8]))
FDR3


fourwaysFP=NULL
where4way=NULL
where4way=apply(FP.int,1,is.4way)
fourwaysFP=FP.int[which(where4way),]
fourwaysFP
fourwaysTP=NULL
int4=NULL
where4way=NULL
where4way=apply(TP.int,1,is.4way)
fourwaysTP=TP.int[which(where4way),]
fourwaysTP
if(!is.matrix(fourwaysFP)) fourwaysFP=t(as.matrix(fourwaysFP))
if(!is.matrix(fourwaysTP)) fourwaysTP=t(as.matrix(fourwaysTP))



FDR4=NULL
FDR4= sum(fourwaysFP[,8])/(sum(fourwaysFP[,8])+sum(fourwaysTP[,8]))
FDR4


adds=NULL

addsFP=NULL
where1way=NULL
where1way=apply(FP.int,1,is.1way)
addsFP=FP.int[which(where1way),]
addsFP
addsTP=NULL
int1=NULL
where1way=NULL
where1way=apply(TP.int,1,is.1way)
addsTP=TP.int[which(where1way),]
addsTP
# FDR1=sum(addsFP[,8])/(sum(addsFP[,8])+sum(addsTP[,8]))
# FDR1
# E.fp=sum(FP.int[,8])/n.sym
Efp1=sum(FP.add[,2])/n.sym
Efp1
Efp2=sum(twowaysFP[,8])/n.sym
Efp2
Efp3=sum(threewaysFP[,8])/n.sym
Efp3
Efp4=sum(fourwaysFP[,8])/n.sym
Efp4

EFp=c(Efp1,Efp2,Efp3,Efp4,E.fp)
efp=cbind(paste(c("$E(FP)_1$","$E(FP)_2$","$E(FP)_3$","$E(FP)_4$","$E(FP)$ Total"),sep=""),EFp)
efp=apply(efp,1,tex.tab)

FDr=cbind(paste(c("FDR1","FDR2","FDR3","FDR4","FDR Total"),sep=""),c(FDR1,FDR2,FDR3,FDR4,FDR))
FDr=apply(FDr,1,tex.tab)

write.table(paste("\\begin{table}[!h]\\caption{False discovery rate :CETA }\\label{}\\centering\\begin{tabular}{|c|c|}\\hline"), file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(FDr,file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"), file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))

write.table(paste("\\begin{table}[!h]\\caption{Expected number of false positives :CETA }\\label{}\\centering\\begin{tabular}{|c|c|}\\hline"), file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(efp,file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
write.table(paste("\\end{tabular}\\end{table}"), file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))


# write.table(paste("E(FP)=",E.fp,sep=""), file = paste(OutFile,"_Power_CETA_tab.txt",sep=""), append = TRUE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
# # FDR.POWER.CETA(true.int,tab.add.p,tab.int.p,max.thres=0.5,OutFile=OutFile)



}
}

##########################################################################
do.list=function(VEC,int)
{
return(list(v=remove.levels(VEC[3]),n=VEC[1],p=as.numeric(int[VEC[2][[1]],2])))
}
##########################################################################
tvariables.from.model<-function(vec)
#function to return indices of variables in the model
{
vec0=vec[which(vec!=0)]
vecOp=vec0[which(vec0<1000)]
return(unique(vecOp))
}
##########################################################################
do.list.cumP=function(VEC,int)
{
return(list(v=remove.levels(VEC[3]),n=VEC[1],p2=as.numeric(int[as.numeric(VEC[2][[1]]),2]),p3=as.numeric(int[as.numeric(VEC[2][[1]]),3]),p4=as.numeric(int[as.numeric(VEC[2][[1]]),4])))
}
##########################################################################
compare.expr=function(expr1,exprList)
{
sum1=0
for(jj in 1:nrow(exprList))
if(sum(expr1%in%exprList[jj,])==7) sum1=sum1+1
if(sum1>0){return(TRUE)}else{return(FALSE)}
}
########################
##################################################################
make.tree=function(lop,lnull,n0v,lmv)
{
return(vec=c(rep(n0v[1],lop), n0v[2:(lmv+1)],rep(0,lnull)))
}	
##########################################################################
rem.null=function(vec)
## function to remove elements of a tree (vec) which are nor leaves nor operators
{
return(vec[which(vec!=0)])
}
##########################################################################
do.list.cumP=function(VEC,int)
{
return(list(v=remove.levels(VEC[3]),n=VEC[1],p2=as.numeric(int[as.numeric(VEC[2][[1]]),2]),p3=as.numeric(int[as.numeric(VEC[2][[1]]),3]),p4=as.numeric(int[as.numeric(VEC[2][[1]]),4])))
}
##########################################################################
do.list.cumTP=function(VEC,int)
{
return(list(v=remove.levels(VEC[3]),n=VEC[1],p2=as.numeric(int[as.numeric(VEC[2][[1]]),2]),p3=as.numeric(int[as.numeric(VEC[2][[1]]),3]),p4=as.numeric(int[as.numeric(VEC[2][[1]]),4]),p5=as.numeric(int[as.numeric(VEC[2][[1]]),5])))
}
##########################################################################
sum.pT.list=function(List1,f.i)
{
L2=NULL 
L2=lapply(seq_along(f.i),function(i)
        unlist(List1[[f.i[i]]]$p2))
p2=unlist(L2)
L3=NULL	
L3=lapply(seq_along(f.i),function(i)
        unlist(List1[[f.i[i]]]$p3))
p3=unlist(L3)	 
L4=NULL	
L4=lapply(seq_along(f.i),function(i)
        unlist(List1[[f.i[i]]]$p4))
p4=unlist(L4)	

L5=NULL	
L5=lapply(seq_along(f.i),function(i)
        unlist(List1[[f.i[i]]]$p5))
p5=unlist(L5)

N=NULL	
N=lapply(seq_along(f.i),function(i)
        unlist(List1[[f.i[i]]]$n$counts))
n=sum(unlist(N))	
 
return(list(n=n,p2=p2,p3=p3,p4=p4,p5=p5))
 }
##########################################################################
sum.p.list=function(List1,f.i)
{
L2=NULL 
L2=lapply(seq_along(f.i),function(i)
        unlist(List1[[f.i[i]]]$p2))
p2=unlist(L2)
L3=NULL	
L3=lapply(seq_along(f.i),function(i)
        unlist(List1[[f.i[i]]]$p3))
p3=unlist(L3)	 
L4=NULL	
L4=lapply(seq_along(f.i),function(i)
        unlist(List1[[f.i[i]]]$p4))
p4=unlist(L4)	

N=NULL	
N=lapply(seq_along(f.i),function(i)
        unlist(List1[[f.i[i]]]$n$counts))
n=sum(unlist(N))	
 
return(list(n=n,p2=p2,p3=p3,p4=p4))
 }
  
  
  
#  

Is.Op.vec=function(vec)
{
return(unlist(lapply(vec,is.operator)))
}



un.op.list=function(i,mat=PIi,ind.list=op.ind)
{
return(mat[i,ind.list[[i]]])
}




##########
data.sim.proc<-function(object,n.sym,n.snp,n.obs,n.lvs,n.trs,n.itr,a.start,a.end,lambda,p.param,n.tree.mod,tree.lvs,eff.sizes,slope,tree.ind,op,mypath=file.path("H:\\modMCLR\\package_update\\new_lib\\M1",fsep =""))
## whole simulation procedure for different settings of parameters
## for each  a summary files and plots are created
## makes summaries for MCLR and mMCLR for each 'a' at the same simulated data set
{
		#dir.create(mypath)##### creating a folder 
		nm<-n.snp
		X=snpmatrix(nsnp=n.snp,nrows=n.obs,p=p.param)
		MM=model_int(X=X,nlvs==n.lvs,n.tree.mod=n.tree.mod,tree.lvs=tree.lvs,eff.sizes=eff.sizes,slope=slope,tree.ind=tree.ind,op=op)
		true.int<-MM$m.int	
if(!is.null(true.int)){
 write.table(true.int, file = "M.int.txt", append = FALSE, quote = FALSE, sep = "\t ",eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"))
  true.int1=read.table(file = "M.int.txt")
  fn <- "M.int.txt"
  if (file.exists(fn)) file.remove(fn)
  }
true.int=apply(t(apply(as.matrix(true.int1),1,tree.merge)),1,intpaste)
		A<-a.start
		P<-p.param
while(A<=a.end){
		Out.file<-file.path(mypath,'//result_p=0',10*P,fsep="")
		sim.res=simulate(nsym=n.sym,nsnp=n.snp,nrows=n.obs,nitr=n.itr,ntrs=n.trs,nlvs=n.lvs,p=P,a=A,lambda=lambda,outFile=Out.file,n.tree.mod=n.tree.mod,tree.lvs=tree.lvs,eff.sizes=eff.sizes,slope=slope,tree.ind=tree.ind,op=op)	
A<-A+0.1
}
return(sim.res)
}
#  



				
mc.chain.summary<-function(filename,OutFile)
## function to summarize Markov-chains from a file "slogiclisting.tmp" being  a result of an original Kooperberg-Ruczinski 'logreg' function 
{
stats=read.table(filename)
model.weight=stats[,3]
stats=stats[,4:ncol(stats)]
stats.seq=data.frame(t(stats))
stats.seq=as.list(stats.seq)
ntrs=length(stats.seq$X1)/nlvs
T=lapply(stats.seq,vec.trees,nlvs=nlvs,ntrs=ntrs,indic=1)
nmlvs=apply(stats,1,vec.leaves)
nmtrs=apply(stats,1,vec.trees,nlvs=nlvs,ntrs=ntrs,indic=0)
nm=length(T)# number of models on the list 
#########################################
Iind=lapply(T,Int.tree.ind,ntrs=ntrs,nlvs=nlvs,lev=1)## a list with trees
PIi=do.call(rbind, T)#unlisted- as a matrix with trees in rows
var.list=apply(PIi,1,variables.from.model)
tree.sizes=unlist(lapply(var.list,length))
## histogram unweighted
hist(tree.sizes)
jpeg(file=file.path(OutFile,paste("Hist_mc-TreeSizes_k",ntrs,".jpeg",sep=""),fsep=""),height=12,width=12,unit="cm",res=300)
hist(tree.sizes,col = "lightblue", border = "pink",main=paste("Tree sizes for max ",ntrs,"  trees allowed",sep=""),xlab="Tree size")
# legend(x=max(tree.sizes)+1,y=max(count.rows(tree.sizes)[,1]),c(paste("Model: ",object$snpmodel,sep=""),paste("No. of snps=", n.snp),paste("No. of MC iterations=",n.itr),paste("'a' parameter=",A),paste("'p' parameter=", p.param),paste("No of simulations=",n.sym),paste("No of individuals=",n.obs),paste("Max no. of leaves per tree=", (n.lvs+1)/2)),pch=c(1,1,1,1,1),col=c("black","black","black","black","black"),inset = .01,cex=0.4,bty="n", title.adj=0.15,title="Simulation settings")

dev.off()
### weights (no of repetitions)
tree.weight=rep(model.weight,each=ntrs)
weighted.tree.sizes=unlist(mapply(rep,tree.sizes,each=tree.weight))
## histogram weighted
hist(weighted.tree.sizes)
jpeg(file=file.path(OutFile,paste("Hist_mc-TreeSizes_weighted_k",ntrs,".jpeg",sep=""),fsep=""),height=12,width=12,unit="cm",res=300)
hist(weighted.tree.sizes,col = "lightblue", border = "pink",main=paste("Tree sizes for max ",ntrs," trees allowed",sep=""),xlab="Weighted tree size")
dev.off()

}



