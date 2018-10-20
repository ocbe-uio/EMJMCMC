#load a required library
library("RCurl")
#read in the package most recent version
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")


#a function for postproceeding the features at the end of the run
simplifyposteriors=function(X,posteriors,th=0.0001,thf=0.2, resp)
{
  posteriors=posteriors[-which(posteriors[,2]<th),]
  rhash=hash()
  for(i in 1:length(posteriors[,1]))
  {
    expr=posteriors[i,1]
    res=model.matrix(data=X,object = as.formula(paste0(resp,"~",expr)))
    ress=c(stri_flatten(round(res[,2],digits = 4),collapse = ""),stri_flatten(res[,2],collapse = ""),posteriors[i,2],expr)
    if(!((ress[1] %in% values(rhash))))
      rhash[[ress[1]]]=ress
    else
    {
      if(ress[1] %in% keys(rhash))
      {
        rhash[[ress[1]]][3]= (as.numeric(rhash[[ress[1]]][3]) + as.numeric(ress[3]))
        if(stri_length(rhash[[ress[1]]][4])>stri_length(expr))
          rhash[[ress[1]]][4]=expr
      }
    }

  }
  res=as.data.frame(t(values(rhash)[c(3,4),]))
  res$V1=as.numeric(as.character(res$V1))
  res=res[which(res$V1>thf),]
  res=res[order(res$V1, decreasing = T),]
  clear(rhash)
  rm(rhash)
  res[which(res[,1]>1),1]=1
  colnames(res)=c("posterior","tree")
  return(res)
}


#a set of nonlinearities that will be used in the DBRM model
cosi=function(x)cos(x/180*pi)
sini=function(x)sin(x/180*pi)
expi=function(x)
{
  r=exp(x)
  if(r==Inf)
    return(10000000)
  else
    return(r)
}

InvX=function(x)
{
  if(x==0)
    return(10000000)
  else
    return(1/x)

}
troot=function(x)abs(x)^(1/3)

#specify the number of runs
MM = 100
#specify the number of threads used
M = 12
NM= 1000
#specify the number of features + 1 per model
compmax = 16
#specify some preliminary filtration of the features' tresholds
th=(10)^(-5)
thf=0.05

j=1
tryCatch({

set.seed(j)

#read the epigenetic data in and preprocess it
data.example = read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/Epigenetic%20Data/epigen.txt"),sep = ",",header = T)[,2:30]
data.example=data.example[sample.int(dim(data.example)[1],500),]
data.example=data.example[,c(2,5,6,8:10,12:17,21,23,24,29)]
data.example$eg3000=data.example$express_noisy>3000
data.example$eg10000=data.example$express_noisy>10000
data.example$express_noisy=NULL
data.example$pos1 = data.example$pos
data.example$pos2 = data.example$pos
data.example$pos3 = data.example$pos
fparams =colnames(data.example )[-c(1,2,3,18,19,20)]#c(colnames(data.example )[c(8:10,12:17,21:24,29)],"f(data.example$pos,model=\"ar1\")","f(data.example$pos1,model=\"rw1\")","f(data.example$pos2,model=\"iid\")","f(data.example$pos3,model=\"ou\")")
fobservs = colnames(data.example)[2]

#specify the initial formula
formula1 = as.formula(paste(fobservs,"~ 1 +",paste0(fparams,collapse = "+")))
#specify tuning parameters of the algorithm for exploring DBRM of interest
#notice that allow_offsprings=3 corresponds to the GMJMCMC runs and
#allow_offsprings=4 -to the RGMJMCMC runs
vect=list(formula = formula1,outgraphs=F,max.time = 180,data = data.example,latnames = c("f(data.example$pos,model=\"ar1\")","f(data.example$pos1,model=\"rw1\")","f(data.example$pos2,model=\"iid\")","f(data.example$pos3,model=\"ou\")","offset(log(total_bases))"),estimator = estimate.inla.poisson,estimator.args =  list(data = data.example),recalc_margin = 199, save.beta = F,interact = T,relations=c("cos","sigmoid","tanh","atan","sin","erf"),relations.prob =c(0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=3,mutation_rate = 200, last.mutation = 2000,max.tree.size = 200000, Nvars.max = (compmax-1),p.allow.replace=0.7,p.allow.tree=0.1,p.nor=0.3,p.and = 0.7),n.models = 10000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,print.freq = 1000,advanced.param = list(
  max.N.glob=as.integer(10),
  min.N.glob=as.integer(5),
  max.N=as.integer(3),
  min.N=as.integer(1),
  printable = F))


params = list(vect)[rep(1,M)]
#specify additional information to be used on different threads (e.g. seeds)
for(i in 1:M)
{

  params[[i]]$cpu=i*j
  params[[i]]$simul="scenario_epi_"
  params[[i]]$simid=j
  params[[i]]$NM=1000
  params[[i]]$simlen=23
}
gc()
#perform garbage collection
gc()
print(paste0("begin simulation ",j))
#explore DBRM on M threads in parallel using the GMJMCMC/RGMJMCMC algorithm
results=parall.gmj(X = params, M = 3)

#print the raw results from the threads
print(results)

#prepare the data structures for final analysis of the runs
resa=array(data = 0,dim = c(compmax,M*3))
post.popul = array(0,M)
max.popul = array(0,M)
nulls=NULL

#check which threads had non-zero exit status
not.null=1
for(k in 1:M)
{
  if(is.character(results[[k]]))
  {
    nulls=c(nulls,k)
    next
  }
  if(length(results[[k]])==0)
  {
    nulls=c(nulls,k)
    next
  }
  else
  {
    not.null = k
  }

}

#for all of the successful runs collect the results into the data structures
for(k in 1:M)
{
  if(k %in% nulls)
  {
    results[[k]]=results[[not.null]]
  }
  max.popul[k]=results[[k]]$cterm
  post.popul[k]=results[[k]]$post.populi
  if(length(resa[,k*3-2])==(length(results[[k]]$fparam)+1))
  {
    resa[,k*3-2]=c(results[[k]]$fparam,"Post.Gen.Max")
    resa[,k*3-1]=c(results[[k]]$p.post,results[[k]]$cterm)
    resa[,k*3]=rep(post.popul[k],length(results[[k]]$p.post)+1)
  }else
  {
    #idsx=order(results[[k]]$p.post,decreasing = T,na.last = T)
    resa[,k*3-2]=rep(results[[k]]$fparam[1],length(resa[,k*3-2]))
    resa[,k*3-1]=rep(0,length(resa[,k*3-1]))
    resa[,k*3]=rep(-10^9,length(resa[,k*3]))
    max.popul[k]= -10^9
    post.popul[k]= -10^9
  }

}

#delete the unused further variables and perfrom garbage collection
gc()
rm(results)

#renormalize estimates of the marginal inclusion probabilities
#based on all of the runs
ml.max=max(max.popul)
post.popul=post.popul*exp(-ml.max+max.popul)
p.gen.post=post.popul/sum(post.popul)
hfinal=hash()
for(ii in 1:M)
{
  resa[,ii*3]=p.gen.post[ii]*as.numeric(resa[,ii*3-1])
  resa[length(resa[,ii*3]),ii*3]=p.gen.post[ii]
  if(p.gen.post[ii]>0)
  {
    for(jj in 1:(length(resa[,ii*3])-1))
    {
      if(resa[jj,ii*3]>0)
      {
        if(as.integer(has.key(hash = hfinal,key =resa[jj,ii*3-2]))==0)
          hfinal[[resa[jj,ii*3-2]]]=as.numeric(resa[jj,ii*3])
        else
          hfinal[[resa[jj,ii*3-2]]]=hfinal[[resa[jj,ii*3-2]]]+as.numeric(resa[jj,ii*3])
      }

    }
  }
}

posteriors=values(hfinal)

#delete the unused variables
clear(hfinal)
rm(hfinal)
rm(resa)
rm(post.popul)
rm(max.popul)
posteriors=as.data.frame(posteriors)
posteriors=data.frame(X=row.names(posteriors),x=posteriors$posteriors)
posteriors$X=as.character(posteriors$X)
#simplify the found features and their posteriors
tryCatch({
  res1=simplifyposteriors(X = data.example,posteriors = posteriors, th,thf,resp = "methylated_bases")
  row.names(res1)=1:dim(res1)[1]
  write.csv(x =res1,row.names = F,file = paste0("postEPIGEN_",j,".csv"))
},error = function(err){
  print("error")
  print(err)
  write.csv(x =posteriors,row.names = F,file = paste0("postEPIGENERR_",j,".csv"))
},finally = {
})
#delete the unused further variables and perfrom garbage collection
rm(X)
rm(data.example)
rm(vect)
rm(params)
gc()
print(paste0("end simulation ",j))
},error = function(err){
print("error")
j=j-1
print(paste0("repeat  simulation ",j))
},finally = {

print(paste0("end simulation ",j))
#delete the unused further variables and perfrom garbage collection
rm(X4)
rm(data.example)
rm(vect)
rm(params)
gc()
})


