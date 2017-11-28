
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")
library("survival")

load("/storage/edelmand/IMPACT2/Daten/betas_clinic_for_analysis.RData")

clinical<-cbind(clinic2[,-1], t(betas_ges))


if(TRUE){
  #function to be used in screening
  coxscreen<-function(x) tryCatch(return(summary(coxph(Surv(V1,V2) ~1+V3,data = as.data.frame(cbind(clinical$fudays_35fu,clinical$death_event_35fu
  ,clinical[[x]]))))$logtest[3]),error = function(err) {print(err)
    return(1)})
  
  screen<-mclapply(FUN = coxscreen, X = names(clinical)[-c(14:18)],mc.cores = 63)
}else
{
 screen<-read.csv(file = "pvalues.csv",header = F,row.names = F,col.names = F) 
}



ids<-NULL
for(i in 1:length(screen))
  if(screen[[i]]<=0.000001)
    ids<-c(ids,i)

#ids=read.csv("indeces.csv")

clinical<-clinical[,c(ids,14:18)]

clinical<-clinical[-which(clinic2$stage4==1),]
clinical<-clinical[-which(clinic2$neoadth==1),]

#some clustering
sel <- which(anno$chr=="chr9")
clinical <- clinical[,c(1:18,(sel+18))]
anno <- anno[sel,]
sel1 <- order(anno$pos)
clinical[,18:dim(clinical)[2]]<- clinical[,sel+18]
anno <- anno[sel,]
corri <- cor((clinical[,19:dim(clinical)[2]]),method="spearman")
disti <- as.dist(1-abs(corri))
clusti <- hclust(disti,method="complete")
clustering <- cutree(clusti,h=0.85)

data.example <- clinical[,1:18]

idclust<-unique(clustering)
for(i in idclust)
{
  avg.id <- which(clustering==i)
  if(length(avg.id)>1)
    data.example[[paste0("clust",i)]]<- rowMeans(clinical[,18+avg.id])
  else if(length(avg.id)==1)
    data.example[[paste0("clust",i)]]<- clinical[,18+avg.id]
}



#options("expressions"=500000)
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")
library("survival")
#fla= paste0("Surv(fudays_35fu,crc_death_35fu)~1 + ",paste(names[-c(16680:16684)],collapse = "+"))

formula=as.formula(paste0("Surv(fudays_35fu,crc_death_35fu)~1 + ",paste(names(clinical)[-c((length(clinical)-4):length(clinical))],collapse = "+")))


#estimator function (for speed inla is not adressed here currently)

#estimator function (for speed inla is not adressed here currently)
estimate.lm.BIC <- function(formula, data, n = 1125, m = 1000, c = 16,u=150)
{
  size<-stri_count_fixed(str = as.character(formula)[3],pattern = "+")
  
  if(size>u||size ==0)
  {
    return(list(mlik = (-50000 + rnorm(1,0,1)),waic = 50000 + rnorm(1,0,1), dic =  50000+ rnorm(1,0,1),summary.fixed =list(mean = array(0,size+1))))
  }
  out<-NULL
  
  capture.output(tryCatch({out <- (coxph(formula, data=data))}))
  
  if(is.null(out))
  {
    return(list(mlik = (-50000 + rnorm(1,0,1)),waic = 50000 + rnorm(1,0,1), dic =  50000+ rnorm(1,0,1),summary.fixed =list(mean = array(0,size+1))))
  }
  
  
  logmarglik <- (2*out$loglik[1]- size*log(n))/2# for mbic would be out$loglik[1] -  size*log(m*m*n/c) ) + 2*log(factorial(size)))/2#
  return(list(mlik = logmarglik,waic = -logmarglik , dic = -logmarglik, summary.fixed = c(0,out$coefficients)))
}


estimate.lm.MBIC2 <- function(formula, data, n = dim(data.example)[1], m =dim(data.example)[2], c = 16,u=300)
{
  size<-stri_count_fixed(str = as.character(formula)[3],pattern = "+")
  
  if(size>u)
  {
    return(list(mlik = (-50000 + rnorm(1,0,1) - size*log(m*m*n/c) + 2*log(factorial(size+1))),waic = 50000 + rnorm(1,0,1), dic =  50000+ rnorm(1,0,1),summary.fixed =list(mean = array(0,size+1))))
  }else{
    out<-NULL
    
    capture.output(tryCatch({out <- (coxph(formula, data=data))}))
    
    if(is.null(out))
    {
      return(list(mlik = (-50000 + rnorm(1,0,1)),waic = 50000 + rnorm(1,0,1), dic =  50000+ rnorm(1,0,1),summary.fixed =list(mean = array(0,size+1))))
    }
    logmarglik <- (2*out$loglik[1]- size*log(n) - size*log(m*m/c) + 2*log(factorial(size)))/2
    # use dic and aic as bic and aic correspondinly
    return(list(mlik = logmarglik,waic =  -logmarglik , dic =   -logmarglik,summary.fixed =list(mean = c(0,out$coefficients))))
  }
}



do.call.emjmcmc<-function(vect)
{
  library("survival")  
  set.seed(as.integer(vect$cpu))
  do.call(runemjmcmc, vect[1:vect$simlen])
  vals<-values(hashStat)
  fparam<-mySearch$fparam
  cterm<-max(vals[1,],na.rm = T)
  ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
  post.populi<-sum(exp(values(hashStat)[1,][1:vect$NM]-cterm),na.rm = T)
  clear(hashStat)
  rm(hashStat)
  rm(vals)
  gc()
  return(list(post.populi = post.populi, p.post =  ppp$p.post, cterm = cterm, fparam = fparam))
}

#data.example = clinical
rm(clinical)
gc()


formula=as.formula(paste0("Surv(fudays_35fu,crc_death_35fu)~1 + ",paste(names(data.example)[-c(14:18)],collapse = "+")))



MM = 63
M = 63
size.init=1000
NM= 5000
compmax = 101
th<-(10)^(-5)
thf<-0.001
gc()
M.cpu<-63

vect<-list(formula = formula, locstop.nd = T, keep.origin = F,p.add = 0.1,max.time = 75, p.add.default = 0.1, pool.cor.prob = T,secondary <- NULL, outgraphs=F,data = data.example,estimator = estimate.lm.MBIC2,presearch=F, locstop =F,estimator.args =  list(data = data.example),recalc_margin = 999,gen.prob = c(1,0,0,0,0), save.beta = F,interact = F,relations=c("cos"),relations.prob =c(0.1),interact.param=list(allow_offsprings=3,mutation_rate = 1000, last.mutation = 15000, max.tree.size = 4, Nvars.max =(compmax-1),p.allow.replace=0.7,p.allow.tree=0.25,p.nor=0,p.and = 0.9),n.models = 25000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,print.freq = 1000,advanced.param = list(
  max.N.glob=as.integer(130),
  min.N.glob=as.integer(10),
  max.N=as.integer(5),
  min.N=as.integer(1),
  printable = F))

params <- list(vect)[rep(1,M)]

for(i in 1:M)
{
  #cov.names<-names[sample.int(n = length(names),size = size.init,prob = abs(cors))]
  #sum<-summary(lm(as.formula(paste0("Y~1+",paste(cov.names,collapse = "+"))),data = geno))
  #cov.names<-names(sum$coefficients[-1,4])
  #params[[i]]$
  params[[i]]$cpu<-i
  params[[i]]$simul<-"scenario_JM_"
  params[[i]]$simid<-1
  params[[i]]$NM<-NM
  params[[i]]$simlen<-31
}


gc()
print(paste0("begin simulation ",1))
results<-parall.gmj(X = params, M = M.cpu)


compmax=length(results[[3]]$fparam)+1

#res<-do.call(runemjmcmc,args = params[[3]][1:27])
#res$p.post
#length(which(!is.na(res$m.post)))
#detected<-mySearch$fparam[which(res$p.post>0.1)]


#print(results)

#wait()

resa<-array(data = 0,dim = c(compmax,M*3))
post.popul <- array(0,M)
max.popul <- array(0,M)
nulls<-NULL

not.null<-1
for(k in 1:M)
{
  if(is.character(results[[k]]))
  {
    nulls<-c(nulls,k)
    next
  }
  if(length(results[[k]])==0)
  {
    nulls<-c(nulls,k)
    next
  }
  else
  {
    not.null <- k
  }
  
}


for(k in 1:M)
{
  if(k %in% nulls)
  {
    results[[k]]<-results[[not.null]]
  }
  max.popul[k]<-results[[k]]$cterm
  post.popul[k]<-results[[k]]$post.populi
  if(length(resa[,k*3-2])==(length(results[[k]]$fparam)+1))
  {
    resa[,k*3-2]<-c(results[[k]]$fparam,"Post.Gen.Max")
    resa[,k*3-1]<-c(results[[k]]$p.post,results[[k]]$cterm)
    resa[,k*3]<-rep(post.popul[k],length(results[[k]]$p.post)+1)
  }else
  {
    #idsx<-order(results[[k]]$p.post,decreasing = T,na.last = T)
    resa[,k*3-2]<-rep(results[[k]]$fparam[1],length(resa[,k*3-2]))
    resa[,k*3-1]<-rep(0,length(resa[,k*3-1]))
    resa[,k*3]<-rep(-10^9,length(resa[,k*3]))
    max.popul[k]<- -10^9
    post.popul[k]<- -10^9
  }
  
}


gc()

ml.max<-max(max.popul)
post.popul<-post.popul*exp(-ml.max+max.popul)
p.gen.post<-post.popul/sum(post.popul)
hfinal<-hash()
for(ii in 1:M)
{
  resa[,ii*3]<-p.gen.post[ii]*as.numeric(resa[,ii*3-1])
  resa[length(resa[,ii*3]),ii*3]<-p.gen.post[ii]
  if(p.gen.post[ii]>0)
  {
    for(jj in 1:(length(resa[,ii*3])-1))
    {
      if(resa[jj,ii*3]>0)
      {
        #print(paste0(ii,"  and ",jj))
        if(as.integer(has.key(hash = hfinal,key =resa[jj,ii*3-2]))==0)
          hfinal[[resa[jj,ii*3-2]]]<-as.numeric(resa[jj,ii*3])
        else
          hfinal[[resa[jj,ii*3-2]]]<-hfinal[[resa[jj,ii*3-2]]]+as.numeric(resa[jj,ii*3])
      }
      
    }
  }
}

posteriors<-values(hfinal)

#print(posteriors)
clear(hfinal)
rm(hfinal)
rm(resa)
rm(post.popul)
rm(max.popul)
posteriors<-as.data.frame(posteriors)
posteriors<-data.frame(X=row.names(posteriors),x=posteriors$posteriors)
posteriors$X<-as.character(posteriors$X)


write.csv(x =posteriors,row.names = F,file = paste0("postFull_",9,".csv"))


posteriors <- read.csv("postFull_1.csv",stringsAsFactors = F)
X = data.example
nms<-names(X)
names(X)=paste0(rep("I(",length(names(X))),names(X),rep(")",length(names(X))))
X=X[,which(names(X)%in%posteriors[,1])]
nam<-names(X)
nms<-nms[which(names(X)%in%posteriors[,1])]
names(X)<-nms
X[,2]<-as.numeric(as.character(X[,2]))
X[,3]<-as.numeric(as.character(X[,3]))
cors<-cor(X)
rhash<-hash()
ready<-NULL
for(i in 1:length(posteriors[,1]))
{
  if(i %in% ready)
    next
  else
  {
    id<-which(nam==posteriors[,1])
    ids<-which(posteriors[,1]%in%rownames(cors)[which(abs(cors[,id])>0.2)])
    #print(ids)
    ready<-c(ready,i,which(abs(cors[,id])>0.2))
    
    posteriors[i,2]<-sum(posteriors[ids,2])
    print( posteriors[i,2])
    expr<-posteriors[i,1]
    
    #res<-model.matrix(data=X,object = as.formula(paste0(resp,"~",expr)))
    ress<-c(stri_flatten(round(rnorm(1,0,1),digits = 4),collapse = ""),stri_flatten(rnorm(1,0,1),collapse = ""),posteriors[i,2],expr)
    if(!((ress[1] %in% values(rhash))))
      rhash[[ress[1]]]<-ress
    else
    {
      if(ress[1] %in% keys(rhash))
      {
        rhash[[ress[1]]][3]<- (as.numeric(rhash[[ress[1]]][3]) + as.numeric(ress[3]))
        if(stri_length(rhash[[ress[1]]][4])>stri_length(expr))
          rhash[[ress[1]]][4]<-expr
      }
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



write.csv(x =res,row.names = F,file = paste0("postAnal_","survival",".csv"),sep = ';' )



#posteriors <- read.csv("postFull_1.csv",stringsAsFactors = F)
