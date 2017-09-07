#library(dplyr)

library(plyr)
library(dplyr)
library(magrittr)

options("exppression" = 500000)
##########################################
#: Calculate True h^2 for each Scenario :#
##########################################
simPars <- 
  read.table("/home/michaelh/SIMULATION_paper/scripts/simulationParameters_new.txt", 
	     header = TRUE, 
	     stringsAsFactors = FALSE) %>%
  select(SNPId = SNP_causal, Pos = BP_causal, MAF = MAF_causal, S1, S2, S3, S4)

causSNPsS1 <- with(simPars, SNPId[S1 != 0])
causSNPsS2 <- with(simPars, SNPId[S2 != 0]) 
causSNPsS3 <- with(simPars, SNPId[S3 != 0]) 
causSNPsS4 <- with(simPars, SNPId[S4 != 0]) 

genoData <- read.table("/home/michaelh/SIMULATION_paper/data/CHR1_NFBC.raw",
                       header = TRUE,
                       stringsAsFactors = FALSE)
names(genoData) <- gsub("_.$", "", names(genoData)) %>% gsub("\\.", "-", .)

genoData <- genoData[c("FID", "IID", simPars$SNPId)]

expectedData <- data.frame(as.matrix(genoData[-(1:2)]) %*% as.matrix(simPars[c('S1', 'S2', 'S3', 'S4')])) 

explainedVar <- apply(expectedData, 2, function(x) var(x)/(var(x)+1)) # True h^2


############################################
#: Calculate Neigbourhoods of causal SNPs :#
#: Definition: Neigbourhood 		  :#
#	       corr > .3 and +/- 1(2) MBP :#
############################################
distance_cutoff <- 1e6
corr_cutoff <- .3

physMap <- read.table("/home/michaelh/SIMULATION_paper/data/CHR1_NFBC.bim",
		      header = FALSE,
		      stringsAsFactors = FALSE)[c(2, 4)] %>%
	   select(SNPid = V2, pos = V4) %>%
           filter(!grepl('^cnv', SNPid))	   
	       		  
genoData <- read.table("/home/michaelh/SIMULATION_paper/data/CHR1_NFBC.raw", 
			header = TRUE, 
			stringsAsFactors = FALSE)
names(genoData) <- gsub("_.$", "", names(genoData)) %>% gsub("\\.", "-", .)	 

findNeigbour <- function(causSNP) {
  #causSNP <- "rs1498308"
  #print(causSNP)
    physMap %>%
    mutate(distance = abs(pos - pos[SNPid == causSNP])) %>%
    filter(distance <= distance_cutoff) %>%
    mutate(corr = sapply(SNPid, function(x) abs(cor(genoData[causSNP], genoData[x]))),
	   causSNPid = causSNP) %>%
    filter(corr >= corr_cutoff)
}

dataNeigbourhoodS1 <- plyr::ldply(causSNPsS1, findNeigbour)
dataNeigbourhoodS2 <- plyr::ldply(causSNPsS2, findNeigbour)
dataNeigbourhoodS3 <- plyr::ldply(causSNPsS3, findNeigbour)
dataNeigbourhoodS4 <- plyr::ldply(causSNPsS4, findNeigbour)

source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")

pheno<-read.csv(paste0("/home/michaelh/SIMULATION_paper/data_S2_nocausal_5402/pimass/data.recode.pheno_",1,".txt"),header = F)
geno<-t(read.csv("/home/michaelh/SIMULATION_paper/data_S2_nocausal_5402/pimass/data.recode.mean.geno.txt",header = F,stringsAsFactors = F))
names<-geno[1,]
geno<-as.data.frame(geno[-c(1,2,3),])
names(geno)<-names
geno$Y<-pheno$V1
geno <- as.data.frame(mclapply(geno, as.numeric))

cors<-cor(geno$Y,geno[,1:24602])

estimate.lm.MBIC2 <- function(formula, data, n = 5402, m = 24602, c = 4,u=150)
{
  size<-stri_count_fixed(str = as.character(formula)[3],pattern = "+")
  
  if(size>u)
  {
    return(list(mlik = (-50000 + rnorm(1,0,1) - size*log(m*m*n/c) + 2*log(factorial(out$rank))),waic = 50000 + rnorm(1,0,1), dic =  50000+ rnorm(1,0,1),summary.fixed =list(mean = array(0,size+1))))
  }
  
  out <- lm(formula = formula,data = data)
  logmarglik <- (2*logLik(out) - out$rank*log(m*m*n/c) + 2*log(factorial(out$rank)))/2
  # use dic and aic as bic and aic correspondinly
  return(list(mlik = logmarglik,waic = AIC(out) , dic =  BIC(out),summary.fixed =list(mean = coef(out))))
  
}


MM = 10
M = 16
size.init=100
NM= 1000
compmax = 41
th<-(10)^(-5)
thf<-0.05

for(j in 1:100)
{
  tryCatch({
    
    set.seed(j)
    
    pheno<-read.csv(paste0("/home/michaelh/SIMULATION_paper/data_S2_nocausal_5402/pimass/data.recode.pheno_",j,".txt"),header = F)
    geno$Y<-as.numeric(pheno$V1)
                    
    
    cov.names<-names[sample.int(n = length(names),size = size.init,prob = abs(cors))]
    sum<-summary(lm(as.formula(paste0("Y~1+",paste(cov.names,collapse = "+"))),data = geno))
    cov.names<-names(sum$coefficients[-1,4])
    
    formula1 <- as.formula(paste0("Y~1+",paste(cov.names,collapse = "+")))
    
    secondary <-names[-which(names %in% cov.names)] 
    
  
    vect<-list(outgraphs=F,data = geno,estimator = estimate.lm.MBIC2,presearch=F, locstop =T,estimator.args =  list(data = geno),recalc_margin = 249,gen.prob = c(1,0,0,0,0), save.beta = F,interact = T,relations=c("cos"),relations.prob =c(0.1),interact.param=list(allow_offsprings=3,mutation_rate = 250,last.mutation = 1000, max.tree.size = 4, Nvars.max =40,p.allow.replace=0.7,p.allow.tree=0.2,p.nor=0,p.and = 0.9),n.models = 10000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,print.freq = 100,advanced.param = list(
      max.N.glob=as.integer(10),
      min.N.glob=as.integer(5),
      max.N=as.integer(3),
      min.N=as.integer(1),
      printable = F))
    
    #res<-do.call(runemjmcmc,args = vect)
    
    params <- list(vect)[rep(1,M)]
    
    for(i in 1:M)
    {
      cov.names<-names[sample.int(n = length(names),size = size.init,prob = abs(cors))]
      sum<-summary(lm(as.formula(paste0("Y~1+",paste(cov.names,collapse = "+"))),data = geno))
      cov.names<-names(sum$coefficients[-1,4])
      
      formula1 <- as.formula(paste0("Y~1+",paste(cov.names,collapse = "+")))
      
      params[[i]]$secondary <-names[-which(names %in% cov.names)] 
      params[[i]]$formula = formula1
      params[[i]]$cpu<-i*j
      params[[i]]$simul<-"scenario_JM_"
      params[[i]]$simid<-j
      params[[i]]$simlen<-25
    }
    
    gc()
    print(paste0("begin simulation ",j))
    results<-parall.gmj(X = params)
    print(results)
    
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
        resa[,k*3-2]<-rep(results[[k]]$fparam[1],length(resa[,k*3-2]))
        resa[,k*3-1]<-rep(0,length(resa[,k*3-1]))
        resa[,k*3]<-rep(-10^9,length(resa[,k*3]))
      }
      
    }
    
    
    gc()
    rm(results)
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
    
    print(posteriors)
    clear(hfinal)
    rm(hfinal)
    rm(resa)
    rm(post.popul)
    rm(max.popul)
    posteriors<-as.data.frame(posteriors)
    posteriors<-data.frame(X=row.names(posteriors),x=posteriors$posteriors)
    posteriors$X<-as.character(posteriors$X)
    tryCatch({
      res1<-simplifyposteriors(X = X4,posteriors = posteriors, th,thf)
      row.names(res1)<-1:dim(res1)[1]
      write.csv(x =res1,row.names = F,file = paste0("postGMJSIM_",j,".csv"))
    },error = function(err){
      print("error")
      write.csv(x =posteriors,row.names = F,file = paste0("postGMJSIM_",j,".csv"))
    },finally = {
      
      print(paste0("end simulation ",j))
      
    })
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
    rm(X4)
    rm(data.example)
    rm(vect)
    rm(params)
    gc()
  })
  
}







detected<-sort(mySearch$fparam[which(res$p.post>0.5)])
detected<-stri_replace(str = detected,fixed = "I(",replacement = "")
detected<-stri_replace(str = detected,fixed = ")",replacement = "")

detect.true.unique<-unique(dataNeigbourhoodS2$causSNPid[which(dataNeigbourhoodS2$SNPid %in% detected)])
detect.true<-which(detected %in% dataNeigbourhoodS2$SNPid)

detlen<-length(detect.true.unique)
totlen<-length(detected)-length(detect.true)+length(detect.true.unique)

detlen/20
(totlen-detlen)/totlen



