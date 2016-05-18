rm(list = ls(all = TRUE))

install.packages("INLA", repos="http://www.math.ntnu.no/inla/R/testing")
install.packages("bigmemory")
install.packages("snow")
install.packages("Rmpi")
install.packages("ade4")
install.packages("sp")
install.packages("BAS")


library(sp)
library(INLA)
library(parallel)
library(bigmemory)
library(snow)
library(MASS)
library(ade4)
library(copula)
library(compiler)
library(BAS)
require(stats)
#compile INLA 

#closeAllConnections()

#inla.compiled <- cmpfun(inla)

estimate.bas.glm <- function(formula, data, Y, family, prior, logn)
{
  
  #only poisson and binomial families are currently adopted
  X <- model.matrix(object = formula,data = data)
  out <- bayesglm.fit(X, Y, family=family,coefprior=prior)
  # use dic and aic as bic and aic correspondinly
  return(list(mlik = out$logmarglik,waic = -(out$deviance + 2*out$rank) , dic =  -(out$deviance + logn*out$rank)))
  
}

parallelize<-function(X,FUN)
{
  max.cpu <- length(X)
  cl <-makeCluster(max.cpu ,type = paral.type,outfile = "")#outfile = ""
  clusterEvalQ(cl = cl,expr = c(library(INLA),library(bigmemory)))
  if(exists("statistics"))
  {
    clusterExport(cl=cl, "statistics")
    clusterEvalQ(cl,{statistics <- attach.resource(statistics);1})
  }
  res.par <- parLapply(cl = cl, X, FUN)
  stopCluster(cl)
  return(res.par)
}

estimate.bas.lm <- function(formula, data, prior, n, g = 0)
{
  
  out <- lm(formula = formula,data = data)
  # 1 for aic, 2 bic prior, else g.prior
  
  p <- out$rank
  if(prior == 1)
  {
    ss<-sum(out$residuals^2)
    logmarglik <- -0.5*(log(ss)+2*p) 
  }
  else if(prior ==2)
  {
    ss<-sum(out$residuals^2)
    logmarglik <- -0.5*(log(ss)+log(n)*p) 
  }
  else
  {
    Rsquare <- summary(out)$r.squared
    logmarglik =  .5*(log(1.0 + g) * (n - p -1)  - log(1.0 + g * (1.0 - Rsquare)) * (n - 1))*(p!=1)
  }
  
  # use dic and aic as bic and aic correspondinly
  return(list(mlik = logmarglik,waic = AIC(out) , dic =  BIC(out)))
  
}

# covariates are thus
# 1) A constant +1
# 2) Base structures (CHG, CGH, CHH) +2
# 3) Gene class (5 classes) +4
# 4) Position within a gene (4 positions including promoters, after coding parts, intrones, exones) +3 
# 5) Gene expression (continious) +1
# 6) A CpG island factor +1
# 7) A tranposone factor +1
# 8) combinations (if reasonable)
#closeAllConnections()

remove(EMJMCMC2016)
remove(mySearch)

EMJMCMC2016 <- setRefClass(Class = "EMJMCMC2016", inheritPackage = TRUE , 
  fields = list(estimator.args = "list",
                max.cpu = "integer",
                objective = "integer",
                p.prior = "numeric",
                min.N = "integer",
                estimator = "function",
                parallelize = "function",
                parallelize.global = "function",
                parallelize.hyper  = "function",
                min.N.randomize = "integer",
                max.N.randomize = "integer",
                type.randomize = "integer", 
                aa = "numeric", 
                cc = "numeric",
                printable.opt = "logical",
                fobserved = "vector",
                switch.type = "integer",
                n.size ="integer",
                LocImprove = "array",
                max.N = "integer",
                recalc.margin = "integer",
                max.N.glob = "integer",
                min.N.glob = "integer",
                max.cpu.glob = "integer",
                max.cpu.hyper = "integer",
                switch.type.glob = "integer",
                isobsbinary = "array",
                fparam = "vector",
                p.add = "array",
                latent.formula = "character",
                Nvars = "integer", 
                seed = "integer", 
                M.nd = "integer",
                locstop.nd = "logical",
                M.mcmc = "integer",
                SA.param = "list",
                g.results = "big.matrix"
  ), 
  methods = list(
  #class constructor   
  initialize = function(estimator.function = inla, estimator.args.list = list(family = "gaussian",data = data.example, control.fixed=list(prec=list(default= 0.00001),prec.intercept = 0.00001,mean=list(default= 0),mean.intercept = 0),
                                                                              control.family = list(hyper = list(prec = list(prior = "loggamma",param = c(0.00001,0.00001),initial = 0))),
                                                                              control.compute = list(dic = TRUE, waic = TRUE, mlik = TRUE)), search.args.list = NULL,latent.formula = "")
  {
    estimator <<- estimator.function
    estimator.args <<- estimator.args.list
    latent.formula <<- latent.formula   
    g.results <<- big.matrix(nrow = 4,ncol = 2)
    g.results[1,1]<- -Inf
    g.results[1,2]<- 1
    g.results[2,1]<- Inf
    g.results[2,2]<- 1
    g.results[3,1]<- Inf
    g.results[3,2]<- 1
    g.results[4,1]<- 0
    g.results[4,2]<- 0
    
    if(is.null(search.args.list))
    {
      max.cpu <<- as.integer(2)
      objective <<- as.integer(1)
      if(Sys.info()['sysname']=="Windows")
      {
        parallelize <<- lapply
        parallelize.global <<- lapply
        parallelize.hyper  <<- lapply
      }
      else
      {
        parallelize <<- mclapply
        parallelize.global <<- mclapply
        parallelize.hyper  <<- mclapply
      }
      Nvars  <<- as.integer(length(fparam.example)) 
      min.N <<- as.integer(Nvars/6)
      min.N.glob <<- as.integer(Nvars/3)
      max.N.glob <<- as.integer(Nvars/2)
      max.N <<- as.integer(Nvars/5)
      switch.type.glob <<- as.integer(2)
      min.N.randomize <<- as.integer(1)
      max.N.randomize <<- as.integer(1)
      type.randomize <<- as.integer(3)
      max.cpu.glob <<- as.integer(5)
      max.cpu.hyper <<- as.integer(2)
      printable.opt <<- FALSE
      locstop.nd <<-FALSE
      aa <<- 0.9
      cc <<- 0.0
      M.nd <<- as.integer(Nvars)
      M.mcmc <<- as.integer(5)
      SA.param <<- list(t.min = 0.0001,t.init = 10, dt = 3, M = 3)
      fobserved <<- fobserved.example
      switch.type <<- as.integer(2)
      n.size <<- as.integer(10)
      LocImprove <<- as.array(c(50,50,50,50,150))
      isobsbinary <<- as.array(0:(length(fparam.example)-1))
      fparam <<- fparam.example
      p.add <<- array(data = 0.5,dim = Nvars)
      if(exists("statistics"))
      {
        recalc.margin <<- as.integer(100)
      }else if(exists("statistics1"))
      {
        recalc.margin <<- as.integer(100)
      }else
      {
        recalc.margin <<- as.integer(-1)
      }
      seed <<- as.integer(runif(n = 1,min = 1,max = 10000))
      p.prior <<- runif(n = Nvars, min = 0.5,max = 0.5)
    }
    else
    {
      max.cpu <<- as.integer(search.args.list$max.cpu)
      objective <<- as.integer(search.args.list$objective)
      parallelize <<- search.args.list$parallelize
      parallelize.global <<- search.args.list$parallelize.global
      parallelize.hyper  <<- search.args.list$parallelize.hyper
      p.prior <<-  search.args.list$p.prior
      min.N <<- as.integer(search.args.list$min.N)
      printable.opt <<- search.args.list$printable.opt
      min.N.glob <<- as.integer(search.args.list$min.N.glob)
      max.N.glob <<- as.integer(search.args.list$max.N.glob)
      switch.type.glob <<- as.integer(search.args.list$switch.type.glob)
      min.N.randomize <<- as.integer(search.args.list$min.N.randomize)
      max.N.randomize <<- as.integer(search.args.list$max.N.randomize)
      type.randomize <<- as.integer(search.args.list$type.randomize)
      max.cpu.glob <<- as.integer(search.args.list$max.cpu.glob)
      locstop.nd <<-search.args.list$locstop.nd
      max.cpu.hyper <<- as.integer(search.args.list$max.cpu.hyper)
      aa <<- search.args.list$lambda.a
      cc <<- search.args.list$lambda.c
      M.nd <<- as.integer(search.args.list$stepsGreedy)
      M.mcmc <<- as.integer(search.args.list$stepsLocMCMC)
      SA.param <<- search.args.list$SA.params
      fobserved <<- search.args.list$fobserved
      switch.type <<- as.integer(search.args.list$fswitch.type)
      n.size <<- as.integer(search.args.list$n.size)
      LocImprove <<-search.args.list$prior.optimizer.freq
      max.N <<- as.integer(search.args.list$max.N)
      fparam <<- search.args.list$fparam
      isobsbinary <<- as.array(0:(length(fparam)-1))
      p.add <<- as.array(search.args.list$p.add)
      recalc.margin <<- as.integer(search.args.list$recalc.margin)
      Nvars  <<- as.integer(length(fparam)) 
      seed <<-  search.args.list$seed
    }
    
  },
  #transform binary numbers to decimal
  bittodec = function(bit) #transform a binary vector into a natural number to correspond between vector of solutions and storage array
  {
    n<-length(bit)
    dec <- 0
    for(i in 1:n)
    {
      j<-n-i
      dec <- dec + ((2)^j)*bit[i] 
    }
    return(dec)
  },
  #transform decimal numbers to binary
  dectobit = function(dec) #transform a natural number into a binary vector to correspond between vector of solutions and storage array
  {
    q<-dec
    bin<-NULL
    while(q!=0)
    {
      r<-q/2
      q=floor(r)
      bin<-c(as.integer(2*(r-q)),bin)
    }
    return(bin)
  },
  #calculate move probabilities
  calculate.move.logprobabilities = function(varold, varnew, switch.type,min.N,max.N)
  {
    if(switch.type == 1) # random size random N(x)
    {
      warning("This option should not be chosen for randomization unless p.add == 0.5 ", call. = FALSE)
      min.N = max.N
      
      log.mod.switch.prob <- log(1/(max.N - min.N +1)) # probability of having that many differences
      KK<-sum(abs(varold-varnew))
      
      log.mod.switch.prob <- 0 #always the same probabilities for moves within thenighbourhood for p=0.5
      log.mod.switchback.prob <- 0
    } else if(switch.type == 2) #fixed N(x) inverse operator
    {  
      if(min.N!=max.N)
      {
        warning("min.N should be equal to max.N in swap type neighbourhoods min.N:=max.N", call. = FALSE)
        min.N = max.N
      }
      log.mod.switch.prob <- log(1/(max.N - min.N +1))
      KK<-max.N
      log.mod.switch.prob <- log.mod.switch.prob + KK*log(factorial(Nvars - KK + 1)/factorial(Nvars)) 
      log.mod.switchback.prob <- log.mod.switch.prob 
    }else if(switch.type == 3)  # random sized inverse N(x)
    {
      log.mod.switch.prob <- log(1/(max.N - min.N +1))
      KK<-sum(abs(varold-varnew))
      log.mod.switch.prob <- log.mod.switch.prob + KK*log(factorial(Nvars - KK + 1)/factorial(Nvars)) 
      log.mod.switchback.prob <- log.mod.switch.prob
    }else if(switch.type == 4)  # fixed N(x) for reverse from type 2 swaps
    {
      if(min.N!=max.N)
      {
        warning("min.N should be equal to max.N in swap type neighbourhoods min.N:=max.N", call. = FALSE)
        min.N = max.N
      }
      log.mod.switch.prob <- log(1/(max.N - min.N +1))
      KK<-max.N
      log.mod.switch.prob <- log.mod.switch.prob + KK*log(factorial(Nvars - KK + 1)/factorial(Nvars)) 
      log.mod.switchback.prob <- log.mod.switch.prob 
    }else if(switch.type >  4)
    {
      
      log.mod.switch.prob <- log(x = 1)
      log.mod.switchback.prob <-log(x = 1)
    }
    return(list(log.switch.forw.prob = log.mod.switch.prob, log.switch.back.prob = log.mod.switchback.prob))
  },
  #build a new model to draw from the old one return selection probabilities
  buildmodel= function(varcur.old,statid, shift.cpu,max.cpu,switch.type,min.N,max.N,changeble.coord = NULL)
  { 
    vect<- vector(length = max.cpu,mode = "list")
    shift<-0
    for(cpu in 1:(max.cpu))
    {
      set.seed(runif(1,1,10000), kind = NULL, normal.kind = NULL)
      if(!is.null(varcur.old))
      {
        varcur<-varcur.old
      }else
      {
        varcur<-rbinom(n = (Nvars),size = 1,0.5)
        varcur.old<-varcur
      }
      
      changeble <-FALSE
      if(!is.null(changeble.coord))
      {
        changeble <- TRUE
      }
      if(switch.type == 1) # random size random N(x) # try avoiding when doing mcmc rather than optimization
      {
        log.mod.switch.prob <- log(1/(max.N - min.N +1))
        KK<-floor(runif(n=1, min.N, max.N + 0.999999999))
        log.mod.switch.prob <- log.mod.switch.prob + KK*log(factorial(Nvars - KK + 1)/factorial(Nvars)) 
        log.mod.switchback.prob <- log.mod.switch.prob
        change.buf <- array(data = 0,dim = Nvars)
        if(changeble){
          for(ttt in 1:KK)
          {
            iid<-floor(runif(n = 1,min = 1,max = Nvars+0.999999999))
            if(changeble.coord[iid]==1)
              next
            change.buf[iid] = 1
            varcur[iid] <-rbinom(n = 1,size = 1,prob = p.add[iid])
            log.mod.switch.prob = 0#log.mod.switch.prob + log(dbinom(x = varcur[iid],size = 1,prob = p.add[iid])) # this is one of the pathes to get there in general
            log.mod.switchback.prob = 0# log.mod.switchback.prob + log(dbinom(x = 1-varcur[iid],size = 1,prob = p.add[iid])) # but this is one of the ways to get back only
          }
        }else
        {
          iid<-floor(runif(n = KK,min = 1,max = Nvars+0.999999999))
          change.buf[iid] = 1
          varcur[iid] <-rbinom(n = KK,size = 1,prob = p.add)
        }
      }else if(switch.type == 2) # fixed sized inverse N(x)
      {
        if(min.N!=max.N)
        {
          warning("min.N should be equal to max.N in swap type neighbourhoods min.N:=max.N", call. = FALSE)
          min.N <- max.N
        }
        log.mod.switch.prob <- log(1/(max.N - min.N +1))
        KK<-max.N
        log.mod.switch.prob <- log.mod.switch.prob + KK*log(factorial(Nvars - KK + 1)/factorial(Nvars)) 
        log.mod.switchback.prob <- log.mod.switch.prob
        
        change.buf <- array(data = 0,dim = Nvars)
        if(changeble){
          for(ttt in 1:KK)
          {
            iid<-floor(runif(n = 1,min = 1,max = Nvars+0.999999999))
            if(change.buf[iid]==1 || changeble.coord[iid]==1)
            {
              KK<-KK+1
              next
            }
            change.buf[iid] = 1
            varcur[iid] = 1-varcur[iid]
          }
        }else
        {
          iid<-floor(runif(n = KK,min = 1,max = Nvars+0.999999999))
          change.buf[iid] = 1
          varcur[iid] = 1-varcur[iid]
        }
        
      }else if(switch.type == 3)  # random sized inverse N(x)
      {
        log.mod.switch.prob <- log(1/(max.N - min.N +1))
        KK<-floor(runif(n=1, min.N, max.N + 0.999999999))
        log.mod.switch.prob <- log.mod.switch.prob + KK*log(factorial(Nvars - KK + 1)/factorial(Nvars)) 
        log.mod.switchback.prob <- log.mod.switch.prob
        change.buf <- array(data = 0,dim = Nvars)
        if(changeble){
          for(ttt in 1:KK)
          {
            iid<-floor(runif(n = 1,min = 1,max = Nvars+0.999999999))
            if(change.buf[iid]==1 || changeble.coord[iid]==1)
            {
              KK<-KK+1
              next
            }
            change.buf[iid] = 1
            varcur[iid] = 1-varcur[iid]
          }
        }else
        {
          iid<-floor(runif(n = KK,min = 1,max = Nvars+0.999999999))
          change.buf[iid] = 1
          varcur[iid] = 1-varcur[iid]
        }
      }else if(switch.type == 4)  # fixed N(x) for reverse from type 2 swaps
      {
        if(min.N!=max.N)
        {
          warning("min.N should be equal to max.N in swap type neighbourhoods min.N:=max.N", call. = FALSE)
          min.N <- max.N
        }
        log.mod.switch.prob <- log(1/(max.N - min.N +1))
        KK<-max.N
        log.mod.switch.prob <- log.mod.switch.prob + KK*log(factorial(Nvars - KK + 1)/factorial(Nvars)) 
        log.mod.switchback.prob <- log.mod.switch.prob
        ids <- which(changeble.coord == 0)
        change.buf <- array(data = 0,dim = Nvars)
        if(changeble){
          for(ttt in KK)
          {
            iid<-floor(runif(n = 1,min = 1,max = Nvars+0.999999999))
            if(change.buf[iid]==1 || changeble.coord[iid]==1)
            {
              KK<-KK+1
              next
            }
            change.buf[iid] = 1
            varcur[iid] = 1-varcur[iid]
          }}
        else
        {
          iid<-floor(runif(n = KK,min = 1,max = Nvars+0.999999999))
          change.buf[iid] = 1
          varcur[iid] = 1-varcur[iid]
        }
      }else if(switch.type == 5)
      {
        change.buf <- array(data = 0,dim = (Nvars))
        log.mod.switch.prob<-0
        log.mod.switchback.prob <-0
        add<-0
        if(cpu+shift>Nvars)
        {
          shift<-1 - cpu
          varcur.old<-rep(0,times = Nvars)
        }
        if(varcur.old[cpu+shift]!=1)
        {
          varcur<-varcur.old
          varcur[cpu+shift]<-1
        }else
        {
          if(cpu+shift<Nvars)
            shift<-shift+1
          while(varcur.old[cpu+shift]==1)
          {
            shift<-shift+1
            if(cpu+shift>=Nvars)
            {
              shift <- (Nvars - cpu)
              break
            }
          }
          varcur<-varcur.old
          varcur[cpu+shift]<-1
        }
        
      }else if(switch.type == 6)
      {
        change.buf <- array(data = 0,dim = (Nvars))
        log.mod.switch.prob<-0
        log.mod.switchback.prob <-0
        add<-0
        if(cpu+shift>Nvars)
        {
          shift<-1 - cpu
          varcur.old<-rep(1,times = Nvars)
        }
        
        if(varcur.old[cpu+shift]!=0)
        {
          varcur<-varcur.old
          varcur[cpu+shift]<-0
        }else
        {
          
          if(cpu+shift<Nvars)
            shift<-shift+1
          while(varcur.old[cpu+shift]==0)
          {
            shift<-shift+1
            if(cpu+shift>=Nvars)
            {
              shift <- (Nvars - cpu)
              break
            }
          }
          varcur<-varcur.old
          varcur[cpu+shift]<-0
        }
        
      }else if(switch.type == 7)
      {
        change.buf <- array(data = 0,dim = (Nvars))
        log.mod.switch.prob<-0
        log.mod.switchback.prob <-0
        vec<-dectobit(cpu + shift.cpu)
        varcur<-c(array(0,dim = (Nvars -length(vec))),vec) #issues here
      }else if(switch.type == 8)
      {
        
        log.mod.switch.prob <-0
        log.mod.switchback.prob <-0
        varcur<-varcur.old
        change.buf <- array(data = 0,dim = Nvars)
        #if(printable.opt)print("type 8 invoked")
        #if(printable.opt)print(varcur)
      }else
      {
        
        log.mod.switch.prob <- 0
        log.mod.switchback.prob <-0
        change.buf <- array(data = 1,dim = Nvars)
        varcur <-rbinom(n = Nvars,size = 1,prob = p.add)
        #if(printable.opt)print("type 8 invoked")
        #if(printable.opt)print(varcur)
      }
      
      #     for(g in 1:max(isobsbinary))
      #     {
      #       if(length(varcur[which(isobsbinary == g && varcur %in% c(1,3))])==length(which(isobsbinary == g)))
      #         varcur[which(isobsbinary == g && varcur %in% c(1,3))[1]]=0
      #       if(length(varcur[which(isobsbinary == g && varcur %in% c(2,3))])==length(which(isobsbinary == g)))
      #         varcur[which(isobsbinary == g && varcur %in% c(2,3))[1]]=0
      #     }
      
      
      
      covobs <- if(fparam[1]=="Const")fparam[which(varcur[-1] == 1)+1]else fparam[which(varcur==1)]
      
      
      #obsconst<-2*as.integer((varcur[1] %in% c(1,3)) || length(covobs) == 0 || (length(covobs)==length(which(isobsbinary != 0))) ) -1 
      
      obsconst<-ifelse(fparam[1]=="Const",2*as.integer((varcur[1])) -1,1)
      
      
      id<-bittodec(varcur)
      id<-id+1
      
      
      
      if(ifelse(exists("statistics1"),is.na(statistics1[id,4]),ifelse(exists("statistics"),is.na(statistics[id,4]),TRUE)))#||TRUE)
      {
        formula <- NULL
        capture.output({withRestarts(tryCatch(capture.output({formula <- as.formula(paste(paste(fobserved[1]), " ~ ",obsconst,ifelse(length(covobs)>0," + ",""), paste(covobs, collapse=" + "), latent.formula)) })), abort = function(){onerr<-TRUE;fm<-NULL})}) ## not considered currently in RJMCMC, is only valid for model selection
        if(is.null(formula)){
          formula <- as.formula(paste(paste(fobserved[1]), " ~ ",obsconst,ifelse(length(covobs)>0," + ",""), paste(covobs, collapse=" + ")))
          
        }
        
      }else
      {
        formula <- NULL
      }
      
      vect[[cpu]]<-list(formula = formula, varcur = varcur, statid = statid, changed = change.buf, log.mod.switch.prob = log.mod.switch.prob , log.mod.switchback.prob = log.mod.switchback.prob )
      
    }
    return(vect)
  },
  #fit the given model by means of the specified estimator
  fitmodel = function(model)
  { 
    if(!is.null(model))
    {
      
      fm<-NULL
      capture.output({withRestarts(tryCatch(capture.output({id<-bittodec(model$varcur)})), 
                                   abort = function(){onerr<-TRUE})})
      
      if(is.null(id))
        id = 0
      id<-id+1
      
      if(exists("statistics1")){
        if(is.na(statistics1[id,4]))
        { 
          
          if(printable.opt)print("Invoked from EMJMCMC environment")
          onerr<-FALSE
          #if(printable.opt)print("INLA internal error")
          statistics1[id,c(2,3)]<-100000
          statistics1[id,1]<--100000
          statistics1[id,4:14]<-0
          #prior = "normal",param = c(model$beta.mu.prior,model$beta.tau.prior))                                                                       
          
          #       capture.output({withRestarts(tryCatch(capture.output({fm<-inla(formula = model$formula,family = "binomial",Ntrials = data$total_bases,data = data,control.fixed = list(mean = list(default = model$beta.mu.prior),mean.intercept = model$beta.mu.prior, prec = list(default = model$beta.tau.prior), prec.intercept = model$beta.tau.prior) ,control.compute = list(dic = model$dic.t, waic = model$waic.t, mlik = model$mlik.t))
          #       })), abort = function(){onerr<-TRUE;fm<-NULL})}) # fit the modal, get local improvements
          #       
          capture.output({withRestarts(tryCatch(capture.output({fm<-do.call(estimator, c(estimator.args, model$formula))
          })), abort = function(){onerr<-TRUE;fm<-NULL})}) # fit the modal, get local improvements
          
          
          
          if(!is.null(fm))
          {
            statistics1[id,2]<-fm$waic[[1]]
            statistics1[id,1]<-fm$mlik[[1]]
            statistics1[id,3]<-fm$dic[[1]]
            #                                                          if(id>1)
            #                                                          {
            #                                                            inxx<-which(model$varcur==1)
            #                                                            if(length(inxx)==length(fm$summary.fixed$mean))
            #                                                              statistics1[id,15+inxx]<-fm$summary.fixed$mean 
            #                                                          }
            
            if(fm$waic[[1]]<g.results[2,1] && !is.na(fm$waic[[1]]))
            {
              g.results[2,1]<-fm$waic[[1]]
              g.results[2,2]<-as.integer(id)
            }
            if(fm$mlik[[1]]>g.results[1,1] && !is.na(fm$mlik[[1]]))
            {
              g.results[1,1]<-fm$mlik[[1]]
              g.results[1,2]<-as.integer(id)
            }
            
            if(fm$dic[[1]]<g.results[3,1]&& !is.na(fm$dic[[1]]))
            {
              g.results[3,1]<-fm$dic[[1]]
              g.results[3,2]<-as.integer(id)
            }
            
            g.results[4,2] <- g.results[4,2]+1
            if(g.results[4,2]%%recalc.margin == 0)
            {
              p.add <<- as.array(post_proceed_results(statistics1)$p.post)
            }
            
          }
          
          
        }
        if(model$statid!=-1)
          statistics1[id,model$statid+1]<-statistics1[id,model$statid+1] + 1
        g.results[4,1] <- g.results[4,1]+1
        return(list(mlik=statistics1[id,1],waic=statistics1[id,2],dic=statistics1[id,3]))
      }else  if(exists("statistics")){
        if(is.na(statistics[id,4]))
        { 
          if(printable.opt)print("Invoked from EMJMCMC SUB environment")
          onerr<-FALSE
          #if(printable.opt)print("INLA internal error")
          statistics[id,c(2,3)]<-100000
          statistics[id,1]<--100000
          statistics[id,4:14]<-0
          #prior = "normal",param = c(model$beta.mu.prior,model$beta.tau.prior))                                                                       
          
          #       capture.output({withRestarts(tryCatch(capture.output({fm<-inla(formula = model$formula,family = "binomial",Ntrials = data$total_bases,data = data,control.fixed = list(mean = list(default = model$beta.mu.prior),mean.intercept = model$beta.mu.prior, prec = list(default = model$beta.tau.prior), prec.intercept = model$beta.tau.prior) ,control.compute = list(dic = model$dic.t, waic = model$waic.t, mlik = model$mlik.t))
          #       })), abort = function(){onerr<-TRUE;fm<-NULL})}) # fit the modal, get local improvements
          #       
          capture.output({withRestarts(tryCatch(capture.output({fm<-do.call(estimator, c(estimator.args, model$formula))
          })), abort = function(){onerr<-TRUE;fm<-NULL})}) # fit the modal, get local improvements
          
          
          
          if(!is.null(fm))
          {
            statistics[id,2]<-fm$waic[[1]]
            statistics[id,1]<-fm$mlik[[1]]
            statistics[id,3]<-fm$dic[[1]]
            #                                                          if(id>1)
            #                                                          {
            #                                                            inxx<-which(model$varcur==1)
            #                                                            if(length(inxx)==length(fm$summary.fixed$mean))
            #                                                              statistics[id,14+inxx]<-fm$summary.fixed$mean 
            #                                                          }
            if(fm$waic[[1]]<g.results[2,1] && !is.na(fm$waic[[1]]))
            {
              g.results[2,1]<-fm$waic[[1]]
              g.results[2,2]<-as.integer(id)
            }
            if(fm$mlik[[1]]>g.results[1,1] && !is.na(fm$mlik[[1]]))
            {
              g.results[1,1]<-fm$mlik[[1]]
              g.results[1,2]<-as.integer(id)
            }
            
            if(fm$dic[[1]]<g.results[3,1]&& !is.na(fm$dic[[1]]))
            {
              g.results[3,1]<-fm$dic[[1]]
              g.results[3,2]<-as.integer(id)
            }
            
            g.results[4,2] <- g.results[4,2]+1
            if(g.results[4,2]%%recalc.margin == 0)
            {
              proceeeded <- post_proceed_results(statistics)
              p.add <<- as.array(proceeeded$p.post)
              #g.results[4,2] <- 
            }
          }
          
          
        }
        if(model$statid!=-1)
          statistics[id,model$statid+1]<-statistics[id,model$statid+1] + 1
        g.results[4,1] <- g.results[4,1]+1
        return(list(mlik=statistics[id,1],waic=statistics[id,2],dic=statistics[id,3]))
        
      }
      else
      {
        capture.output({withRestarts(tryCatch(capture.output({fm<-do.call(estimator, c(estimator.args, model$formula))
        })), abort = function(){onerr<-TRUE;fm<-NULL})}) # fit the modal, get local improvements
        
        if(!is.null(fm)){
          
          if(fm$waic[[1]]<g.results[2,1] && !is.na(fm$waic[[1]]))
          {
            g.results[2,1]<-fm$waic[[1]]
            g.results[2,2]<-as.integer(id)
          }
          if(fm$mlik[[1]]>g.results[1,1] && !is.na(fm$mlik[[1]]))
          {
            g.results[1,1]<-fm$mlik[[1]]
            g.results[1,2]<-as.integer(id)
          }
          
          if(fm$dic[[1]]<g.results[3,1]&& !is.na(fm$dic[[1]]))
          {
            g.results[3,1]<-fm$dic[[1]]
            g.results[3,2]<-as.integer(id)
          }
          g.results[4,1] <- g.results[4,1]+1
          g.results[4,2] <- g.results[4,2]+1
          return(list(mlik=fm$mlik[[1]],waic=fm$waic[[1]],dic=fm$dic[[1]]))
        }
        else
        {
          g.results[4,1] <- g.results[4,1]+1
          return(list(mlik=-Inf,waic=Inf,dic=Inf))
        }
      }
    }
    g.results[4,1] <- g.results[4,1]+1
    g.results[4,2] <- g.results[4,2]+1
    return(list(mlik=-Inf,waic=Inf,dic=Inf))
  },
  #lambda function for mtmcmc
  lambda = function(c,alpha,g1,g2,g.domain.pos) # simmetric choice driving function
  {
    if((c!=0))
    {
      res<-((((1/c)*(1+(g1+g2)*(g.domain.pos))/(1-(g1+g2)*(1-g.domain.pos)))^alpha))
      #!#if(printable.opt)print(res)
      return(res)
    }else
    {
      #!#if(printable.opt)print(1)
      return(1)
    }
  },
  #norm between probabilities vectors
  normprob=function(p1,p2)
  {
    nn<-abs(1-sum((p1+0.1)/(p2+0.1))/length(p1))
    if(is.na(nn))
      nn<-Inf
    return(nn)
  },
  #local mcmc procedure
  learnlocalMCMC=function(model)
  {
    M<-M.mcmc
    mlikcur<- model$mlikcur
    waiccur<-model$waiccur
    varcand<-model$varcur
    varcur<-model$varcur
    varcurb<-model$varcur
    varglob<-model$varcur
    
    modglob<-NULL
    fm<-NULL
    fmb<-NULL
    
    # estimate large jump in a reverse move
    if(model$reverse || is.infinite(mlikcur))
    {
      vectbg<-buildmodel(max.cpu = 1,varcur.old = varcur,statid = model$statid,switch.type=8, min.N = min.N,max.N = max.N) 
      if(!is.null(vectbg[[1]]$formula))
      {
        bgmod <- lapply(X = vectbg,FUN = .self$fitmodel)  
        waiccur<-bgmod[[1]]$waic
        mlikcur<-bgmod[[1]]$mlik
      }
      else
      {
        iidd<-bittodec(varcur)+1
        waiccur<-statistics1[iidd,2]
        mlikcur<-statistics1[iidd,1]
      } 
      
      
      
    }
    
    
    if(printable.opt)print(paste("Begin with ", mlikcur))
    
    mlikglob<- mlikcur
    mlikcand<- mlikcur
    waiccand<- waiccur
    waicglob<- waiccur
    waiccur<-  waiccur
    
    for(m in 1:M)
    {
      withRestarts(tryCatch({
        
        # statistics <- describe(statistics)
        vect<-buildmodel(varcur.old = varcur,statid = model$statid,max.cpu = max.cpu,switch.type=switch.type, min.N = min.N, max.N = max.N)
        cluster<-TRUE
        flag1<-0
        
        for(mod_id in 1:max.cpu)
        {
          if(is.null(vect[[mod_id]]$formula))
          {
            flag1<-flag1+1
          }
          
        }
        
        #flag1<-sum(is.null(vect[[]]$formula))
        
        if(flag1==max.cpu)
        {
          cluster<-FALSE
          if(printable.opt)print("!!!!MTMCMC models already estimated!!!!")
        }else
        {
          
          res.par <- parallelize(X = vect,FUN = .self$fitmodel)
          
        }
        p.select.y <- array(data = 0, dim = max.cpu)
        for(mod_id in 1:max.cpu)
        {
          if(cluster)
          {
            fm<-res.par[[mod_id]]
            
            if(is.null(fm)&&(is.na(res.par[[mod_id]]$waic)))
            {
              varcand<-varcurb
              if(printable.opt)print("locMTMCMC Model Fit Error!?")
              next
            }
          }
          
          varcand<-vect[[mod_id]]$varcur
          
          if(cluster)
          {
            waiccand<-res.par[[mod_id]]$waic
            mlikcand<-res.par[[mod_id]]$mlik
          }else
          {
            iidd<-bittodec(varcand)+1
            waiccand<-statistics1[iidd,2]
            mlikcand<-statistics1[iidd,1]
          }
          
          if((mlikcand>mlikglob)) #update the parameter of interest
          {
            if(printable.opt)print(paste("locMTMCMC update waic.glob = ", waiccand))
            if(printable.opt)print(paste("locMTMCMC update waic.glob.mlik = ",  mlikglob))
            mlikglob<-mlikcand
            waicglob<-waiccand
            varglob<-varcand
            if(cluster)
              modglob<-fm
          }
          
          
          g1 <- waiccur
          
          if(waiccur == Inf)
          {
            g1 = 1
          }
          
          p.select.y[mod_id]<-(mlikcand + vect[[mod_id]]$log.mod.switchback.prob+log(lambda(c = cc, alpha = aa, g1 = -g1, g2 = -waiccand,g.domain.pos =  FALSE))) # correct for different criteria later
          
          if(is.na(p.select.y[mod_id]))
            p.select.y[mod_id] <- 0
          if(is.infinite(p.select.y[mod_id]) || p.select.y[mod_id]>100000000)
          { 
            #if(printable.opt)print(paste("very large log.w.y detected ",p.select.y[mod_id]))
            p.select.y[mod_id] <- 100000000
          }
          
        }
        
        max.p.select.y <- max(p.select.y)
        p.select.y<-p.select.y-max.p.select.y
        
        #if(printable.opt)print(paste("max log.w.y is ",max.p.select.y,"normilized log.w.n.y is ", paste(p.select.y,collapse = ", ")))
        
        
        ID<-sample(x = max.cpu,size = 1,prob = exp(p.select.y))
        
        if(printable.opt)print(paste("cand ",ID," selected"))
        
        varcand<-vect[[ID]]$varcur
        
        if(cluster)
        {
          waiccand<-res.par[[ID]]$waic
          mlikcand<-res.par[[ID]]$mlik
        }else
        {
          iidd<-bittodec(varcand)+1
          waiccand<-statistics1[iidd,2]
          mlikcand<-statistics1[iidd,1]
        }
        
        #p.Q.cand<- p.select.y[ID]/sum(p.select.y)
        
        if(printable.opt)print("do reverse step")
        
        p.select.z <- array(data = 0.01, dim = max.cpu)
        
        
        if(max.cpu!=1)
        {
          vect1<-buildmodel(max.cpu = max.cpu -1,varcur.old = varcand,statid = model$statid,switch.type = switch.type, min.N = min.N, max.N = max.N)
          
          cluster<-TRUE
          
          flag1<-0
          
          for(mod_id in 1:(max.cpu-1))
          {
            if(is.null(vect1[[mod_id]]$formula))
            {
              flag1<-flag1+1
            }
            
          }
          
          if(flag1==(max.cpu-1))
          {
            cluster<-FALSE
            if(printable.opt)print("!!!!MTMCMC reverse models already estimated!!!!")
          }else
          {
            res.par.back <- parallelize(X = vect1,FUN = .self$fitmodel)
          }
          
          for(mod_id in 1:(max.cpu-1))
          {
            
            if(cluster)
            {
              if(is.null(fm)&&(is.na(res.par.back[[mod_id]]$waic)))
              {
                if(printable.opt)print("locMTMCMC Model Fit Error!?")
                next
              }
            }
            
            varcand.b<-vect1[[mod_id]]$varcur
            
            if(cluster)
            {
              waiccand.b<-res.par.back[[mod_id]]$waic
              mlikcand.b<-res.par.back[[mod_id]]$mlik
            }else
            {
              iidd<-bittodec(varcand.b)+1
              waiccand.b<-statistics1[iidd,2]
              mlikcand.b<-statistics1[iidd,1]
            }
            
            if((mlikcand.b>mlikglob))
            {
              if(printable.opt)print(paste("locMTMCMC update waic.glob = ", waiccand.b))
              if(printable.opt)print(paste("locMTMCMC update waic.glob.mlik = ", mlikcand.b))
              mlikglob<-mlikcand.b
              waicglob<-waiccand.b
              varglob<-varcand.b
              if(cluster)
                modglob<-fm
            }
            
            g1 = waiccand
            
            if(waiccand == Inf)
            {
              g1 = 1
            }
            
            p.select.z[mod_id]<-(mlikcand.b+vect1[[mod_id]]$log.mod.switchback.prob+(lambda(c = cc, alpha = aa, g1 = -g1, g2 = -waiccand.b,g.domain.pos = FALSE))) # correct for different criteria later
            
            if(is.na(p.select.z[mod_id]))
              p.select.z[mod_id]=0
            if(is.infinite(p.select.z[mod_id]) || p.select.z[mod_id] > 100000000)
            {
              #if(printable.opt)print(paste("very large log.w.y detected ",p.select.z[mod_id]))
              p.select.z[mod_id] <- 100000000
            }
          }
        }
        
        if( waiccur == Inf)
        {
          g1 = 1
        }
        p.select.z[max.cpu] <- (mlikcur+vect[[ID]]$log.mod.switch.prob+(lambda(c = cc, alpha = aa, g1 = -g1, g2 = -waiccand,g.domain.pos = FALSE)))
        
        if(is.na(p.select.z[mod_id]))
          p.select.z[mod_id]=0
        if(is.infinite(p.select.z[mod_id]) || p.select.z[mod_id] > 100000000)
        {
          #if(printable.opt)print(paste("very large log.w.y detected ",p.select.z[mod_id]))
          p.select.z[mod_id] <- 100000000
        }
        
        max.p.select.z <- max(p.select.z)
        p.select.z<-p.select.z-max.p.select.z
        
        if(printable.opt)print(paste("max log.w.z is ",max.p.select.z,"normilized log.w.n.z is ", paste(p.select.z,collapse = ", ")))
        
        if(log(runif(n = 1,min = 0,max = 1)) < (log(sum(exp(p.select.y)))-log(sum(exp(p.select.z)))) + max.p.select.y - max.p.select.z )
        {
          mlikcur<-mlikcand
          if(printable.opt)print(paste("locMTMCMC update ratcur = ", mlikcand))
          if(printable.opt)print(paste("locMTMCMC accept move with ", waiccand))
          varcur<-varcand
          waiccur<-waiccand
          
        }
        
      }),abort = function(){fm<-fmb;closeAllConnections();options(error=traceback);  onerr<-TRUE})
    }
    
    #if(printable.opt)print("FINISH LOCAL MTMCMC")
    
    #!#if(printable.opt)print(points)
    
    if(model$reverse == FALSE)
    {
      if(is.null(varcur))
      {
        #if(printable.opt)print("No moves acceoted in the procedure")
        varcur<-varcand
        waiccur<-waiccand
        varcur<-varcand
      }
      
      vect<-buildmodel(max.cpu = 1,varcur.old = varcur,statid = model$statid,switch.type = type.randomize,min.N = min.N.randomize,max.N = max.N.randomize)
      
      
      varcur<-vect[[1]]$varcur
      #if(printable.opt)print(varcur)
      
      cluster<-TRUE
      
      
      if(is.null(vect[[1]]$formula))
      {
        cluster<-FALSE
        if(printable.opt)print("!!!!MTMCMC reverse models already estimated!!!!")
      }else
      {
        mod<-fitmodel(vect[[1]])
      }
      
      if(cluster)
      {
        waiccur<-mod$waic
        mlikcur<-mod$mlik
      } else
      {
        iidd<-bittodec(varcur)+1
        waiccur<-statistics1[iidd,2]
        mlikcur<-statistics1[iidd,1]
      }
      # incorporate what happens for the backward optimization
      
      model.prob<-vect[[1]]$log.mod.switch.prob
      model.prob.fix<-vect[[1]]$log.mod.switchback.prob
      
    }else  # incorporate what happens for the reverse move
    {
      if(is.null(varcur))
      {
        #if(printable.opt)print("No moves acceoted in the reverse procedure")
        varcur<-varcand
        waiccur<-waiccand
      }
      model.probs<-calculate.move.logprobabilities(varold = varcur, varnew = model$varold,switch.type = type.randomize,min.N = min.N.randomize,max.N = max.N.randomize) 
      model.prob<-model.probs$log.switch.forw.prob
      model.prob.fix<-model.probs$log.switch.back.prob
    }
    
    if(is.null(varcur))
    {
      #if(printable.opt)print("NO VARCUR OBTAINED")
      varcur<-i
      waiccur<-Inf
      varcur<-model$varold
    }
    
    return(list(varcur = varcur, waiccur = waiccur, mlikcur = mlikcur, log.prob.cur = model.prob,log.prob.fix = model.prob.fix, varglob = varglob, waicglob = waicglob, mlikglob = mlikglob, modglob = modglob))
  },
  #local simulated annealing optimization
  learnlocalSA=function(model)
  {
    t.min<-SA.param$t.min
    t<-SA.param$t.init
    dt<-SA.param$dt
    M<-SA.param$M
    varcand<-model$varcur
    varcurb<-model$varcur
    varcur<-model$varcur
    varglob<-model$varcur
    varglob<-NULL
    modglob<-NULL
    probcur<-1
    probrevcur<-1
    first.prob <- 1
    fm<-NULL
    fmb<-NULL
    
    mlikcur<- model$mlikcur
    waiccur<-model$waiccur
    # estimate large jump in a reverse move
    # estimate large jump in a reverse move
    if((model$reverse && !model$sa2) || is.infinite(mlikcur))
    {
      vectbg<-buildmodel(max.cpu = 1,varcur.old = varcur,statid = model$statid,switch.type=8, min.N = min.N,max.N = max.N) 
      if(!is.null(vectbg[[1]]$formula))
      {
        bgmod <- lapply(X = vectbg,FUN = .self$fitmodel)  
        waiccur<-bgmod[[1]]$waic
        mlikcur<-bgmod[[1]]$mlik
      }
      else
      {
        iidd<-bittodec(varcur)+1
        waiccur<-statistics1[iidd,2]
        mlikcur<-statistics1[iidd,1]
      } 
      
    }
    
    if(printable.opt)print(paste("Begin with ", mlikcur))
    mlikglob<- mlikcur
    mlikcand<- mlikcur
    waiccand<- waiccur
    waicglob<- waiccur
    waiccur<-  waiccur
    
    while(t>t.min)
    {
      if(printable.opt)print(paste("anneal to ",t))
      t.new<-t*exp(-dt)
      if(model$reverse == TRUE && t.new <= t.min)
      {
        M<-M -1
      }
      for(m in 1:M)
      {
        withRestarts(tryCatch({
          
          mmax.cpu = max.cpu
          if(model$switch.type == 5)
          {
            mmax.cpu = Nvars - sum(varcur)
          }else if(model$switch.type == 6)
          {
            mmax.cpu = sum(varcur)
          }
          if(mmax.cpu == 0)
            mmax.cpu = 1
          
          vect<-buildmodel(max.cpu = mmax.cpu,varcur.old = varcur,statid = model$statid,switch.type=model$switch.type, min.N = min.N,max.N = max.N) 
          cluster<-TRUE
          flag1<-0
          for(mod_id in 1:mmax.cpu)
          {
            if(is.null(vect[[mod_id]]$formula))
            {
              flag1<-flag1+1
            }
            
          }
          if(flag1==mmax.cpu)
          {
            cluster<-FALSE
            if(printable.opt)print("!!!!SA Models already estimated!!!!")
          }else
          {
            res.par <- parallelize(X = vect,FUN = .self$fitmodel)      
          }
          for(mod_id in 1:mmax.cpu)
          {
            if(cluster)
            {
              fmb<-fm
              fm<-res.par[[mod_id]]
              waiccand<-Inf
              if(is.null(fm)&&(is.na(res.par[[mod_id]]$waic)))
              {
                varcand<-varcurb
                if(printable.opt)print("SA Model Fit Error!?")
                next
              }
            }
            
            varcand<-vect[[mod_id]]$varcur
            if(cluster)
            {
              waiccand<-res.par[[mod_id]]$waic
              mlikcand<-res.par[[mod_id]]$mlik
            }else
            {
              iidd<-bittodec(varcand)+1
              waiccand<-statistics1[iidd,2]
              mlikcand<-statistics1[iidd,1]
              #if(is.na(waiccand))
              #if(printable.opt)print(iidd)
            }
            
            if(objective==0)
            {
              objcand<-waiccand
              objcur<-waiccur
              objglob<-waicglob
            }else
            {
              objcand<- -mlikcand
              objcur<-  -mlikcur
              objglob<- -mlikglob
            }
            
            if(t == SA.param$t.init && mod_id == 1 && m == 2)
            {
              delta<-objcand - objcur
              first.prob <- vect[[mod_id]]$log.mod.switchback.prob + log(punif(q = exp(x = delta/t),min = 0,max = 1))
            }
            
            if(objcand<objcur)
            {
              if(printable.opt)print(paste("SA accept move with ", objcand))
              waiccur<-waiccand
              varcur<-varcand
              mlikcur<-mlikcand
              if(objcand<objglob)
              {
                waicglob<-waiccand
                varglob<-varcand
                mlikglob<-mlikcand
                if(cluster)
                  modglob<-fm
                
                if(printable.opt)print(paste("SA update global optima with", objcand))
              }
            }else
            { 
              delta<-objcand - objcur
              if(runif(n = 1,min = 0,max = 1) <= exp(x = -delta/t))
              {
                
                model.probs<-calculate.move.logprobabilities(varold = varcur, varnew = varcand,switch.type = model$switch.type,min.N = min.N,max.N = max.N) 
                probcur<-model.probs$log.switch.forw.prob
                probrevcur<-model.probs$log.switch.back.prob
                waiccur<-waiccand
                varcur<-varcand
                mlikcur<-mlikcand
                if(printable.opt)print(paste("SA accept move with ", objcand))
              }
            }
            
          }
          
        }),abort = function(){varcur<-varcurb; fm<-fmb;closeAllConnections();options(error=traceback);  onerr<-TRUE})
      }
      t<-t.new   
    }
    t<-t/exp(-dt)
    
    if(model$reverse == FALSE)
    {
      model.prob<-log(punif(q = exp(x = -delta/t),min = 0,max = 1)) +  probcur # log(P(Mk,Mk-1))
      model.prob.fix<-log(punif(q = exp(x = delta/t),min = 0,max = 1)) + probrevcur # log(P(Mk-1,Mk))
      
      if(model$sa2 == TRUE)
      {
        model.prob.fix<-model.prob.fix + first.prob # correcting for the term for local improvements of type 3.
      }
    }else  # incorporate what happens for the reverse move
    {
      
      if(is.null(varcur))
      {
        if(printable.opt)print("No moves accepted in the reverse procedure")
        varcur<-model$varcur
        objcur<-model$objcur
      }
      
      delta<-objcur-model$objold
      
      
      model.probs<-calculate.move.logprobabilities(varold = varcur, varnew = model$varold, switch.type = model$switch.type,min.N = min.N,max.N = max.N)  
      model.prob<-punif(q = exp(x = -delta/t),min = 0,max = 1,log.p = TRUE) +  model.probs$log.switch.forw.prob
      model.prob.fix<-punif(q = exp(x = delta/t),min = 0,max = 1,log.p = TRUE) + model.probs$log.switch.back.prob
      
      if(model.prob==-Inf)
      {
        model.prob<- -100000000
      }
      if(model.prob.fix==-Inf)
      {   
        model.prob.fix<- -100000000
      }
    }
    
    return(list(varcur = varcur, waiccur = waiccur, mlikcur = mlikcur, log.prob.cur = model.prob,log.prob.fix = model.prob.fix, varglob = varglob, waicglob = waicglob, mlikglob = mlikglob, modglob = modglob))
  },
  #forward selection procedure
  forward_selection=function(model)
  {
    if(printable.opt)print("begin forward selection procedure")
    varcand<-model$varcur
    varcurb<-model$varcur
    varglob<-varcand
    mlikglob<- model$mlikcur
    mlikcur<- model$mlikcur
    waiccand<- model$waiccur
    waicglob<-model$waiccur
    waiccur<-model$waiccur
    waiccurb<-model$waiccur
    varglob<-NULL
    modglob<-NULL
    
    
    fm<-NULL
    fmb<-NULL
    
    ub<-bittodec(array(1,length(varcurb)))
    layer<-length(which(varcurb == 0))
    
    mlikcur<- model$mlikcur
    waiccur<-model$waiccur
    # estimate large jump in a reverse move
    # estimate large jump in a reverse move
    if(is.infinite(mlikcur))
    {
      vectbg<-buildmodel(max.cpu = 1,varcur.old = varcurb,statid = model$statid,switch.type=8, min.N = min.N,max.N = max.N) 
      if(!is.null(vectbg[[1]]$formula))
      {
        bgmod <- lapply(X = vectbg,FUN = .self$fitmodel)  
        waiccur<-bgmod[[1]]$waic
        mlikcur<-bgmod[[1]]$mlik
      }
      else
      {
        iidd<-bittodec(varcand)+1
        waiccur<-statistics1[iidd,2]
        mlikcur<-statistics1[iidd,1]
      } 
      if(!is.na(mlikcur) &&  !is.na(waiccur) )
      {
        mlikglob<- mlikcur
        mlikcur<-  mlikcur
        waiccand<- waiccur
        waicglob<- waiccur
        waiccur<-  waiccur
        waiccurb<- waiccur
      }
    }
    
    
    while(layer>0)
    {
      withRestarts(tryCatch({
        
        if(printable.opt)print(paste("proceed with layer",layer))
        if(printable.opt)print(paste("current solution is",as.character(varcand)))
        
        vect<-buildmodel(max.cpu = layer,varcur.old = varcurb,statid = model$statid, switch.type = 5,min.N = min.N, max.N = max.N)
        
        if(printable.opt)print(paste("finish preparing models at layer",layer))
        
        cluster<-TRUE
        flag1<-0
        for(mod_id in 1:layer)
        {
          if(is.null(vect[[mod_id]]$formula))
          {
            flag1<-flag1+1
          }
          
        }
        if(flag1==layer)
        {
          cluster<-FALSE
          if(printable.opt)print("!!!!forward Models already estimated!!!!")
        }else
        {
          res.par <- parallelize(X = vect,FUN = .self$fitmodel)       
        }
        
        if(printable.opt)print(paste("end forward optimizing at layer",layer))
        for(mod_id in 1:layer)
        {
          if(cluster){
            fmb<-fm
            fm<-res.par[[mod_id]]
            waiccand<-Inf
            if(is.null(fm)&&(is.na(res.par[[mod_id]]$waic)))
            {
              varcand<-varcurb
              if(printable.opt)print("forward Model Fit Error!?")
              next
            }
          }
          
          varcand<-vect[[mod_id]]$varcur
          if(cluster)
          {
            waiccand<-res.par[[mod_id]]$waic
            mlikcand<-res.par[[mod_id]]$mlik
          }else
          {
            iidd<-bittodec(varcand)+1
            waiccand<-statistics1[iidd,2]
            mlikcand<-statistics1[iidd,1]
            #if(is.na(waiccand))
            #if(printable.opt)print(iidd)
          }
          
          if(objective==0)
          {
            objcand<-waiccand
            objcur<-waiccur
            objglob<-waicglob
          }else
          {
            objcand<- -mlikcand
            objcur<-  -mlikcur
            objglob<- -mlikglob
          }
          
          
          if(objcand<=objcur || mod_id ==1)
          {
            if(printable.opt)print(paste("forward accept with ", objcand))
            objcur<-objcand
            varcurb<-varcand
            waiccur<-waiccand
            varcur<-varcand
            mlikcur<-mlikcand
            
            if(objcur<objglob)
            {
              objglob<-objcur
              waicglob<-waiccand
              varglob<-varcand
              mlikglob<-mlikcand
              if(!is.null(fm))
                modglob<-fm
              
              if(printable.opt)print(paste("forward global optima with ", objcand))
            }
          }
          #                                                         else
          #                                                          { 
          #                                                            if(waiccand<=waiccur)
          #                                                            {
          #                                                              waiccur<-waiccand
          #                                                              varcur<-varcand
          #                                                            }
          #                                                            
          #                                                          }
          
        }
        if(objcur!=objglob)
        {
          if(model$locstop)
          { 
            break 
          }else
          {
            varcurb<-varcur
          }
        }
        
      }),abort = function(){varcur<-varcurb; fm<-fmb;closeAllConnections();options(error=traceback);  onerr<-TRUE})
      
      
      layer<-layer-1
    }
    
    model.prob<-1
    
    
    model.prob.fix<-1
    
    
    return(list(varcur = varglob, waiccur = waicglob, mlikcur = mlikglob, log.prob.cur = model.prob,log.prob.fix = model.prob.fix, varglob = varglob, waicglob = waicglob, mlikglob = mlikglob, modglob = modglob))
  },
  #backward selection procedure
  backward_selection=function(model)
  {
    #if(printable.opt)print("begin backward selection procedure")
    varcand<-model$varcur
    varcurb<-model$varcur
    varglob<-varcand
    mlikglob<- model$mlikcur
    mlikcur<- model$mlikcur
    waiccur<- model$waiccur
    waicglob<-model$waiccur
    varglob<-NULL
    modglob<-NULL
    waiccurb<- model$waiccur
    
    fm<-NULL
    fmb<-NULL
    
    ub<-bittodec(array(1,length(varcurb)))
    layer<-length(which(varcurb == 1))
    
    if(is.infinite(mlikcur))
    {
      vectbg<-buildmodel(max.cpu = 1,varcur.old = varcurb,statid = model$statid,switch.type=8, min.N = min.N,max.N = max.N) 
      if(!is.null(vectbg[[1]]$formula))
      {
        bgmod <- lapply(X = vectbg,FUN = .self$fitmodel)  
        waiccur<-bgmod[[1]]$waic
        mlikcur<-bgmod[[1]]$mlik
      }
      else
      {
        iidd<-bittodec(varcand)+1
        waiccur<-statistics1[iidd,2]
        mlikcur<-statistics1[iidd,1]
      } 
      if(!is.na(mlikcur) &&  !is.na(waiccur) )
      {
        mlikglob<- mlikcur
        mlikcur<-  mlikcur
        waiccand<- waiccur
        waicglob<- waiccur
        waiccur<-  waiccur
        waiccurb<- waiccur
      }
    }
    
    
    while(layer>0)
    {
      withRestarts(tryCatch({
        
        if(printable.opt)print(paste("backward proceed with layer",layer))
        if(printable.opt)print(paste("current backward solution is",as.character(varcand)))
        vect<-buildmodel(max.cpu = layer,varcur.old = varcurb,statid = model$statid, switch.type = 6,min.N = min.N, max.N = max.N) 
        
        if(printable.opt)print(paste("finish backward preparing models at layer",layer))
        
        cluster<-TRUE
        flag1<-0
        for(mod_id in 1:layer)
        {
          if(is.null(vect[[mod_id]]$formula))
          {
            flag1<-flag1+1
          }
          
        }
        if(flag1==layer)
        {
          cluster<-FALSE
          if(printable.opt)print("!!!!backward Models already estimated!!!!")
        }else
        {
          res.par <- parallelize(X = vect,FUN = .self$fitmodel)       
        }
        if(printable.opt)print(paste("end backward optimizing at layer",layer))
        
        for(mod_id in 1:layer)
        {
          if(cluster){
            fmb<-fm
            fm<-res.par[[mod_id]]
            waiccand<-Inf
            if(is.null(fm)&&(is.na(res.par[[mod_id]]$waic)))
            {
              varcand<-varcurb
              if(printable.opt)print("backward Model Fit Error!?")
              next
            }
          }
          
          varcand<-vect[[mod_id]]$varcur
          if(cluster)
          {
            waiccand<-res.par[[mod_id]]$waic
            mlikcand<-res.par[[mod_id]]$mlik
          }else
          {
            iidd<-bittodec(varcand)+1
            waiccand<-statistics1[iidd,2]
            mlikcand<-statistics1[iidd,1]
            #if(is.na(waiccand))
            #if(printable.opt)print(iidd)
          }
          
          if(objective==0)
          {
            objcand<-waiccand
            objcur<-waiccur
            objglob<-waicglob
          }else
          {
            objcand<- -mlikcand
            objcur<-  -mlikcur
            objglob<- -mlikglob
          }
          
          
          if(objcand<=objcur|| mod_id ==1)
          {
            if(printable.opt)print(paste("backward accept with ", objcand))
            objcur<-objcand
            varcurb<-varcand
            waiccur<-waiccand
            varcur<-varcand
            mlikcur<-mlikcand
            
            if(objcur<objglob)
            {
              objglob<-objcur
              waicglob<-waiccand
              varglob<-varcand
              mlikglob<-mlikcand
              if(!is.null(fm))
                modglob<-fm
              
              if(printable.opt)print(paste("backward global optima with ", objcand))
            }
          }
          
        }
        if(objcur!=objglob)
        {
          if(model$locstop)
          { 
            break 
          }else
          {
            varcurb<-varcur
          }
        }
        
        
      }),abort = function(){varcur<-varcurb; fm<-fmb;closeAllConnections();options(error=traceback);  onerr<-TRUE})
      
      
      layer<-layer-1
    }
    
    
    
    model.prob<-1
    
    
    model.prob.fix<-1
    
    
    return(list(varcur = varglob, waiccur = waicglob, mlikcur = mlikglob, log.prob.cur = model.prob,log.prob.fix = model.prob.fix, varglob = varglob, waicglob = waicglob, mlikglob = mlikglob, modglob = modglob))
  },
  #full selection procedure 
  full_selection=function(model)
  {
    if(printable.opt)print(paste("begin full selection procedure!","Careful, ",2^Nvars," models have to be estimated"))
    if(Nvars>30)
      if(printable.opt)print("Finishing the procedure might well take forever!")
    
    varcand<-array(0,Nvars)
    varcurb<-varcand
    varglob<-varcand
    varglob<-NULL
    modglob<-NULL
    mlikglob<- model$mlikcur
    mlikcur<- model$mlikcur
    waiccand<- model$waiccur
    waicglob<-model$waiccur
    waiccur<-model$waiccur
    waiccurb<-model$waiccur
    
    
    fm<-NULL
    fmb<-NULL
    
    ubs<-as.integer(bittodec(array(1,Nvars)) + 1)
    
    ub<-model$ub
    
    totit<-as.integer(ubs/ub) + 1
    
    if(model$totalit<totit)
    {
      totit<-model$totalit
    }
    
    if(is.infinite(mlikcur))
    {
      vectbg<-buildmodel(max.cpu = 1,varcur.old = varcurb,statid = model$statid,switch.type=8, min.N = min.N,max.N = max.N) 
      if(!is.null(vectbg[[1]]$formula))
      {
        bgmod <- lapply(X = vectbg,FUN = .self$fitmodel)  
        waiccur<-bgmod[[1]]$waic
        mlikcur<-bgmod[[1]]$mlik
      }
      else
      {
        iidd<-bittodec(varcand)+1
        waiccur<-statistics1[iidd,2]
        mlikcur<-statistics1[iidd,1]
      }
      if(!is.na(waiccur)&&!is.na(mlikcur))
      {
        mlikglob<- mlikcur
        mlikcand<- mlikcur
        waiccand<- waiccur
        waicglob<- waiccur
        waiccur<-  waiccur
        waiccurb<- waiccur
      }
    }
    
    
    for(i in 1:totit)
    {
      if(ub*i>ubs)
      {
        
        ub<-ubs - ub*(i-1)-1
        if(printable.opt)print(paste("last ",ub," iterations to complete"))
        varcurb<-varcand
      }
      withRestarts(tryCatch({
        
        
        vect<-buildmodel(max.cpu = ub,varcur.old = varcurb,statid = model$statid,switch.type = 7, shift.cpu = model$ub*(i-1),min.N = min.N,max.N = max.N) 
        if(printable.opt)print(paste("proceed with full ecumeration"))
        if(printable.opt)print(paste("current solution is",as.character(varcand)))
        cluster<-TRUE
        flag1<-0
        for(mod_id in 1:ub)
        {
          if(is.null(vect[[mod_id]]$formula))
          {
            flag1<-flag1+1
          }
          
        }
        if(flag1==ub)
        {
          cluster<-FALSE
          if(printable.opt)print("!!!!full models already estimated!!!!")
        }else
        {
          res.par <- parallelize(X = vect,FUN = .self$fitmodel)       
        }
        
        if(printable.opt)print(paste("end optimizing full ecumeration"))
        
        for(mod_id in 1:ub)
        {
          if(cluster){
            fmb<-fm
            fm<-res.par[[mod_id]]
            waiccand<-Inf
            if(is.null(fm)&&(is.na(res.par[[mod_id]]$waic)))
            {
              varcand<-varcurb
              if(printable.opt)print("full Model Fit Error!?")
              next
            }
          }
          
          varcand<-vect[[mod_id]]$varcur
          if(cluster)
          {
            waiccand<-res.par[[mod_id]]$waic
            mlikcand<-res.par[[mod_id]]$mlik
          }else
          {
            iidd<-bittodec(varcand)+1
            waiccand<-statistics1[iidd,2]
            mlikcand<-statistics1[iidd,1]
            #if(is.na(waiccand))
            #if(printable.opt)print(iidd)
          }
          
          if(objective==0)
          {
            objcand<-waiccand
            objcur<-waiccur
            objglob<-waicglob
          }else
          {
            objcand<- -mlikcand
            objcur<-  -mlikcur
            objglob<- -mlikglob
          }
          
          
          if(objcand<=objcur)
          {
            if(printable.opt)print(paste("full accept with ", objcand))
            objcur<-objcand
            varcurb<-varcand
            waiccur<-waiccand
            varcur<-varcand
            mlikcur<-mlikcand
            
            if(objcur<objglob)
            {
              objglob<-objcur
              waicglob<-waiccand
              varglob<-varcand
              mlikglob<-mlikcand
              if(!is.null(fm))
                modglob<-fm
              
              if(printable.opt)print(paste("full global optima with ", objcand))
            }
          }
          varcurb<-varcand  
        }
        
        #waiccurb<-waiccur
        
      }),abort = function(){varcur<-varcurb; fm<-fmb;closeAllConnections();options(error=traceback);  onerr<-TRUE})
    }
    
    
    #                                                     if(length(which(varcur == 1))==(Nvars-1))
    #                                                     {
    #                                                       vectbg<-buildmodel(max.cpu = 1,varcur.old = varcurb,statid = model$statid,switch.type=8, min.N = min.N,max.N = max.N) 
    #                                                       if(!is.null(vectbg[[1]]$formula))
    #                                                       {
    #                                                         bgmod <- lapply(X = vectbg,FUN = .self$fitmodel)  
    #                                                         waiccand<-bgmod[[1]]$waic
    #                                                         mlikcand<-bgmod[[1]]$mlik
    #                                                       }
    #                                                       else
    #                                                       {
    #                                                         iidd<-bittodec(varcur)+1
    #                                                         waiccand<-statistics1[iidd,2]
    #                                                         mlikcand<-statistics1[iidd,1]
    #                                                       } 
    #                                                       if(objective==0)
    #                                                       {
    #                                                         objcand<-waiccand
    #                                                         objcur<-waiccur
    #                                                         objglob<-waicglob
    #                                                       }else
    #                                                       {
    #                                                         objcand<- -mlikcand
    #                                                         objcur<-  -mlikcur
    #                                                         objglob<- -mlikglob
    #                                                       }
    #                                                       
    #                                                       
    #                                                       if(objcand<=objcur)
    #                                                       {
    #                                                         if(printable.opt)print(paste("full accept with ", objcand))
    #                                                         objcur<-objcand
    #                                                         varcurb<-varcand
    #                                                         waiccur<-waiccand
    #                                                         varcur<-varcand
    #                                                         mlikcur<-mlikcand
    #                                                         
    #                                                         if(objcur<objglob)
    #                                                         {
    #                                                           objglob<-objcur
    #                                                           waicglob<-waiccand
    #                                                           varglob<-varcand
    #                                                           mlikglob<-mlikcand
    #                                                           if(!is.null(fm))
    #                                                             modglob<-fm
    #                                                           
    #                                                           if(printable.opt)print(paste("full global optima with ", objcand))
    #                                                         }
    #                                                       }
    #                                                     }
    #                                                     
    model.prob<-1
    
    
    model.prob.fix<-1
    
    
    return(list(varcur = varglob, waiccur = waicglob, mlikcur = mlikglob, log.prob.cur = model.prob,log.prob.fix = model.prob.fix, varglob = varglob, waicglob = waicglob, mlikglob = mlikglob, modglob = modglob))
  },
  #forward backward random dance
  forw_backw_walk = function(model)
  {
    varcur<-rbinom(n = Nvars,size = 1,prob = runif(n = 1,min = 0,max = model$p1))
    mlikcur<- -Inf
    waiccur<- Inf
    for(i in 1:model$steps)
    {
      fff<-forward_selection(list(varcur=rbinom(n = Nvars,size = 1,prob = runif(n = 1,min = 0,max = model$p1)),mlikcur=-Inf,waiccur =Inf,locstop = FALSE,statid=-1))
      if(objective ==1)
      {
        if(fff$mlikglob>mlikcur)
        {
          mlikcur<-fff$mlikglob
          varcur<-fff$varglob
          waiccur<-fff$waicglob
        }
      }
      else
      {
        if(fff$waicglob<waiccur)
        {
          mlikcur<-fff$mlikglob
          varcur<-fff$varglob
          waiccur<-fff$waicglob
        }
        
      }
      set.seed(i*model$steps)
      bbb<-backward_selection(list(varcur=rbinom(n = Nvars,size = 1,prob =  runif(n = 1,min = model$p1,max = 1)),mlikcur=-Inf,waiccur =Inf,locstop = FALSE,statid=-1))
      if(objective ==1)
      {
        if(bbb$mlikglob>mlikcur)
        {
          mlikcur<-bbb$mlikglob
          varcur<-bbb$varglob
          waiccur<-bbb$waicglob
        }
      }
      else
      {
        if(bbb$waicglob<waiccur)
        {
          mlikcur<-bbb$mlikglob
          varcur<-bbb$varglob
          waiccur<-bbb$waicglob
        }
        
      }
    }
    if(model$reverse == FALSE)
    {
      
      vect<-buildmodel(max.cpu = 1,varcur.old = varcur,statid = -1,switch.type = type.randomize,min.N = min.N.randomize,max.N = max.N.randomize)
      
      varcur<-vect[[1]]$varcur
      #if(printable.opt)print(varcur)
      
      cluster<-TRUE
      
      
      
      if(is.null(vect[[1]]$formula))
      {
        cluster<-FALSE
        if(printable.opt)print("!!!!Back Forw reverse model already estimated!!!!")
      }else
      {
        
        mod <- lapply(X = vect, FUN = fitmodel)
        
      }
      
      if(cluster)
      {
        waiccur<-mod[[1]]$waic
        mlikcur<-mod[[1]]$mlik
      }else
      {
        iidd<-bittodec(varcur)+1
        waiccur<-statistics1[iidd,2]
        mlikcur<-statistics1[iidd,1]
      }
      
      # incorporate what happens for the backward optimization
      
      model.prob<-vect[[1]]$log.mod.switch.prob
      model.prob.fix<-vect[[1]]$log.mod.switchback.prob
      
    }else  # incorporate what happens for the reverse move
    {
      
      model.probs<-calculate.move.logprobabilities(switch.type = type.randomize,varold = varcur, varnew = model$varold,min.N = min.N.randomize,max.N = max.N.randomize) 
      model.prob<-model.probs$log.switch.forw.prob
      model.prob.fix<-model.probs$log.switch.back.prob
      
    }
    
    return(list(varcur = varcur, waiccur = waiccur, mlikcur = mlikcur, log.prob.cur = model.prob,log.prob.fix = model.prob.fix, varglob = varcur, waicglob = waiccur, mlikglob = mlikcur))
  },  
  #local greedy optimization
  learnlocalND=function(model)
  {
    
    # Step.nd
    varcand<-model$varcur
    varglob<-model$varcur
    varcurb<-model$varcur
    mlikcand<- model$mlikcur
    waiccand<-model$waiccur
    modglob<-NULL
    fm<-NULL
    fmb<-NULL
    opt.achieved<-FALSE
    
    
    # estimate large jump in a reverse move
    if(model$reverse || is.infinite(mlikcand))
    {
      vectbg<-buildmodel(max.cpu = 1,varcur.old = varcand,statid = model$statid,switch.type=8, min.N = min.N,max.N = max.N) 
      if(!is.null(vectbg[[1]]$formula))
      {
        bgmod <- lapply(X = vectbg,FUN = .self$fitmodel)  
        waiccand<-bgmod[[1]]$waic
        mlikcand<-bgmod[[1]]$mlik
      }
      else
      {
        iidd<-bittodec(varcand)+1
        waiccand<-statistics1[iidd,2]
        mlikcand<-statistics1[iidd,1]
      } 
      
      
    }
    
    
    if(printable.opt)print(paste("Begin with ",mlikcand))
    
    mlikglob<- mlikcand
    mlikcand<- mlikcand
    waiccand<- waiccand
    waicglob<- waiccand
    waiccur<-  waiccand
    
    buf.M.nd <- M.nd
    if(model$switch.type == 5)
    {
      buf.M.nd = Nvars - sum(varcurb)
    }else if(model$switch.type == 6)
    {
      buf.M.nd = sum(varcurb)
    }
    if(buf.M.nd == 0)
      buf.M.nd = 1
    if(M.nd<buf.M.nd)
      buf.M.nd<-M.nd
    
    for(iterat in 1:buf.M.nd)
    {
      withRestarts(tryCatch({
        # statistics <- describe(statistics)
        mmax.cpu = max.cpu
        if(model$switch.type == 5)
        {
          mmax.cpu = Nvars - sum(varcurb)
        }else if(model$switch.type == 6)
        {
          mmax.cpu = sum(varcurb)
        }
        if(mmax.cpu == 0)
          mmax.cpu = 1
        vect<-buildmodel(max.cpu = mmax.cpu,varcur.old = varcurb,statid = model$statid, switch.type = model$switch.type,min.N = min.N,max.N = max.N)
        
        cluster<-TRUE
        
        flag1<-0
        
        for(mod_id in 1:mmax.cpu)
        {
          if(is.null(vect[[mod_id]]$formula))
          {
            flag1<-flag1+1
          }
          
        }
        
        if(flag1==mmax.cpu)
        {
          cluster<-FALSE
          if(printable.opt)print("!!!!Greedy Models already estimated!!!!")
        }else
        {
          res.par <- parallelize(X = vect,FUN = .self$fitmodel)   
        }
        
        for(mod_id in 1:mmax.cpu)
        {
          varcand1<-vect[[mod_id]]$varcur
          if(cluster)
          {
            waiccand1<-res.par[[mod_id]]$waic
            mlikcand1<-res.par[[mod_id]]$mlik
          }else
          {
            iidd<-bittodec(varcand1)+1
            waiccand1<-statistics1[iidd,2]
            mlikcand1<-statistics1[iidd,1]
          }
          
          if(objective==0)
          {
            objcand<-waiccand1
            objcur<-waiccand
            objglob<-waicglob
          }
          else
          {
            objcand<- -mlikcand1
            objcur<-  -mlikcand
            objglob<- -mlikglob
          }
          
          if(objcand<objcur || mod_id ==1)
          { 
            varcand<-varcand1
            waiccand<-waiccand1
            mlikcand<-mlikcand1
            if(printable.opt)print(paste("GREEDY update local optima with ", objcand))
            if(cluster)
              fm<-res.par[[mod_id]]
          }
        }
        #if(printable.opt)print(waiccand)
        
        if(objective==0)
        {
          objcand<-waiccand1
          objcur<-waiccand
          objglob<-waicglob
        }
        else
        {
          objcand<- -mlikcand1
          objcur<-  -mlikcand
          objglob<- -mlikglob
        }
        
        if(objcur<objglob)
        {
          waicglob<-waiccand
          varglob<-varcand
          varcurb<-varcand
          mlikglob<-mlikcand
          if(cluster)
            modglob<-fm
          
          if(printable.opt)print(paste("GREEDY update global optima with ", objcur))
          
        }
        
      }),abort = function(){opt.achieved <- TRUE; fm<-fmb;closeAllConnections();options(error=traceback);  onerr<-TRUE})
      
      if(objcur!=objglob)
      {
        if(locstop.nd)
        { 
          break 
        }else
        {
          varcurb<-varcand
        }
      }
      
      
    }
    
    #!#if(printable.opt)print(points)
    
    if(model$reverse == FALSE)
    {
      
      vect<-buildmodel(max.cpu = 1,varcur.old = varcand,statid = model$statid,switch.type = type.randomize,min.N = min.N.randomize,max.N = max.N.randomize)
      
      varcur<-vect[[1]]$varcur
      #if(printable.opt)print(varcur)
      
      cluster<-TRUE
      
      
      
      if(is.null(vect[[1]]$formula))
      {
        cluster<-FALSE
        if(printable.opt)print("!!!!Greedy reverse model already estimated!!!!")
      }else
      {
        
        mod <- lapply(X = vect, FUN = fitmodel)
        
      }
      
      if(cluster)
      {
        waiccur<-mod[[1]]$waic
        mlikcur<-mod[[1]]$mlik
      }else
      {
        iidd<-bittodec(varcur)+1
        waiccur<-statistics1[iidd,2]
        mlikcur<-statistics1[iidd,1]
      }
      
      # incorporate what happens for the backward optimization
      
      model.prob<-vect[[1]]$log.mod.switch.prob
      model.prob.fix<-vect[[1]]$log.mod.switchback.prob
      
    }else  # incorporate what happens for the reverse move
    {
      
      model.probs<-calculate.move.logprobabilities(switch.type = type.randomize,varold = varcand, varnew = model$varold,min.N = min.N.randomize,max.N = max.N.randomize) 
      model.prob<-model.probs$log.switch.forw.prob
      model.prob.fix<-model.probs$log.switch.back.prob
      varcur <- varcand
      waiccur <- waiccand
      mlikcur <- mlikcand
    }
    
    return(list(varcur = varcur, waiccur = waiccur, mlikcur = mlikcur, log.prob.cur = model.prob,log.prob.fix = model.prob.fix, varglob = varglob, waicglob = waicglob, mlikglob = mlikglob, modglob = modglob))
  },
  #global emjmcmc procedure for model selection (hyper heuristic logic in terms of COP)
  modejumping_mcmc=function(glob.model)
  {
    
    stm <- proc.time()
    if(printable.opt)print("Begin model selection EMJMCMC2016 procedure")
    set.seed(runif(n = 1, min = 1, max = seed), kind = NULL, normal.kind = NULL)
    acc_moves<-1
    distrib_of_proposals <- glob.model$distrib_of_proposals
    distrib_of_neighbourhoods <- glob.model$distrib_of_neighbourhoods
    # do the search and simulations accross the modes
    g.results[4,1]<- 0
    g.results[4,2]<- 0
    forward_selection(list(varcur=rep(0,length(fparam.example)),mlikcur=-Inf,waiccur =Inf,locstop = FALSE,statid=-1))
    backward_selection(list(varcur=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),mlikcur=-Inf,waiccur =Inf,locstop = FALSE,statid=-1))
    
    if(exists("statistics1"))
    {
      p.add <<- as.array(post_proceed_results(statistics1)$p.post)
    }
    
    #set up initial parameters
    if(is.null(glob.model$varcur))
    {
      vec<-dectobit(runif(n = 1,min = 1,max = 2 ^(Nvars))) # generate an initial solution
      varcur<-c(array(0,dim = (Nvars -length(vec))),vec)
    }else if(length(glob.model$varcur[which(glob.model$varcur %in% c(0,1))])==Nvars)
    {
      varcur<-glob.model$varcur
    }
    else
    {
      if(printable.opt)print("Incorrect initial solution set be the user, a random one is generated")
      vec<-dectobit(runif(n = 1,min = 1,max = 2 ^(Nvars))) # generate an initial solution
      varcur<-c(array(0,dim = (Nvars -length(vec))),vec)
      
    }
    
    varcurb<-varcur
    varglob<-varcur
    modglob<-NULL
    
    p1 = array(data = 0,dim = Nvars)
    p2 = array(data = 1,dim = Nvars)
    j<-0
    p.post<-array(data = 1,dim = Nvars)
    waiccur<-Inf
    waicglob<-Inf
    mlikcur<- -Inf
    mlikglob<- -Inf
    ratcur<- -Inf
    fm<-NULL
    eps.emp<-normprob(p1,p2)
    max.cpu.buf<-max.cpu
    delta.time <- 0
    LocImprove<-0
    LocNeighbor<-0
    max.cpu.buf<-max.cpu.glob
 
    while((eps.emp>=glob.model$eps || j<= glob.model$maxit || j <= glob.model$burnin) && delta.time < glob.model$max.time && g.results[4,1]<= glob.model$trit && g.results[4,2]<= glob.model$trest)
    {
      p1<-p.post/acc_moves
      set.seed(runif(n = 1, min = 1, max = seed*100), kind = NULL, normal.kind = NULL)
      LocImprove <- (sample(x = 5,size = 1,prob = distrib_of_proposals) - 1)
      LocNeighbor<-(sample(x = 7,size = 1,prob = distrib_of_neighbourhoods[LocImprove+1,]))
      switch.type.glob.buf = LocNeighbor
      switch.type.buf = LocNeighbor
      if(LocNeighbor == 7)
      {
        switch.type.glob.buf = 9
        switch.type.buf = 9
        
      }
      
      #if(printable.opt)print(LocImprove)
      j<-j+1
      if(j%%glob.model$print.freq == 0)
      {
        print(paste(j," iterations completed up to now after ",delta.time," cpu minutes"," best MLIK found ",mlikglob ," current mlik found ",mlikcur,  "current acceptance ratio ",acc_moves/j))             
      }
      if(j%%100==0)
        seed = runif(n = 1,min = 0,max = 100000)
      
      #withRestarts(tryCatch({
      
      varcur<-varcurb
      if(LocImprove<=3)
      {  
        vect<-buildmodel(max.cpu = 1,varcur.old = varcurb,statid = 4 + LocImprove,min.N = min.N.glob,max.N = max.N.glob,switch.type = switch.type.glob.buf)
        max.cpu.buf = 1
      }else
      {
        vect<-buildmodel(max.cpu = max.cpu.glob,varcur.old = varcurb,statid = 4 + LocImprove,min.N = min.N,max.N = max.N,switch.type = switch.type.glob.buf)
        max.cpu.buf = max.cpu.glob
      }
      
      cluster<-TRUE
      
      flag1<-0
      
      for(mod_id in 1:max.cpu.buf)
      {
        if(is.null(vect[[mod_id]]$formula))
        {
          flag1<-flag1+1
        }
        
      }
      
      if(flag1==max.cpu.glob)
      {
        cluster<-FALSE
        if(printable.opt)print("!!!!Models already estimated!!!!")
      }else
      {
        if(max.cpu.glob > 1)
          res.par <- parallelize(X = vect,FUN = .self$fitmodel)
        else
          res.par <- lapply(X = vect,FUN = .self$fitmodel)
      }
      
      if(LocImprove>3)
      {
        
        if(printable.opt)print("!!!!Proceed with no local improvements!!!!")
        p.select.y <- array(data = 0, dim = max.cpu.glob)
        for(mod_id in 1:max.cpu.glob)
        {
          if(cluster)
          {
            fm<-res.par[[mod_id]]
            
            if(is.null(fm)&&(is.na(res.par[[mod_id]]$waic)))
            {
              varcand<-varcurb
              if(printable.opt)print("GlobMTMCMC Model Fit Error!?")
              next
            }
          }
          
          varcand<-vect[[mod_id]]$varcur
          
          if(cluster)
          {
            waiccand<-res.par[[mod_id]]$waic
            mlikcand<-res.par[[mod_id]]$mlik
          }else
          {
            iidd<-bittodec(varcand)+1
            waiccand<-statistics1[iidd,2]
            mlikcand<-statistics1[iidd,1]
          }
          
          if((mlikcand>mlikglob)) #update the parameter of interest
          {
            if(printable.opt)print(paste("GlobMTMCMC update waic.glob = ", waiccand))
            if(printable.opt)print(paste("GlobMTMCMC update waic.glob.mlik = ",  mlikglob))
            mlikglob<-mlikcand
            waicglob<-waiccand
            varglob<-varcand
            if(cluster)
              modglob<-fm
          }
          
          
          g1 <- waiccur
          
          if(waiccur == Inf)
          {
            g1 = 1
          }
          
          p.select.y[mod_id]<-(mlikcand + vect[[mod_id]]$log.mod.switchback.prob+log(lambda(c = cc, alpha = aa, g1 = -g1, g2 = -waiccand,g.domain.pos =  FALSE))) # correct for different criteria later
          
          if(is.na(p.select.y[mod_id]))
            p.select.y[mod_id] <- 0
          if(is.infinite(p.select.y[mod_id]) || p.select.y[mod_id]>100000000)
          { 
            #if(printable.opt)print(paste("very large log.w.y detected ",p.select.y[mod_id]))
            p.select.y[mod_id] <- 100000000
          }
          
        }
        
        max.p.select.y <- max(p.select.y)
        p.select.y<-p.select.y-max.p.select.y
        
        #if(printable.opt)print(paste("max log.w.y is ",max.p.select.y,"normilized log.w.n.y is ", paste(p.select.y,collapse = ", ")))
        
        
        ID<-sample(x = max.cpu.glob,size = 1,prob = exp(p.select.y))
        
        if(printable.opt)print(paste("cand ",ID," selected"))
        
        varcand<-vect[[ID]]$varcur
        
        if(cluster)
        {
          waiccand<-res.par[[ID]]$waic
          mlikcand<-res.par[[ID]]$mlik
        }else
        {
          iidd<-bittodec(varcand)+1
          waiccand<-statistics1[iidd,2]
          mlikcand<-statistics1[iidd,1]
        }
        
        #p.Q.cand<- p.select.y[ID]/sum(p.select.y)
        
        if(printable.opt)print("do reverse step")
        
        p.select.z <- array(data = 0.01, dim = max.cpu.glob)
        
        
        if(max.cpu.glob!=1)
        {
          if(switch.type.glob.buf==5)
            cstm<-6
          else if(switch.type.glob.buf==6)
          {
            cstm<-5
          }
          else
          {
            cstm <- switch.type.glob.buf
          }
          vect1<-buildmodel(max.cpu = max.cpu.glob -1,varcur.old = varcand,statid = 4 + LocImprove,switch.type = cstm, min.N = min.N, max.N = max.N)
          
          cluster<-TRUE
          
          flag1<-0
          
          for(mod_id in 1:(max.cpu.glob-1))
          {
            if(is.null(vect1[[mod_id]]$formula))
            {
              flag1<-flag1+1
            }
            
          }
          
          if(flag1==(max.cpu.glob-1))
          {
            cluster<-FALSE
            if(printable.opt)print("!!!!MTMCMC reverse models already estimated!!!!")
          }else
          {
            res.par.back <- parallelize(X = vect1,FUN = .self$fitmodel)
          }
          
          for(mod_id in 1:(max.cpu.glob-1))
          {
            
            if(cluster)
            {
              if(is.null(fm)&&(is.na(res.par.back[[mod_id]]$waic)))
              {
                if(printable.opt)print("locMTMCMC Model Fit Error!?")
                next
              }
            }
            
            varcand.b<-vect1[[mod_id]]$varcur
            
            if(cluster)
            {
              waiccand.b<-res.par.back[[mod_id]]$waic
              mlikcand.b<-res.par.back[[mod_id]]$mlik
            }else
            {
              iidd<-bittodec(varcand.b)+1
              waiccand.b<-statistics1[iidd,2]
              mlikcand.b<-statistics1[iidd,1]
            }
            
            if((mlikcand.b>mlikglob))
            {
              if(printable.opt)print(paste("GlobMTMCMC update waic.glob = ", waiccand.b))
              if(printable.opt)print(paste("GlobMTMCMC update waic.glob.mlik = ", mlikcand.b))
              mlikglob<-mlikcand.b
              waicglob<-waiccand.b
              varglob<-varcand.b
              if(cluster)
                modglob<-fm
            }
            
            g1 = waiccand
            
            if(waiccand == Inf)
            {
              g1 = 1
            }
            
            p.select.z[mod_id]<-(mlikcand.b+vect1[[mod_id]]$log.mod.switchback.prob+(lambda(c = cc, alpha = aa, g1 = -g1, g2 = -waiccand.b,g.domain.pos = FALSE))) # correct for different criteria later
            
            if(is.na(p.select.z[mod_id]))
              p.select.z[mod_id]=0
            if(is.infinite(p.select.z[mod_id]) || p.select.z[mod_id] > 100000000)
            {
              #if(printable.opt)print(paste("very large log.w.y detected ",p.select.z[mod_id]))
              p.select.z[mod_id] <- 100000000
            }
          }
        }
        
        if( waiccur == Inf)
        {
          g1 = 1
        }
        p.select.z[max.cpu.glob] <- (mlikcur+vect[[ID]]$log.mod.switch.prob+(lambda(c = cc, alpha = aa, g1 = -g1, g2 = -waiccand,g.domain.pos = FALSE)))
        
        if(is.na(p.select.z[mod_id]))
          p.select.z[mod_id]=0
        if(is.infinite(p.select.z[mod_id]) || p.select.z[mod_id] > 100000000)
        {
          #if(printable.opt)print(paste("very large log.w.y detected ",p.select.z[mod_id]))
          p.select.z[mod_id] <- 100000000
        }
        
        max.p.select.z <- max(p.select.z)
        p.select.z<-p.select.z-max.p.select.z
        
        if(printable.opt)print(paste("max log.w.z is ",max.p.select.z,"normilized log.w.n.z is ", paste(p.select.z,collapse = ", ")))
        
        if(log(runif(n = 1,min = 0,max = 1)) < (log(sum(exp(p.select.y)))-log(sum(exp(p.select.z)))) + max.p.select.y - max.p.select.z )
        {
          mlikcur<-mlikcand
          ratcur<-mlikcand
          if(printable.opt)print(paste("global MTMCMC update ratcur = ", mlikcand))
          if(printable.opt)print(paste("global MTMCMC accept move with ", waiccand))
          varcurb<-varcand
          waiccur<-waiccand
          
          id<-bittodec(varcurb)+1
          
          acc_moves<-acc_moves+1
          if(j<glob.model$burnin)
          {
            distrib_of_proposals[5]<-distrib_of_proposals[5]+1
          }else
          {
            if(exists("statistics1"))
            {
              statistics1[id,14]<-statistics1[id,14] + 1
              statistics1[id,4]<-statistics1[id,4] + 1
            }
            p.post<- (p.post + varcand) 
          }
        }else if(j<glob.model$burnin && distrib_of_proposals[5]>0)
        {
          distrib_of_proposals[5]<-distrib_of_proposals[5] - 1
        }
        
      }else
      {
        
        if(cluster)
        {
          fm<-res.par[[mod_id]]
          
          if(is.null(fm)&&(is.na(res.par[[mod_id]]$waic)))
          {
            varcand<-varcurb
            if(printable.opt)print("EMJMCMC Model Fit Error!?")
            next
          }
        }
        
        varcand<-vect[[mod_id]]$varcur
        
        if(cluster)
        {
          waiccand<-res.par[[mod_id]]$waic
          mlikcand<-res.par[[mod_id]]$mlik
        }else
        {
          iidd<-bittodec(varcand)+1
          waiccand<-statistics1[iidd,2]
          mlikcand<-statistics1[iidd,1]
        }
        
        varcur<-varcand
      }
      
      # try local improvements
      if(LocImprove<=3)
      {
        # sa improvements
        if(LocImprove == 0 || LocImprove == 3)
        {
          if(printable.opt)print("Try SA imptovements")
          buf.change <- array(data = 1,dim = Nvars)
          buf.change[which(varcur - varcurb!=0)]=0
          #buf.opt[floor(runif(n = n.size,min = 1,max = Nvars+0.999999))] = 1
          
          if(objective ==0)
          {
            objold<-waiccur
            objcur<-waiccand
            
          }else
          {
            objold<- -mlikcur
            objcur<- -mlikcand
          }
          model = list(statid = 4 + LocImprove, switch.type = switch.type.buf, change = buf.change,mlikcur = mlikcur,
                       varcur = varcur,varold = varcurb, objcur = objcur,objold = objold, sa2 = ifelse(LocImprove == 3,TRUE,FALSE) ,reverse=FALSE)
          
          SA.forw<-learnlocalSA(model)          
          ratcand<-SA.forw$mlikcur
          
          if(LocImprove == 0)
          {
            
            ids = which(buf.change == 1)
            model$varcur <- SA.forw$varcur
            model$waiccur<-SA.forw$waiccur
            model$varcur[ids] = 1 - model$varcur[ids]
            
            model$mlikcur<- -Inf
            model$waiccur<- Inf
            #estimate the jump!!!
            if(objective ==0)
            {
              model$objold<-waiccur
              
            }else
            {
              model$objold<- -ratcand
            }
            
            model$reverse = TRUE
            if(switch.type.buf==5)
            {
              model$switch.type <- 6
            }
            else if(switch.type.buf==6)
            {
              model$switch.type <- 5
            }
            SA.back<-learnlocalSA(model)      
            
          }
          
          if(LocImprove == 0)
          {
            
            if(log(runif(n = 1,min = 0,max = 1))<=(ratcand - ratcur - SA.forw$log.prob.cur + SA.forw$log.prob.fix + SA.back$log.prob.cur - SA.back$log.prob.fix))
            {
              ratcur<-ratcand
              mlikcur<-ratcand
              if(printable.opt)print(paste("update ratcur through SA = ", ratcur))
              varcurb<-SA.forw$varcur
              acc_moves<-acc_moves+1
              id<-bittodec(varcurb)+1
              
              if(j<glob.model$burnin)
              {
                distrib_of_proposals[1]<-distrib_of_proposals[1]+1
              }else
              {
                if(exists("statistics1")){
                  statistics1[id,10]<-statistics1[id,10] + 1
                  statistics1[id,4]<-statistics1[id,4] + 1
                }
                p.post<- (p.post + varcand)
              }
              
            }else if(j<glob.model$burnin && distrib_of_proposals[1]>0)
            {
              distrib_of_proposals[1]<-distrib_of_proposals[1] - 1
            }
            if((SA.back$mlikglob<mlikglob))
            {
              if(printable.opt)print(paste("update waic.glob sa.back = ", SA.back$waicglob))
              if(printable.opt)print(paste("update waic.glob.mlik sa.back= ", SA.back$mlikglob))
              waicglob<-SA.back$waicglob
              varglob<-SA.back$varglob
              modglob<-SA.back$modglob
            }
            
          }else
          {
            if(log(runif(n = 1,min = 0,max = 1))<=(ratcand - ratcur - SA.forw$log.prob.cur + SA.forw$log.prob.fix + vect[[mod_id]]$log.mod.switchback.prob - vect[[mod_id]]$log.mod.switch.prob))
            {
              ratcur<-ratcand
              mlikcur<-ratcand
              if(printable.opt)print(paste("update ratcur through SA = ", ratcur))
              varcurb<-SA.forw$varcur
              
              id<-bittodec(varcurb)+1
              
              acc_moves<-acc_moves+1
              if(j<glob.model$burnin)
              {
                distrib_of_proposals[4]<-distrib_of_proposals[4]+1
              }else
              {
                if(exists("statistics1")){
                  statistics1[id,13]<-statistics1[id,13] + 1
                  statistics1[id,4]<-statistics1[id,4] + 1
                }
                p.post<- (p.post + varcand)
                
              }
              
            }else if(j<glob.model$burnin && distrib_of_proposals[4]>0)
            {
              distrib_of_proposals[4]<-distrib_of_proposals[4] - 1
            }
            
          }
          if((SA.forw$mlikglob<mlikglob))
          {
            if(printable.opt)print(paste("update waic.glob sa.forw = ", SA.forw$waicglob))
            if(printable.opt)print(paste("update waic.glob.mlik sa.forw = ", SA.forw$mlikglob))
            waicglob<-SA.forw$waicglob
            varglob<-SA.forw$varglob
            modglob<-SA.forw$modglob
          }
          
        }else  if(LocImprove == 1)
        {
          if(printable.opt)print("Try MTMCMC imptovements")
          buf.change <- array(data = 1,dim = Nvars)
          buf.change[which(varcur - varcurb!=0)]=0
          #buf.opt[floor(runif(n = n.size,min = 1,max = Nvars+0.999999))] = 1
          
          
          
          model = list(statid = 4 + LocImprove,reverse = FALSE, change = buf.change,mlikcur = mlikcur,waiccur = waiccur,
                       varcur = varcur,varold = varcurb)
          
          MTMCMC.forw<-learnlocalMCMC(model)          
          ids = which(buf.change == 1)
          model$varcur <- MTMCMC.forw$varcur
          
          
          model$varcur[ids] = 1 - model$varcur[ids]
          model$mlikcur<- -Inf
          model$waiccur<- Inf
          model$reverse = TRUE
          
          #if(printable.opt)print("learn reverse local MTMCMC")
          MTMCMC.back<-learnlocalMCMC(model)      
          
          #if(printable.opt)print("finish reverse local MTMCMC")
          
          MTMCMC.p.forw<-MTMCMC.forw$log.prob.cur
          MTMCMC.p.back<-MTMCMC.back$log.prob.fix
          ratcand<-MTMCMC.forw$mlikcur
          
          
          if(log(runif(n = 1,min = 0,max = 1))<=(ratcand - ratcur - MTMCMC.forw$log.prob.cur + MTMCMC.forw$log.prob.fix + MTMCMC.back$log.prob.cur - MTMCMC.back$log.prob.fix))
          {
            ratcur<-ratcand
            mlikcur<-ratcand
            #if(printable.opt)print(paste("update ratcur through MTMCMC = ", ratcur))
            varcurb<-MTMCMC.forw$varcur
            acc_moves<-acc_moves+1
            id<-bittodec(varcurb)+1
            
            if(j<glob.model$burnin)
            {
              distrib_of_proposals[2]<-distrib_of_proposals[2]+1
            }else
            {
              if(exists("statistics1")){
                statistics1[id,11]<-statistics1[id,11] + 1
                statistics1[id,4]<-statistics1[id,4] + 1
              }
              p.post<- (p.post + varcand)
              
            }
          }else if(j<glob.model$burnin && distrib_of_proposals[2]>1)
          {
            distrib_of_proposals[2]<-distrib_of_proposals[2] - 1
          }
          
          if((MTMCMC.forw$mlikglob<mlikglob))
          {
            #if(printable.opt)print(paste("update waic.glob MTMCMC.forw = ", MTMCMC.forw$waicglob))
            #if(printable.opt)print(paste("update waic.glob.mlik MTMCMC.forw = ", MTMCMC.forw$mlikglob))
            waicglob<-MTMCMC.forw$waicglob
            varglob<-MTMCMC.forw$varglob
            modglob<-MTMCMC.forw$modglob
          }
          
          if((MTMCMC.back$mlikglob<mlikglob))
          {
            #if(printable.opt)print(paste("update waic.glob MTMCMC.back = ", MTMCMC.back$waicglob))
            #if(printable.opt)print(paste("update waic.glob.mlik MTMCMC.back= ", MTMCMC.back$mlikglob))
            waicglob<-MTMCMC.back$waicglob
            varglob<-MTMCMC.back$varglob
            modglob<-MTMCMC.back$modglob
          }
          
        }else if(LocImprove == 2)
        {
          if(printable.opt)print("Try greedy heuristic imptovements")
          buf.change <- array(data = 1,dim = Nvars)
          buf.change[which(varcur - varcurb!=0)]=0
          
          #buf.opt <- buf.change
          #buf.opt[floor(runif(n = n.size,min = 1,max = Nvars+0.999999))] = 1
          if(objective ==0)
          {
            objold<-waiccur
            objcur<-waiccand
            
          }else
          {
            objold<- -mlikcur
            objcur<- -mlikcand
          }
          
          model = list(statid = 4 + LocImprove, change = buf.change,
                       varcur = varcur,varold = varcurb, switch.type = switch.type.buf,mlikcur = mlikcur,objcur = objcur, objold = objold, reverse=FALSE)
          GREEDY.forw<-learnlocalND(model)          
          ids = which(buf.change == 1)
          model$varcur <- GREEDY.forw$varcur
          
          model$varcur[ids] = 1 - model$varcur[ids]
          model$mlikcur<- -Inf
          model$waiccur<- Inf
          model$reverse = TRUE
          ratcand<-GREEDY.forw$mlikcur
          
          if(objective ==0)
          {
            model$objold<-waiccur
            
          }else
          {
            model$objold<- -ratcand
          }
          if(switch.type.buf==5)
          {
            model$switch.type <- 6
          }
          else if(switch.type.buf==6)
          {
            model$switch.type <- 5
          }
          
          GREEDY.back<-learnlocalND(model)      
          
          GREEDY.p.forw<-GREEDY.forw$log.prob.cur
          GREEDY.p.back<-GREEDY.back$log.prob.fix
          
          
          
          if(log(runif(n = 1,min = 0,max = 1))<=(ratcand - ratcur - GREEDY.forw$log.prob.cur + GREEDY.forw$log.prob.fix + GREEDY.back$log.prob.cur - GREEDY.back$log.prob.fix))
          {
            ratcur<-ratcand
            mlikcur<-ratcand
            if(printable.opt)print(paste("update ratcur through ND = ", ratcur))
            varcurb<-GREEDY.forw$varcur
            acc_moves<-acc_moves+1
            id<-bittodec(varcurb)+1
            
            if(j<glob.model$burnin)
            {
              distrib_of_proposals[3]<-distrib_of_proposals[3]+1
            }else
            {
              if(exists("statistics1")){
                statistics1[id,12]<-statistics1[id,12] + 1
                statistics1[id,4]<-statistics1[id,4] + 1
              }
              p.post<- (p.post + varcand)
              
            }
          }else if(j<glob.model$burnin && distrib_of_proposals[3]>1)
          {
            distrib_of_proposals[3]<-distrib_of_proposals[3] - 1
          }
          
          if((GREEDY.forw$mlikglob<mlikglob))
          {
            if(printable.opt)print(paste("update waic.glob ND.forw = ", GREEDY.forw$waicglob))
            if(printable.opt)print(paste("update waic.glob.mlik ND.forw = ", GREEDY.forw$mlikglob))
            waicglob<-GREEDY.forw$waicglob
            varglob<-GREEDY.forw$varglob
            modglob<-GREEDY.forw$modglob
          }
          
          if((GREEDY.back$mlikglob<mlikglob))
          {
            if(printable.opt)print(paste("update waic.glob ND.back = ", GREEDY.back$waicglob))
            if(printable.opt)print(paste("update waic.glob.mlik ND.back= ", GREEDY.back$mlikglob))
            waicglob<-GREEDY.back$waicglob
            varglob<-GREEDY.back$varglob
            modglob<-GREEDY.back$modglob
          }
          
        }
        
      }
      
      #}),abort = function(){if(printable.opt)print("error");varcur<-varcurb;closeAllConnections();options(error=traceback);  onerr<-TRUE})
      
      
      p2<-p.post/acc_moves
      if(j>glob.model$burnin && recalc.margin == 2^Nvars)
      {
        p.add <<- p2
      }
      eps.emp<-normprob(p1,p2)
      etm <- proc.time()
      delta.time <- (etm[3] - stm[3])/60.0
    }
    
    if(is.null(modglob))
    {
      vect<-buildmodel(max.cpu = 1,varcur.old = varglob,statid = 3, switch.type = 8)
      res.par <- lapply(X = vect, FUN = fitmodel)
      modglob<-res.par[[1]]$fm
    }
    #!#print the results
    
    if(printable.opt)print(paste(j," iterations completed in total"))
    #if(printable.opt)print(paste("WAIC.glob = ", waicglob))
    etm <- proc.time()
    tt<-(etm[3]-stm[3])/60.0
    
    acc_ratio <-  acc_moves/j
    
    if(printable.opt)print(paste(j," moves proposed in total, ", acc_moves," of them accepted, acceptance ratio is ",acc_ratio))
    
    if(printable.opt)print(paste("posterior distribution ofproposals is",  distrib_of_proposals))
    
    id<-bittodec(varglob)+1
    
    if(exists("statistics1"))
    {
      bayes.res<-post_proceed_results(statistics1)
      m.post<-statistics1[,4]/j
      
    }else
    {
      bayes.res<-NULL
      m.post<-NULL
    }
    
    return(list(model=modglob, vars=varglob, waic = waicglob, mlik = mlikglob,  time=(tt), freqs = distrib_of_proposals, acc_ratio =  acc_ratio, p.post = p.post/acc_moves, m.post =  m.post, p.post.freq = p.post, eps=eps.emp, bayes.results =  bayes.res))
  },
  #save big.data results, if the latter are available
  save_results_csv = function(statistics1, filename)
  {
    write.big.matrix(x = statistics1,filename = filename)
  },
  #vizualize the results graphically if the latter are available
  visualize_results = function(statistics1, template, mds_size, crit,draw_dist = FALSE)
  {
    ids = NULL
    
    if(crit$mlik){
      ids<-c(ids,which(statistics1[,1]==-100000))
      ids<-c(ids,which(statistics1[,1]==-Inf))}
    #ids<-c(ids,which(statistics1[,1]>0))
    if(crit$waic){
      ids<-c(ids,which(statistics1[,2]==Inf))
      ids<-c(ids,which(statistics1[,3]==Inf))}
    if(crit$dic){
      ids<-c(ids,which(is.na(statistics1[,3])))}
    
    if(length(ids)!=0)
    {
      if(crit$mlik)
        mlik.lim<-c(min(statistics1[-ids,1],na.rm = TRUE),max(statistics1[-ids,1],na.rm = TRUE))
      if(crit$waic)
        waic.lim<-c(min(statistics1[-ids,2],na.rm = TRUE),max(statistics1[-ids,2],na.rm = TRUE))
      if(crit$dic)
        dic.lim<-c(min(statistics1[-ids,3],na.rm = TRUE),max(statistics1[-ids,3],na.rm = TRUE))
    }else
    {
      if(crit$mlik)
        mlik.lim<-c(min(statistics1[,1],na.rm = TRUE),max(statistics1[,1],na.rm = TRUE))
      if(crit$waic)
        waic.lim<-c(min(statistics1[,2],na.rm = TRUE),max(statistics1[,2],na.rm = TRUE))
      if(crit$dic)
        dic.lim<-c(min(statistics1[,3],na.rm = TRUE),max(statistics1[,3],na.rm = TRUE))
    }
    
    norm<-0.5*sqrt(sum(statistics1[,4],na.rm = TRUE))
    
    if(crit$mlik)
    {
      
      if(printable.opt)print(paste("drawing ",workdir,template,"_legend.jpg",sep = ""))
      capture.output({withRestarts(tryCatch(capture.output({
        jpeg(file=paste(workdir,template,"_legend.jpg",sep = ""))
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
        
        
        if(printable.opt)print(paste("drawing ",workdir,template,"_mlik.jpg",sep = ""))
        jpeg(file=paste(workdir,template,"_mlik.jpg",sep = ""))
        plot(ylim = mlik.lim, xlab = "model_id", ylab="MLIK" , statistics1[,1],pch=19, col = 7,cex= 1*((statistics1[,9]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
        points(statistics1[,1],pch=8,  col = ifelse(statistics1[,4]>0,5,0),cex= ifelse(statistics1[,4]>0,statistics1[,4]/norm+1,0))
        points(statistics1[,1],pch=2,  col = ifelse(statistics1[,10]>0,2,0),cex= ifelse(statistics1[,10]>0,statistics1[,10]/norm+1,0))
        points(statistics1[,1],pch=3,  col = ifelse(statistics1[,11]>0,3,0),cex= ifelse(statistics1[,11]>0,statistics1[,11]/norm+1,0))
        points(statistics1[,1],pch=4,  col = ifelse(statistics1[,12]>0,4,0),cex= ifelse(statistics1[,12]>0,statistics1[,12]/norm+1,0))
        points(statistics1[,1],pch=6,  col = ifelse(statistics1[,13]>0,6,0),cex= ifelse(statistics1[,13]>0,statistics1[,13]/norm+1,0))
        points(statistics1[,1],pch=1,  col = ifelse(statistics1[,14]>0,1,0),cex= ifelse(statistics1[,14]>0,statistics1[,14]/norm+1,0))
        dev.off()
        
      })), abort = function(){onerr<-TRUE})})
    }
    
    if(crit$waic)
    {
      
      if(printable.opt)print(paste("drawing ",workdir,template,"_waic.jpg",sep = ""))
      capture.output({withRestarts(tryCatch(capture.output({
        jpeg(file=paste(workdir,template,"_waic.jpg",sep = ""))
        plot(ylim = waic.lim, xlab = "model_id", ylab="WAIC" ,statistics1[,2],pch=19, col = 7,cex= 1*((statistics1[,9]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
        points(statistics1[,2],pch=8,  col = ifelse(statistics1[,4]>0,5,0),cex= ifelse(statistics1[,4]>0,statistics1[,4]/norm+1,0))
        points(statistics1[,2],pch=2,  col = ifelse(statistics1[,10]>0,2,0),cex= ifelse(statistics1[,10]>0,statistics1[,10]/norm+1,0))
        points(statistics1[,2],pch=3,  col = ifelse(statistics1[,11]>0,3,0),cex= ifelse(statistics1[,11]>0,statistics1[,11]/norm+1,0))
        points(statistics1[,2],pch=4,  col = ifelse(statistics1[,12]>0,4,0),cex= ifelse(statistics1[,12]>0,statistics1[,12]/norm+1,0))
        points(statistics1[,2],pch=6,  col = ifelse(statistics1[,13]>0,6,0),cex= ifelse(statistics1[,13]>0,statistics1[,13]/norm+1,0))
        points(statistics1[,2],pch=1,  col = ifelse(statistics1[,14]>0,1,0),cex= ifelse(statistics1[,14]>0,statistics1[,14]/norm+1,0))
        dev.off()
      })), abort = function(){onerr<-TRUE})})
    }
    
    if(crit$dic)
    {
      
      if(printable.opt)print(paste("drawing ",workdir,template,"_dic.jpg",sep = ""))
      capture.output({withRestarts(tryCatch(capture.output({
        jpeg(file=paste(workdir,template,"_dic.jpg",sep = ""))
        plot(ylim = dic.lim, xlab = "model_id", ylab="DIC" ,statistics1[,3],pch=19, col = 7,cex= 1*((statistics1[,9]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
        points(statistics1[,3],pch=8,  col = ifelse(statistics1[,4]>0,5,0),cex= ifelse(statistics1[,4]>0,statistics1[,4]/norm+1,0))
        points(statistics1[,3],pch=2,  col = ifelse(statistics1[,10]>0,2,0),cex= ifelse(statistics1[,10]>0,statistics1[,10]/norm+1,0))
        points(statistics1[,3],pch=3,  col = ifelse(statistics1[,11]>0,3,0),cex= ifelse(statistics1[,11]>0,statistics1[,11]/norm+1,0))
        points(statistics1[,3],pch=4,  col = ifelse(statistics1[,12]>0,4,0),cex= ifelse(statistics1[,12]>0,statistics1[,12]/norm+1,0))
        points(statistics1[,3],pch=6,  col = ifelse(statistics1[,13]>0,6,0),cex= ifelse(statistics1[,13]>0,statistics1[,13]/norm+1,0))
        points(statistics1[,3],pch=1,  col = ifelse(statistics1[,14]>0,1,0),cex= ifelse(statistics1[,14]>0,statistics1[,14]/norm+1,0))
        dev.off()
      })), abort = function(){onerr<-TRUE})})
    }
    
    if(crit$waic && crit$mlik)
    {
      
      if(printable.opt)print(paste("drawing ",workdir,template,"_mlik-waic.jpg",sep = ""))
      capture.output({withRestarts(tryCatch(capture.output({
        jpeg(file=paste(workdir,template,"_mlik-waic.jpg",sep = ""))
        plot(ylim = mlik.lim, xlim = waic.lim,xlab = "WAIC", ylab="MLIK" ,statistics1[,2],statistics1[,1], pch=19, col = 7,cex= 1*((statistics1[,9]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
        points(statistics1[,2],statistics1[,1],pch=8,  col = ifelse(statistics1[,4]>0,5,0),cex= ifelse(statistics1[,4]>0,statistics1[,4]/norm+1,0))
        points(statistics1[,2],statistics1[,1],pch=2,  col = ifelse(statistics1[,10]>0,2,0),cex= ifelse(statistics1[,10]>0,statistics1[,10]/norm+1,0))
        points(statistics1[,2],statistics1[,1],pch=3,  col = ifelse(statistics1[,11]>0,3,0),cex= ifelse(statistics1[,11]>0,statistics1[,11]/norm+1,0))
        points(statistics1[,2],statistics1[,1],pch=4,  col = ifelse(statistics1[,12]>0,4,0),cex= ifelse(statistics1[,12]>0,statistics1[,12]/norm+1,0))
        points(statistics1[,2],statistics1[,1],pch=6,  col = ifelse(statistics1[,13]>0,6,0),cex= ifelse(statistics1[,13]>0,statistics1[,13]/norm+1,0))
        points(statistics1[,2],statistics1[,1],pch=1,  col = ifelse(statistics1[,14]>0,1,0),cex= ifelse(statistics1[,14]>0,statistics1[,14]/norm+1,0))
        dev.off()
      })), abort = function(){onerr<-TRUE})})
    }
    
    if(crit$dic && crit$mlik)
    {
      
      if(printable.opt)print(paste("drawing ",workdir,template,"_mlik-dic.jpg",sep = ""))
      capture.output({withRestarts(tryCatch(capture.output({
        jpeg(file=paste(workdir,template,"_mlik-dic.jpg",sep = ""))
        plot(ylim = mlik.lim, xlim = dic.lim,xlab = "DIC", ylab="MLIK" ,statistics1[,3],statistics1[,1], pch=19, col = 7,cex= 1*((statistics1[,9]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
        points(statistics1[,3],statistics1[,1],pch=8,  col = ifelse(statistics1[,4]>0,5,0),cex= ifelse(statistics1[,4]>0,statistics1[,4]/norm+1,0))
        points(statistics1[,3],statistics1[,1],pch=2,  col = ifelse(statistics1[,10]>0,2,0),cex= ifelse(statistics1[,10]>0,statistics1[,10]/norm+1,0))
        points(statistics1[,3],statistics1[,1],pch=3,  col = ifelse(statistics1[,11]>0,3,0),cex= ifelse(statistics1[,11]>0,statistics1[,11]/norm+1,0))
        points(statistics1[,3],statistics1[,1],pch=4,  col = ifelse(statistics1[,12]>0,4,0),cex= ifelse(statistics1[,12]>0,statistics1[,12]/norm+1,0))
        points(statistics1[,3],statistics1[,1],pch=6,  col = ifelse(statistics1[,13]>0,6,0),cex= ifelse(statistics1[,13]>0,statistics1[,13]/norm+1,0))
        points(statistics1[,3],statistics1[,1],pch=1,  col = ifelse(statistics1[,14]>0,1,0),cex= ifelse(statistics1[,14]>0,statistics1[,14]/norm+1,0))
        dev.off()
      })), abort = function(){onerr<-TRUE})})
    }
    
    if(crit$dic && crit$waic)
    {
      
      if(printable.opt)print(paste("drawing ",workdir,template,"_waic-dic.jpg",sep = ""))
      capture.output({withRestarts(tryCatch(capture.output({
        jpeg(file=paste(workdir,template,"_waic-dic.jpg",sep = ""))
        plot(ylim = waic.lim, xlim = dic.lim,xlab = "WAIC", ylab="DIC" ,statistics1[,3],statistics1[,2], pch=19, col = 7,cex= 1*((statistics1[,9]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
        points(statistics1[,3],statistics1[,2],pch=8,  col = ifelse(statistics1[,4]>0,5,0),cex= ifelse(statistics1[,4]>0,statistics1[,4]/norm+1,0))
        points(statistics1[,3],statistics1[,2],pch=2,  col = ifelse(statistics1[,10]>0,2,0),cex= ifelse(statistics1[,10]>0,statistics1[,10]/norm+1,0))
        points(statistics1[,3],statistics1[,2],pch=3,  col = ifelse(statistics1[,11]>0,3,0),cex= ifelse(statistics1[,11]>0,statistics1[,11]/norm+1,0))
        points(statistics1[,3],statistics1[,2],pch=4,  col = ifelse(statistics1[,12]>0,4,0),cex= ifelse(statistics1[,12]>0,statistics1[,12]/norm+1,0))
        points(statistics1[,3],statistics1[,2],pch=6,  col = ifelse(statistics1[,13]>0,6,0),cex= ifelse(statistics1[,13]>0,statistics1[,13]/norm+1,0))
        points(statistics1[,3],statistics1[,2],pch=1,  col = ifelse(statistics1[,14]>0,1,0),cex= ifelse(statistics1[,14]>0,statistics1[,14]/norm+1,0))
        dev.off()
      })), abort = function(){onerr<-TRUE})})
    }
    
    
    norm1<-(sum(statistics1[,4],na.rm = TRUE))
    
    xyz<-which(!is.na(statistics1[,1]))
    xyz<-intersect(xyz,which(statistics1[,1]!=-10000))
    moddee<-which(statistics1[,1]==max(statistics1[,1],na.rm = TRUE))[1]
    zyx<-array(data = NA,dim = 2 ^(Nvars)+1)
    nconsum<-sum(exp(-statistics1[moddee,1]+statistics1[xyz,1]),na.rm = TRUE)
    if( nconsum > 0)
    {
      zyx[xyz]<-exp(statistics1[xyz,1]-statistics1[moddee,1])/nconsum
      y.post.lim<-c(0,max(zyx[xyz]))
    }else{
      
      zyx[xyz]<-statistics1[xyz,3]/norm1
      y.post.lim<-c(0,NaN)
    }
    
    
    if(is.nan(y.post.lim[2]))
      y.post.lim[2]<-max(statistics1[,3]/norm1,na.rm = TRUE)
    if(is.nan(y.post.lim[2]))
      y.post.lim[2]<-1
    
    if(crit$mlik)
    {
      
      if(printable.opt)print(paste("drawing ",workdir,template,"_Pr(M|D).jpg",sep = ""))
      capture.output({withRestarts(tryCatch(capture.output({
        jpeg(file=paste(workdir,template,"_Pr(M|D).jpg",sep = ""))
        plot(xlab = "model_id", ylab="Pr(M(model_id)|D)",ylim = y.post.lim, statistics1[,4]/norm1,pch=19, col = 7,cex= 3*((statistics1[,9]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
        points(statistics1[,4]/norm1,pch=8,  col = ifelse(statistics1[,4]>0,5,0),cex= ifelse(statistics1[,4]>0,statistics1[,4]/norm+1,0))
        points(statistics1[,4]/norm1,pch=2,  col = ifelse(statistics1[,10]>0,2,0),cex= ifelse(statistics1[,10]>0,statistics1[,10]/norm+1,0))
        points(statistics1[,4]/norm1,pch=3,  col = ifelse(statistics1[,11]>0,3,0),cex= ifelse(statistics1[,11]>0,statistics1[,11]/norm+1,0))
        points(statistics1[,4]/norm1,pch=4,  col = ifelse(statistics1[,12]>0,4,0),cex= ifelse(statistics1[,12]>0,statistics1[,12]/norm+1,0))
        points(statistics1[,4]/norm1,pch=6,  col = ifelse(statistics1[,13]>0,6,0),cex= ifelse(statistics1[,13]>0,statistics1[,13]/norm+1,0))
        points(statistics1[,4]/norm1,pch=1,  col = ifelse(statistics1[,14]>0,1,0),cex= ifelse(statistics1[,14]>0,statistics1[,14]/norm+1,0))
        points(zyx,pch=20,  col = 2,cex= 2)
        dev.off()
      })), abort = function(){onerr<-TRUE})})
    }
    
    if(crit$waic && crit$mlik)
    {
      
      if(printable.opt)print(paste("drawing ",workdir,template,"_waic-Pr(M|D).jpg",sep = ""))
      capture.output({withRestarts(tryCatch(capture.output({
        jpeg(file=paste(workdir,template,"_waic-Pr(M|D).jpg",sep = ""))
        plot(xlim = waic.lim, ylim = y.post.lim,xlab = "WAIC(M)", ylab="Pr(M|D)",statistics1[,2],statistics1[,4]/norm1, pch=19, col = 7,cex= 3*((statistics1[,9]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
        points(statistics1[,2],statistics1[,4]/norm1,pch=8,  col = ifelse(statistics1[,4]>0,5,0),cex= ifelse(statistics1[,4]>0,statistics1[,4]/norm+1,0))
        points(statistics1[,2],statistics1[,4]/norm1,pch=2,  col = ifelse(statistics1[,10]>0,2,0),cex= ifelse(statistics1[,10]>0,statistics1[,10]/norm+1,0))
        points(statistics1[,2],statistics1[,4]/norm1,pch=3,  col = ifelse(statistics1[,11]>0,3,0),cex= ifelse(statistics1[,11]>0,statistics1[,11]/norm+1,0))
        points(statistics1[,2],statistics1[,4]/norm1,pch=4,  col = ifelse(statistics1[,12]>0,4,0),cex= ifelse(statistics1[,12]>0,statistics1[,12]/norm+1,0))
        points(statistics1[,2],statistics1[,4]/norm1,pch=6,  col = ifelse(statistics1[,13]>0,6,0),cex= ifelse(statistics1[,13]>0,statistics1[,13]/norm+1,0))
        points(statistics1[,2],statistics1[,4]/norm1,pch=1,  col = ifelse(statistics1[,14]>0,1,0),cex= ifelse(statistics1[,14]>0,statistics1[,14]/norm+1,0))
        points(statistics1[,2],zyx,pch=20,  col = 2,cex= 2)
        dev.off()
      })), abort = function(){onerr<-TRUE})})
    }
    
    if(crit$dic && crit$mlik)
    {
      
      if(printable.opt)print(paste("drawing ",workdir,template,"_dic-Pr(M|D).jpg",sep = ""))
      capture.output({withRestarts(tryCatch(capture.output({
        jpeg(file=paste(workdir,template,"_dic-Pr(M|D).jpg",sep = ""))
        plot(xlim = dic.lim, ylim = y.post.lim,xlab = "dic(M)", ylab="Pr(M|D)",statistics1[,3],statistics1[,4]/norm1, pch=19, col = 7,cex= 3*((statistics1[,9]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
        points(statistics1[,3],statistics1[,4]/norm1,pch=8,  col = ifelse(statistics1[,4]>0,5,0),cex= ifelse(statistics1[,4]>0,statistics1[,4]/norm+1,0))
        points(statistics1[,3],statistics1[,4]/norm1,pch=2,  col = ifelse(statistics1[,10]>0,2,0),cex= ifelse(statistics1[,10]>0,statistics1[,10]/norm+1,0))
        points(statistics1[,3],statistics1[,4]/norm1,pch=3,  col = ifelse(statistics1[,11]>0,3,0),cex= ifelse(statistics1[,11]>0,statistics1[,11]/norm+1,0))
        points(statistics1[,3],statistics1[,4]/norm1,pch=4,  col = ifelse(statistics1[,12]>0,4,0),cex= ifelse(statistics1[,12]>0,statistics1[,12]/norm+1,0))
        points(statistics1[,3],statistics1[,4]/norm1,pch=6,  col = ifelse(statistics1[,13]>0,6,0),cex= ifelse(statistics1[,13]>0,statistics1[,13]/norm+1,0))
        points(statistics1[,3],statistics1[,4]/norm1,pch=1,  col = ifelse(statistics1[,14]>0,1,0),cex= ifelse(statistics1[,14]>0,statistics1[,14]/norm+1,0))
        points(statistics1[,3],zyx,pch=20,  col = 2,cex= 2)
        dev.off()
      })), abort = function(){onerr<-TRUE})})
    }
    
    
    if(crit$mlik)
    {
      
      jpeg(file=paste(workdir,template,"_mlik-Pr(M|D).jpg",sep = ""))
      capture.output({withRestarts(tryCatch(capture.output({
        plot(xlim = mlik.lim,ylim = y.post.lim, xlab = "MLIK", ylab="Pr(M|D)",statistics1[,1],statistics1[,4]/norm1, pch=19, col = 7,cex= 3*((statistics1[,9]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
        points(statistics1[,1],statistics1[,4]/norm1,pch=8,  col = ifelse(statistics1[,4]>0,5,0),cex= ifelse(statistics1[,4]>0,statistics1[,4]/statistics1[,3]*3 +1,0))
        points(statistics1[,1],statistics1[,4]/norm1,pch=2,  col = ifelse(statistics1[,10]>0,2,0),cex= ifelse(statistics1[,10]>0,statistics1[,10]/statistics1[,4]*3 +1,0))
        points(statistics1[,1],statistics1[,4]/norm1,pch=3,  col = ifelse(statistics1[,11]>0,3,0),cex= ifelse(statistics1[,11]>0,statistics1[,11]/statistics1[,4]*3 +1,0))
        points(statistics1[,1],statistics1[,4]/norm1,pch=4,  col = ifelse(statistics1[,12]>0,4,0),cex= ifelse(statistics1[,12]>0,statistics1[,12]/statistics1[,4]*3 +1,0))
        points(statistics1[,1],statistics1[,4]/norm1,pch=6,  col = ifelse(statistics1[,13]>0,6,0),cex= ifelse(statistics1[,13]>0,statistics1[,13]/statistics1[,4]*3 +1,0))
        points(statistics1[,1],statistics1[,4]/norm1,pch=1,  col = ifelse(statistics1[,14]>0,1,0),cex= ifelse(statistics1[,14]>0,statistics1[,14]/statistics1[,4]*3 +1,0))
        points(statistics1[,1],zyx,pch=20,  col = 2,cex= 2)
        dev.off()
      })), abort = function(){onerr<-TRUE})})
    }
    
    if(draw_dist)
    {   
      if(printable.opt)print("Calculating distance matrix, may take a significant amount of time, may also produce errors if your machine does not have enough memory")
      capture.output({withRestarts(tryCatch(capture.output({
        lldd<-2^(Nvars)+1
        
        moddee<-which(zyx==max(zyx,na.rm = TRUE))
        iidr<-which(statistics1[moddee,1]==max(statistics1[moddee,1],na.rm = TRUE))
        iidr<-which(statistics1[moddee[iidr],3]==max(statistics1[moddee[iidr],3],na.rm = TRUE))
        moddee<-moddee[iidr]
        if(length(moddee)>1)
          moddee<-moddee[1]
        
        vec<-dectobit(moddee-1)
        varcur<-c(array(0,dim = (Nvars -length(vec))),vec)
        df = data.frame(varcur)
        
        
        for(i in 1:(lldd-1))
        {
          if(i==moddee)
          {
            next
            
          }else
          {
            vec<-dectobit(i-1)
            varcur<-c(array(0,dim = (Nvars -length(vec))),vec)
            df<-cbind(df,varcur)
            #colnames(x = df)[i] <- paste("solution ",i)
          }
        }
        df<-t(df)
        
        
        x<-dist(x = df,method = "binary")
        
        dists<-c(0,x[1:lldd-1])
        
        #length(dists)
        #which(dists==0)
      })), abort = function(){onerr<-TRUE})})
      
      if(crit$mlik)
      {
        
        if(printable.opt)print(paste("drawing ",workdir,template,"_distance-mlik.jpg",sep = ""))
        capture.output({withRestarts(tryCatch(capture.output({
          jpeg(file=paste(workdir,template,"_distance-MLIK",sep = ""))
          plot(xlab = "x = |M* - M| = distance from the main mode", ylab="MLIK(M)",ylim = mlik.lim,y=c(statistics1[moddee,1],statistics1[-moddee,1]),x=dists,pch=19, col = 7,cex= 1*((statistics1[,9]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
          points(y=c(statistics1[moddee,1],statistics1[-moddee,1]),x=dists,pch=8,  col = ifelse(c(statistics1[moddee,4],statistics1[-moddee,4])>0,5,0),cex= ifelse(c(statistics1[moddee,4],statistics1[-moddee,4])>0,c(statistics1[moddee,4],statistics1[-moddee,4])/norm+1,0))
          points(y=c(statistics1[moddee,1],statistics1[-moddee,1]),x=dists,pch=2,  col = ifelse(c(statistics1[moddee,10],statistics1[-moddee,10])>0,2,0),cex= ifelse(c(statistics1[moddee,10],statistics1[-moddee,10])>0,c(statistics1[moddee,10],statistics1[-moddee,10])/norm+1,0))
          points(y=c(statistics1[moddee,1],statistics1[-moddee,1]),x=dists,pch=3,  col = ifelse(c(statistics1[moddee,11],statistics1[-moddee,11])>0,3,0),cex= ifelse(c(statistics1[moddee,11],statistics1[-moddee,11])>0,c(statistics1[moddee,11],statistics1[-moddee,11])/norm+1,0))
          points(y=c(statistics1[moddee,1],statistics1[-moddee,1]),x=dists,pch=4,  col = ifelse(c(statistics1[moddee,12],statistics1[-moddee,12])>0,4,0),cex= ifelse(c(statistics1[moddee,12],statistics1[-moddee,12])>0,c(statistics1[moddee,12],statistics1[-moddee,12])/norm+1,0))
          points(y=c(statistics1[moddee,1],statistics1[-moddee,1]),x=dists,pch=6,  col = ifelse(c(statistics1[moddee,13],statistics1[-moddee,13])>0,6,0),cex= ifelse(c(statistics1[moddee,13],statistics1[-moddee,13])>0,c(statistics1[moddee,13],statistics1[-moddee,13])/norm+1,0))
          points(y=c(statistics1[moddee,1],statistics1[-moddee,1]),x=dists,pch=1,  col = ifelse(c(statistics1[moddee,14],statistics1[-moddee,14])>0,1,0),cex= ifelse(c(statistics1[moddee,14],statistics1[-moddee,14])>0,c(statistics1[moddee,14],statistics1[-moddee,14])/norm+1,0))
          dev.off()
        })), abort = function(){onerr<-TRUE})})
      }
      if(crit$waic)
      {
        
        if(printable.opt)print(paste("drawing ",workdir,template,"_distance-waic.jpg",sep = ""))
        capture.output({withRestarts(tryCatch(capture.output({
          jpeg(file=paste(workdir,template,"_distance-waic",sep = ""))
          plot(xlab = "x = |M* - M| = distance from the main mode", ylab="WAIC(M)",ylim = waic.lim,y=c(statistics1[moddee,2],statistics1[-moddee,2]),x=dists,pch=19, col = 7,cex= 1*((statistics1[,9]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
          points(y=c(statistics1[moddee,2],statistics1[-moddee,2]),x=dists,pch=8,  col = ifelse(c(statistics1[moddee,4],statistics1[-moddee,4])>0,5,0),cex= ifelse(c(statistics1[moddee,4],statistics1[-moddee,4])>0,c(statistics1[moddee,4],statistics1[-moddee,4])/norm+1,0))
          points(y=c(statistics1[moddee,2],statistics1[-moddee,2]),x=dists,pch=2,  col = ifelse(c(statistics1[moddee,10],statistics1[-moddee,10])>0,2,0),cex= ifelse(c(statistics1[moddee,10],statistics1[-moddee,10])>0,c(statistics1[moddee,10],statistics1[-moddee,10])/norm+1,0))
          points(y=c(statistics1[moddee,2],statistics1[-moddee,2]),x=dists,pch=3,  col = ifelse(c(statistics1[moddee,11],statistics1[-moddee,11])>0,3,0),cex= ifelse(c(statistics1[moddee,11],statistics1[-moddee,11])>0,c(statistics1[moddee,11],statistics1[-moddee,11])/norm+1,0))
          points(y=c(statistics1[moddee,2],statistics1[-moddee,2]),x=dists,pch=4,  col = ifelse(c(statistics1[moddee,12],statistics1[-moddee,12])>0,4,0),cex= ifelse(c(statistics1[moddee,12],statistics1[-moddee,12])>0,c(statistics1[moddee,12],statistics1[-moddee,12])/norm+1,0))
          points(y=c(statistics1[moddee,2],statistics1[-moddee,2]),x=dists,pch=6,  col = ifelse(c(statistics1[moddee,13],statistics1[-moddee,13])>0,6,0),cex= ifelse(c(statistics1[moddee,13],statistics1[-moddee,13])>0,c(statistics1[moddee,13],statistics1[-moddee,13])/norm+1,0))
          points(y=c(statistics1[moddee,2],statistics1[-moddee,2]),x=dists,pch=1,  col = ifelse(c(statistics1[moddee,14],statistics1[-moddee,14])>0,1,0),cex= ifelse(c(statistics1[moddee,14],statistics1[-moddee,14])>0,c(statistics1[moddee,14],statistics1[-moddee,14])/norm+1,0))
          dev.off()
        })), abort = function(){onerr<-TRUE})})
      }
      if(crit$dic)
      {
        
        if(printable.opt)print(paste("drawing ",workdir,template,"_distance-dic.jpg",sep = ""))
        capture.output({withRestarts(tryCatch(capture.output({
          jpeg(file=paste(workdir,template,"_distance-dic",sep = ""))
          plot(xlab = "x = |M* - M| = distance from the main mode", ylab="DIC(M)",ylim = dic.lim,y=c(statistics1[moddee,3],statistics1[-moddee,3]),x=dists,pch=19, col = 7,cex= 1*((statistics1[,9]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
          points(y=c(statistics1[moddee,3],statistics1[-moddee,3]),x=dists,pch=8,  col = ifelse(c(statistics1[moddee,4],statistics1[-moddee,4])>0,5,0),cex= ifelse(c(statistics1[moddee,4],statistics1[-moddee,4])>0,c(statistics1[moddee,4],statistics1[-moddee,4])/norm+1,0))
          points(y=c(statistics1[moddee,3],statistics1[-moddee,3]),x=dists,pch=2,  col = ifelse(c(statistics1[moddee,10],statistics1[-moddee,10])>0,2,0),cex= ifelse(c(statistics1[moddee,10],statistics1[-moddee,10])>0,c(statistics1[moddee,10],statistics1[-moddee,10])/norm+1,0))
          points(y=c(statistics1[moddee,3],statistics1[-moddee,3]),x=dists,pch=3,  col = ifelse(c(statistics1[moddee,11],statistics1[-moddee,11])>0,3,0),cex= ifelse(c(statistics1[moddee,11],statistics1[-moddee,11])>0,c(statistics1[moddee,11],statistics1[-moddee,11])/norm+1,0))
          points(y=c(statistics1[moddee,3],statistics1[-moddee,3]),x=dists,pch=4,  col = ifelse(c(statistics1[moddee,12],statistics1[-moddee,12])>0,4,0),cex= ifelse(c(statistics1[moddee,12],statistics1[-moddee,12])>0,c(statistics1[moddee,12],statistics1[-moddee,12])/norm+1,0))
          points(y=c(statistics1[moddee,3],statistics1[-moddee,3]),x=dists,pch=6,  col = ifelse(c(statistics1[moddee,13],statistics1[-moddee,13])>0,6,0),cex= ifelse(c(statistics1[moddee,13],statistics1[-moddee,13])>0,c(statistics1[moddee,13],statistics1[-moddee,13])/norm+1,0))
          points(y=c(statistics1[moddee,3],statistics1[-moddee,3]),x=dists,pch=1,  col = ifelse(c(statistics1[moddee,14],statistics1[-moddee,14])>0,1,0),cex= ifelse(c(statistics1[moddee,14],statistics1[-moddee,14])>0,c(statistics1[moddee,14],statistics1[-moddee,14])/norm+1,0))
          dev.off()
        })), abort = function(){onerr<-TRUE})})
      }
      
      if(crit$mlik)
      {
        
        if(printable.opt)print(paste("drawing ",workdir,template,"_distance-Pr(M|D).jpg",sep = ""))
        capture.output({withRestarts(tryCatch(capture.output({
          jpeg(file=paste(workdir,template,"_distance-Pr(M|D).jpg",sep = ""))
          plot(xlab = "x = |M* - M| = distance from the main mode",ylim = y.post.lim, ylab="Pr(M|D)",y=c(statistics1[moddee,4],statistics1[-moddee,4])/norm1,x=dists,pch=19, col = 7,cex= 3*((statistics1[,9]+statistics1[,5]+statistics1[,6]+statistics1[,7]+statistics1[,8])>0))
          points(y=c(statistics1[moddee,4],statistics1[-moddee,4])/norm1,x=dists,pch=8,  col = ifelse(c(statistics1[moddee,4],statistics1[-moddee,4])>0,5,0),cex= ifelse(c(statistics1[moddee,4],statistics1[-moddee,4])>0,c(statistics1[moddee,4],statistics1[-moddee,4])/norm+1,0))
          points(y=c(statistics1[moddee,4],statistics1[-moddee,4])/norm1,x=dists,pch=2,  col = ifelse(c(statistics1[moddee,10],statistics1[-moddee,10])>0,2,0),cex= ifelse(c(statistics1[moddee,10],statistics1[-moddee,10])>0,c(statistics1[moddee,10],statistics1[-moddee,10])/norm+1,0))
          points(y=c(statistics1[moddee,4],statistics1[-moddee,4])/norm1,x=dists,pch=3,  col = ifelse(c(statistics1[moddee,11],statistics1[-moddee,11])>0,3,0),cex= ifelse(c(statistics1[moddee,11],statistics1[-moddee,11])>0,c(statistics1[moddee,11],statistics1[-moddee,11])/norm+1,0))
          points(y=c(statistics1[moddee,4],statistics1[-moddee,4])/norm1,x=dists,pch=4,  col = ifelse(c(statistics1[moddee,12],statistics1[-moddee,12])>0,4,0),cex= ifelse(c(statistics1[moddee,12],statistics1[-moddee,12])>0,c(statistics1[moddee,12],statistics1[-moddee,12])/norm+1,0))
          points(y=c(statistics1[moddee,4],statistics1[-moddee,4])/norm1,x=dists,pch=6,  col = ifelse(c(statistics1[moddee,13],statistics1[-moddee,13])>0,6,0),cex= ifelse(c(statistics1[moddee,13],statistics1[-moddee,13])>0,c(statistics1[moddee,13],statistics1[-moddee,13])/norm+1,0))
          points(y=c(statistics1[moddee,4],statistics1[-moddee,4])/norm1,x=dists,pch=1,  col = ifelse(c(statistics1[moddee,14],statistics1[-moddee,14])>0,1,0),cex= ifelse(c(statistics1[moddee,14],statistics1[-moddee,14])>0,c(statistics1[moddee,14],statistics1[-moddee,14])/norm+1,0))
          points(y= c(zyx[moddee],zyx[-moddee]),x=dists,pch=20,  col = 2,cex= 2)
          dev.off()
        })), abort = function(){onerr<-TRUE})})
      }
    }
    if(crit$mlik)
    {
      if(printable.opt)print(paste("drawing ",workdir,template,"_mds-Pr(M|D).jpg",sep = ""))
      if(printable.opt)print("Calculating distance matrix, may take a significant amount of time, may also produce errors if your machine does not have enough memory")
      if(printable.opt)print("Calculating mds coordinates, may take a significant amount of time, may also produce errors if your machine does not have enough memory")
      capture.output({withRestarts(tryCatch(capture.output({
        # further address subset of the set of the best solution of cardinality 1024  
        if(lldd>mds_size)    
        {
          lldd<-mds_size
          quant<-(sort(statistics1[,1],decreasing = TRUE)[lldd+1])
          indmds<-which(statistics1[,1]>quant)
          length(indmds)
          
        }else{
          quant<- -Inf
          indmds<-1:(lldd)
        }
        
        
        
        
        
        vec<-dectobit(moddee-1)
        varcur<-c(array(0,dim = (Nvars -length(vec))),vec)
        df = data.frame(varcur)
        
        
        for(i in 1:(lldd-1))
        {
          if(i==moddee)
          {
            next
            
          }else
          {
            vec<-dectobit(indmds[i]-1)
            varcur<-c(array(0,dim = (Nvars -length(vec))),vec)
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
        jpeg(file=paste(workdir,template,"_mds_map_posteriors.jpg",sep = ""))
        plot(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])>0,c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])/norm1*50,0))
        points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=2,pch = 10,cex= c(zyx[moddee],zyx[setdiff(indmds, moddee)])*50,0)
        points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=2,pch = 10,cex= c(zyx[moddee],zyx[setdiff(indmds, moddee)])*50,0)
        points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])>0,c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])/norm1*50,0))
        points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])>0,c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])/norm1*50,0))
        points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])>0,c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])/norm1*50,0))
        points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])>0,c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])/norm1*50,0))
        points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])>0,c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])/norm1*50,0))
        points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])>0,c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])/norm1*50,0))
        points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])>0,c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])/norm1*50,0))
        points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])>0,c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])/norm1*50,0))
        points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])>0,c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])/norm1*50,0))
        points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])>0,c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])/norm1*50,0))
        points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=5,pch = 8,cex= ifelse(c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])>0,c(statistics1[moddee,4],statistics1[setdiff(indmds, moddee),4])/norm1*50,0))
        points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=7,pch = 19,cex= 0.4)
        points(x.mds[], y.mds[], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",col=1,pch = 19,cex= 0.01)
        dev.off() 
      })), abort = function(){onerr<-TRUE})})
    }
  },
  #calculates posterior probabilities based on a current search
  post_proceed_results = function(statistics1)
  {
    
    xyz<-which(!is.na(statistics1[,1]))
    g.results[4,2] <- length(xyz)
    xyz<-intersect(xyz,which(statistics1[,1]!=-10000))
    moddee<-which(statistics1[,1]==max(statistics1[,1],na.rm = TRUE))[1]
    zyx<-array(data = NA,dim = 2 ^(Nvars)+1)
    nconsum<-sum(exp(-statistics1[moddee,1]+statistics1[xyz,1]),na.rm = TRUE)
    if( nconsum > 0)
    {
      zyx[xyz]<-exp(statistics1[xyz,1]-statistics1[moddee,1])/nconsum
      y.post.lim<-c(0,max(zyx[xyz]))
    }else{
      
      zyx[xyz]<-statistics1[xyz,3]/norm1
      y.post.lim<-c(0,NaN)
    }
    statistics1[,15]<-zyx
    
    lldd<-2^(Nvars)+1
    p.post<-array(data = 0,dim = Nvars)
    for(i in 1:(lldd))
    {
      if(is.na(statistics1[i,15]))
        next
      vec<-dectobit(i-1)
      varcur<-c(array(0,dim = (Nvars -length(vec))),vec)
      p.post <- (p.post + varcur*statistics1[i,15])
      
    }
    
    return(list(p.post = p.post, m.post = zyx, s.mass = sum(exp(statistics1[xyz,1]),na.rm = TRUE)))
  },
  calculate_quality_measures = function(vect,n,truth)
  {
    rmse.pi <- array(data = 0,Nvars)
    bias.pi <- array(data = 0,Nvars)
    for(i in 1:n)
    {
      bias.pi <- (bias.pi + (vect[,i]-truth))
      rmse.pi <- (rmse.pi + (vect[,i]^2+truth^2 - 2*vect[,i]*truth))
    }
    bias.pi <- bias.pi/n
    rmse.pi <- rmse.pi/n
    return(list(bias.pi = bias.pi, rmse.pi = rmse.pi))
  }
  
)

)



# set up the parameters of the simulation or optimization

#define your working directory
workdir<-"/mn/anatu/ansatte-u3/aliaksah/Desktop/Untitled Folder/"
#prepare data
simx <- read.table(paste(workdir,"simcen-x.txt",sep = "",collapse = ""))
simy <- read.table(paste(workdir,"simcen-y.txt",sep = "",collapse = ""))
data.example <- cbind(simy,simx)
names(data.example)[1]="Y"

#fparam <- c("Const",colnames(data)[-1])
fparam.example <- colnames(data.example)[-1]
fobserved.example <- colnames(data.example)[1]

#dataframe for results; n/b +1 is required for the summary statistics
statistics1 <- big.matrix(nrow = 2 ^(length(fparam.example))+1, ncol = 15,init = NA, type = "double")
statistics <- describe(statistics1)

#create MySearch object with default parameters
mySearch = EMJMCMC2016()


# load functions as in BAS article by Clyde, Ghosh and Littman to reproduce their first example
mySearch$estimator = estimate.bas.lm
mySearch$estimator.args = list(data = data.example,prior = 3, g = 100 ,n=100)


#play around with various methods in order to get used to them and see how they work

system.time(
  zzz<-mySearch$learnlocalMCMC(list(varcur=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),mlikcur=-10000,waiccur =10000,statid=5,varold=rep(0,length(fparam.example)),reverse=TRUE))
)


system.time(
  sss<-mySearch$learnlocalSA(list(varcur=rbinom(n = length(fparam.example),size = 1,prob = 0.5),mlikcur=-Inf,waiccur =Inf,statid=5,varold=rep(0,length(fparam.example)),switch.type=5,objcur = -Inf,objold = -Inf, sa2 = FALSE ,reverse=FALSE))
)

mySearch$M.nd = as.integer(1000)

system.time(
  ggg<-mySearch$learnlocalND(list(varcur=rep(1,length(fparam.example)),mlikcur=-Inf,waiccur =Inf,statid=6,varold=rep(1,length(fparam.example)),switch.type=6,objcur = -Inf,objold = -Inf, reverse=TRUE))
)
system.time(
  ggg<-mySearch$learnlocalND(list(varcur=rep(1,length(fparam.example)),mlikcur=-Inf,waiccur =Inf,statid=6,varold=rep(1,length(fparam.example)),switch.type=6,objcur = -Inf,objold = -Inf, reverse=FALSE))
)

system.time(
  mmm<-mySearch$modejumping_mcmc(list(varcur=rbinom(n = length(fparam.example),size = 1,prob = 0.5),statid=5, distrib_of_proposals = c(35,35,35,35,100),distrib_of_neighbourhoods=array(data = 0.5,dim = c(5,3)), eps = 0.0001, trit = 550, burnin = 50, max.time = 30, maxit = 550, print.freq =1))
)

system.time(
  fff<-mySearch$forward_selection(list(varcur=rep(0,length(fparam.example)),mlikcur=-Inf,waiccur =Inf,locstop = FALSE, statid=5))
)

system.time(
  bbb<-mySearch$backward_selection(list(varcur=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),mlikcur=-Inf,waiccur =Inf,locstop = FALSE,statid=5))
)

system.time(
  gamba<-mySearch$forw_backw_walk(list(steps=10,p1=0.5,reverse = TRUE))
)

# this one MUST be completed before moving to the experiments
system.time(
  FFF<-mySearch$full_selection(list(statid=6, totalit =32769, ub = 36,mlikcur=-Inf,waiccur =100000))
)
# check that all models are enumerated during the full search procedure
idn<-which(!is.na(statistics1[,1]))
length(idn)



# see the best current results and total number of iterations
mySearch$g.results[,]
#compare to
which(statistics1[,1]==max(statistics1[,1],na.rm = TRUE))

View(statistics1[2113,])

which(statistics1[868,1]<= -10)

# get graphical output (only makes sence after EMJMCMC procedure carried out)
mySearch$visualize_results(statistics1 = statistics1, template = "test/full1BASclasstest", mds_size = 1024, crit = list(waic=TRUE,mlik=TRUE,dic=TRUE), draw_dist = TRUE)



# once full search is completed, get the truth for the experiment
ppp<-mySearch$post_proceed_results(statistics1 = statistics1)
truth = ppp$p.post # make sure it is equal to Truth column from the article
truth.m = ppp$m.post
truth.prob = ppp$s.mass
ordering = sort(ppp$p.post,index.return=T)
fake500 <- sum(exp(x = sort(statistics1[,1],decreasing = T)[1:700]),na.rm = TRUE)/truth.prob


print("pi truth")
sprintf("%.10f",truth[ordering$ix])

# vecbuf<-vect
# freqs.buf<- freqs
# freqs.p.buf <- freqs.p

#reproduce the 1st experiment as in BAS article


mySearch$switch.type=as.integer(5)
mySearch$switch.type.glob=as.integer(1)
mySearch$max.N.glob=as.integer(6)
mySearch$min.N.glob=as.integer(5)
mySearch$max.N=as.integer(12)
mySearch$min.N=as.integer(5)
mySearch$recalc.margin = as.integer(2^15)
mySearch$max.cpu=as.integer(4)
mySearch$p.add = array(data = 0.5,dim = 15)
mySearch$printable.opt = T
mySearch$p.add = array(data = 0.5,dim = 15)
fff<-mySearch$forward_selection(list(varcur=rep(0,length(fparam.example)),mlikcur=-Inf,waiccur =Inf,locstop = FALSE,statid=-1))
bbb<-mySearch$backward_selection(list(varcur=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),mlikcur=-Inf,waiccur =Inf,locstop = FALSE,statid=-1))


cor(data.example)[1,]

resss<-0
masss<-0
mmm<-0

statistics1 <- big.matrix(nrow = 2 ^(length(fparam.example))+1, ncol = 15,init = NA, type = "double")
statistics <- describe(statistics1)
mySearch$g.results[4,1]<-0
mySearch$g.results[4,2]<-0
distrib_of_neighbourhoods=array(data = 2,dim = c(5,7))
distrib_of_neighbourhoods[,1]<-8
distrib_of_neighbourhoods[,5]<-0
distrib_of_neighbourhoods[,6]<-0

while(mmm<9900)
{

  initsol=rbinom(n = length(fparam.example),size = 1,prob = 0.5)
  mySearch$p.add = array(data = 0.5,dim = 15)
  resm<-mySearch$modejumping_mcmc(list(varcur=initsol,statid=5, distrib_of_proposals = c(0,0,50,0,5500),distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = 0.0001, trit = 3200, trest = 32000 , burnin = 30, max.time = 30, maxit = 100000, print.freq =5))
  resss<-c(resss,mySearch$g.results[4,1])
  mySearch$g.results[4,2]
  mySearch$g.results[4,1]
  resm$bayes.results$s.mass/truth.prob
  mmm<-mySearch$post_proceed_results(statistics1 = statistics1)$s.mass/truth.prob*10000
  masss<-c(masss,mmm)
  plot(resss,ylim = c(0,20000),type = "l", main = "EMJMCMC", ylab = "iterations, mass x 10^3",xlab = "test points (every 50 iterations of EMJCMCMC)" )
  lines(masss,col = 3)
  #   set.seed(length(resss)*20000)
  #   mySearch$seed=as.integer(runif(n = 1,min = length(resss),max = length(resss)*20000))
  #   mySearch$forw_backw_walk(list(steps=2,p1=0.5,reverse = FALSE))
  #   resss<-c(resss,mySearch$g.results[4,2])
  #   mmm<-mySearch$post_proceed_results(statistics1 = statistics1)$s.mass/truth.prob*10000
  #   masss<-c(masss,mmm)
  #   plot(resss,ylim = c(0,10000))
  #   points(masss,col = 4)    
  #   set.seed(length(resss)*10000)
}

mySearch$switch.type=as.integer(1)
mySearch$switch.type.glob=as.integer(1)
mySearch$printable.opt = TRUE


mySearch$max.N.glob=as.integer(4)
mySearch$min.N.glob=as.integer(4)
mySearch$max.N=as.integer(2)
mySearch$min.N=as.integer(2)
mySearch$recalc.margin = as.integer(400)
distrib_of_proposals = c(76.91870,71.25264,87.68184,60.55921,17812.39852)
distrib_of_neighbourhoods=t(array(data = c(7.6651604,16.773326,14.541629,12.839445,2.964227,13.048343,7.165434,
                                           0.9936905,15.942490,11.040131,3.200394,15.349051,5.466632,14.676458,
                                           1.5184551,9.285762,6.125034,3.627547,13.343413,2.923767,15.318774,
                                           14.5295380,1.521960,11.804457,5.070282,6.934380,10.578945,12.455602,
                                           6.0826035,2.453729,14.340435,14.863495,1.028312,12.685017,13.806295),dim = c(7,5)))

Niter <- 100
thining<-1
system.time({
  
  vect <-array(data = 0,dim = c(15,Niter))
  vect.mc <-array(data = 0,dim = c(15,Niter))
  
  inits <-array(data = 0,dim = Niter)
  freqs <-array(data = 100,dim = c(5,Niter))
  freqs.p <-array(data = 100,dim = c(5,7,Niter))
  masses <- array(data = 0,dim = Niter)
  biases.m <- array(data = 0,dim = 2 ^(length(fparam.example))+1)
  biases.m.mc <- array(data = 0,dim = 2 ^(length(fparam.example))+1)
  rmse.m <- array(data = 0,dim = Niter)
  rmse.m.mc <- array(data = 0,dim = Niter)
  iterats <- array(data = 0,dim = c(2,Niter))
  
  for(i in 1:Niter)
  {
    statistics1 <- big.matrix(nrow = 2 ^(length(fparam.example))+1, ncol = 15,init = NA, type = "double")
    statistics <- describe(statistics1)
    mySearch$g.results[4,1]<-0
    mySearch$g.results[4,2]<-0
    mySearch$p.add = array(data = 0.5,dim = 15)
    #distrib_of_neighbourhoods=array(data = runif(n = 5*7,min = 0, max = 20),dim = c(5,7))
    #distrib_of_proposals = runif(n = 5,min = 0, max = 100)
    #distrib_of_proposals[5]=sum(distrib_of_proposals[1:4])*runif(n = 1,min = 50, max = 150)
    print("BEGIN ITERATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print(i)
    set.seed(10*i)
    initsol=rbinom(n = length(fparam.example),size = 1,prob = 0.8)
    inits[i] <- mySearch$bittodec(initsol)
    freqs[,i]<- distrib_of_proposals
    resm<-mySearch$modejumping_mcmc(list(varcur=initsol,statid=5, distrib_of_proposals = distrib_of_proposals,distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = 0.0001, trit = 32700, trest = 3200 , burnin = 50, max.time = 30, maxit = 100000, print.freq =500))
    vect[,i]<-resm$bayes.results$p.post
    vect.mc[,i]<-resm$p.post
    masses[i]<-resm$bayes.results$s.mass/truth.prob
    print(masses[i])
    freqs.p[,,i] <- distrib_of_neighbourhoods
    cur.p.post <- resm$bayes.results$m.post
    cur.p.post[(which(is.na(cur.p.post)))]<-0
    rmse.m[i]<-mean((cur.p.post - truth.m)^2,na.rm = TRUE)
    biases.m<-biases.m + (cur.p.post - truth.m)
    cur.p.post.mc <- resm$m.post
    cur.p.post.mc[(which(is.na(cur.p.post.mc)))]<-0
    rmse.m.mc[i]<-mean((cur.p.post.mc - truth.m)^2,na.rm = TRUE)
    biases.m.mc<-biases.m.mc + (cur.p.post.mc - truth.m)
    iterats[1,i]<-mySearch$g.results[4,1]
    iterats[2,i]<-mySearch$g.results[4,2]
    print("COMPLETE ITERATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! with")
    print(iterats[2,i])
    remove(statistics1)
    remove(statistics)
  }
}
)



Nlim <- 1
order.deviat <- sort(masses,decreasing = TRUE,index.return=T)


print("model bias rm")
sqrt(mean((biases.m/Niter)^2,na.rm = TRUE))*100000
print("model rmse rm")
sqrt(mean(rmse.m))*100000

print("model bias mc")
sqrt(mean((biases.m.mc/Niter)^2,na.rm = TRUE))*100000
print("model rmse mc")
sqrt(mean(rmse.m.mc))*100000


print("model coverages")
mean(masses)
print("mean # of iterations")# even smaller on average than in BAS
mean(iterats[1,])
print("mean # of estimations")# even smaller on average than in BAS
mean(iterats[2,])



# correlation between the MSE and the masses, obviously almost minus 1
cor(rmse.m,masses)
cor(rmse.m.mc,masses)

truth.buf <- array(data = 0,dim = c(15,Niter))
truth.buf[,1:Niter]<-truth
bias <- vect - truth.buf
bias.mc <- vect.mc - truth.buf
rmse <- (vect^2 +truth.buf^2 - 2*vect*truth.buf)
rmse.mc <- (vect.mc^2 +truth.buf^2 - 2*vect.mc*truth.buf)
bias.avg.rm<-rowMeans(bias)
rmse.avg.rm <-sqrt(rowMeans(rmse))
bias.avg.mc<-rowMeans(bias.mc)
rmse.avg.mc <-sqrt(rowMeans(rmse.mc))


print("pi biases rm")
sprintf("%.10f",bias.avg.rm[ordering$ix]*100)
print("pi rmse rm")
sprintf("%.10f",rmse.avg.rm[ordering$ix]*100)

print("pi biases mc")
sprintf("%.10f",bias.avg.mc[ordering$ix]*100)
print("pi rmse mc")
sprintf("%.10f",rmse.avg.mc[ordering$ix]*100)

# # 
View((cbind(bias.avg.rm[ordering$ix],rmse.avg.rm[ordering$ix],bias.avg.mc[ordering$ix],rmse.avg.mc[ordering$ix])*100))


# indexs<-(sort.int(cols.b, index.return=TRUE))
# plot(y = inits[indexs$ix], x = 1:1000,col = 1, type = "p")
# points(cols.b[indexs$ix]*10000,col = 4)
# points(cols.r[indexs$ix]*50000,col = 5)
# points(cols.r[indexs$ix]*50000+cols.b[indexs$ix]*10000,col = 3)
# points(freqs.p[1,indexs$ix]*30,col= 6)
# points(freqs.p[2,indexs$ix]*30,col= 7)
# points(freqs.p[3,indexs$ix]*30,col= 8)
# points(freqs.p[4,indexs$ix]*30,col= 2)
# 

#View((proceeded$p.post-truth)[c(12,14,10,8,6,7,13,11,15,5,9,2,1,3,4)])
#View((proceeded$p.post^2+truth^2-2*proceeded$p.post*truth)[c(12,14,10,8,6,7,13,11,15,5,9,2,1,3,4)])
# cor(cols.b,cols.r)
# cor(cols.b,inits)
# cor(cols.r,inits)
cor(masses,freqs.p[1,1:243])
cor(masses,freqs.p[2,1:243])
cor(masses,freqs.p[3,1:243])
cor(masses,freqs.p[4,1:243])
cor(masses,freqs.p[5,1:243])

#View(statistics1[which(!is.na(statistics1[,1]))])

# build a package
#package.skeleton(name ="EMJMCMC", code_files =paste(workdir,"mode_jumping_package_class.r",sep = "",collapse = ""))
freqs.p[,,order.deviat$ix[1:Nlim]]
freqs[,order.deviat$ix[1:Nlim]]
iterats[1,order.deviat$ix[1:Nlim]]
iterats[2,order.deviat$ix[1:Nlim]]
inits[order.deviat$ix[1:Nlim]]

statistics1 <- big.matrix(nrow = 2 ^(length(fparam.example))+1, ncol = 15,init = NA, type = "double")
statistics <- describe(statistics1)
mySearch$g.results[4,1]<-0
mySearch$g.results[4,2]<-0
mySearch$p.add = array(data = 0.5,dim = 15)


mySearch$max.N.glob=as.integer(4)
mySearch$min.N.glob=as.integer(4)
mySearch$max.N=as.integer(2)
mySearch$min.N=as.integer(2)
mySearch$recalc.margin = as.integer(400)
distrib_of_neighbourhoods=freqs.p[,,order.deviat$ix[1]]
distrib_of_proposals = freqs[,order.deviat$ix[1]]
initsol = inits[order.deviat$ix[1]]
resm<-mySearch$modejumping_mcmc(list(varcur=mySearch$dectobit(inits[order.deviat$ix[1]]),statid=5, distrib_of_proposals = distrib_of_proposals,distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = 0.0001, trit = 3273, trest = 3272 , burnin = 3, max.time = 30, maxit = 100000, print.freq =500))
mySearch$g.results[4,2]
mySearch$g.results[4,1]
resm$bayes.results$s.mass/truth.prob

mean(iterats[2,order.deviat$ix[1:Nlim]])
