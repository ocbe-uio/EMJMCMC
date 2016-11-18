library(data.table)



get_feature_names <- function(train, test){
  var_names <- intersect(names(train), names(test))
  return(var_names)
}


process_tr_and_te <- function(train, test){
  ## Removing features with only one value
  cat("  Removing the constants features.\n")
  var_names <- get_feature_names(train, test)
  for (var in var_names) {
    if (length(unique(train[[var]])) == 1) {
      cat(var, "is constant in train. Removing it.\n")
      train[[var]] <- NULL
      test[[var]] <- NULL
    }
  }
  
  ## Convert Text variables to Numeric
  cat("  Converting text variables to numeric ids\n")
  var_names <- get_feature_names(train, test)
  for (var in var_names) {
    if (class(train[[var]])=="character" || class(train[[var]])=="factor") {
      levels <- unique(c(train[[var]], test[[var]]))
      for(lev in 1:length(levels))
      {
        train[[paste0(var,lev)]]<-0
        test[[paste0(var,lev)]]<-0
        train[[paste0(var,lev)]][which(train[[var]]==levels[lev])]<-1
        test[[paste0(var,lev)]][which(test[[var]]==levels[lev])]<-1
      }
      train[[var]] <- NULL
      test[[var]]  <- NULL
     
    }
  }
  
  # Return cleaned data to global environment
  tr <<- train
  te <<- test
}

X<-fread(input = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/Belarus Diabetus Data/genes_1.txt",sep = ";")[1:350]
Y<-fread(input = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/Belarus Diabetus Data/levels_1.txt",sep = ";")[1:350]

DAT<-cbind(Y,X)
names(DAT)[1]<-"LABL"
DAT$LABL<-DAT$LABL-1
process_tr_and_te(train = DAT,test = DAT)
tr<-as.data.frame(tr)
#build the EMJMCMC package (see the attachment)


nm<-25 # the maximal number of covariates in the model (is not omitted i.i.f. GMJMCMC is run)

gm=T #if true then GMJMCMC (Genetically modified MJMCMC) is addressed otherwise simple MJMCMC


  #now with and AIC informative prior and the mutations of the original variables allowed
  system.time({
    
    formula1 = as.formula(paste(colnames(tr)[1],"~ 1 +",paste0(colnames(tr)[-1],collapse = "+")))
    
    res = runemjmcmc(formula = formula1,data = tr,estimator = estimate.bas.glm.pen ,estimator.args =  list(data = data.example,family=binomial(),prior=aic.prior(), logn = log(200),n = 200, m = 12),recalc_margin = 250, save.beta = F,interact = gm,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=1,mutation_rate = 300,last.mutation = 2500, max.tree.size = (4-1), Nvars.max =nm,p.allow.replace=0.4,p.allow.tree=0.25,p.nor=0.2,p.and = 0.7),n.models = 5000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,print.freq = 100,advanced.param = list(
      max.N.glob=as.integer(10),
      min.N.glob=as.integer(5),
      max.N=as.integer(3),
      min.N=as.integer(1),
      printable = F))
  })
  
  # all model probabilities
  model.prob<-mySearch$post_proceed_results_hash(hashStat)$m.post
  plot(x = 1:length(model.prob),y = model.prob)
  # the highest posterior probability of some combination of genes
  max.prob = max(model.prob)
  print(max.prob)
  # it corresponds to the mlik of 
  print(max(values(hashStat)[1,]))
  # which corresponds to the model based on the following covariates 
  best.com = keys(hashStat)[which(values(hashStat)[1,]==max(values(hashStat)[1,]))]
  print(mySearch$fparam[gregexpr('1', best.com)[[1]][1:length(gregexpr('1', best.com)[[1]])]])
  
  