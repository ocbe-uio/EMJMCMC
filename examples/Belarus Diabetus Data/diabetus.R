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
process_tr_and_te(train = DAT[1:200,],test = DAT[201:350,])
tr<-as.data.frame(tr)
te<-as.data.frame(te)
#build the EMJMCMC package (see the attachment)


nm<-25 # the maximal number of covariates in the model (is not omitted i.i.f. GMJMCMC is run)

gm=T #if true then GMJMCMC (Genetically modified MJMCMC) is addressed otherwise simple MJMCMC


  #now with and AIC informative prior and the mutations of the original variables allowed
  system.time({

    formula1 = as.formula(paste(colnames(tr)[1],"~ 1 +",paste0(colnames(tr)[-1],collapse = "+")))

    res = runemjmcmc(formula = formula1,data = tr,estimator = estimate.bas.glm.pen ,estimator.args =  list(data = data.example,family=binomial(),prior=aic.prior(), logn = log(200),n = 200, m = 12),recalc_margin = 250, save.beta = T,interact = gm,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=1,mutation_rate = 300,last.mutation = 2500, max.tree.size = (4-1), Nvars.max =nm,p.allow.replace=0.4,p.allow.tree=0.25,p.nor=0.2,p.and = 0.7),n.models = 5000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,print.freq = 100,advanced.param = list(
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
  #this probability is extremely low < reasonable alpha, showing no significant models in the explored subspace
  print(max.prob)
  # it corresponds to the mlik of
  print(max(values(hashStat)[1,]))
  # which corresponds to the model based on the following covariates
  best.com = keys(hashStat)[which(values(hashStat)[1,]==max(values(hashStat)[1,]))]
  print(mySearch$fparam[gregexpr('1', best.com)[[1]][1:length(gregexpr('1', best.com)[[1]])]])
  # we can use now GLM to analyze this model
  formula.best<-as.formula(paste(colnames(tr)[1],"~ 1 +",paste0(mySearch$fparam[gregexpr('1', best.com)[[1]][1:length(gregexpr('1', best.com)[[1]])]],collapse = "+")))

  best.frequentest <- glm(formula = formula.best,family = binomial(),data = data.example)
  summary(best.frequentest)
  #deviance resiuals test
  pchisq(q = best.frequentest$deviance,df = best.frequentest$df.residual,lower.tail = F)
  #now let us perfrom some validation based on Bayesian prediction averaging, just to make sure if
  #the data contains any predictive patterns found by the algorithm


  g<-function(x)
  {
    return((x = 1/(1+exp(-x))))
  }


  Nvars<-mySearch$Nvars.max
  linx <-mySearch$Nvars.max+4
  lHash<-length(hashStat)
  mliks <- values(hashStat)[which((1:(lHash * linx)) %% linx == 1)]
  betas <- values(hashStat)[which((1:(lHash * linx)) %% linx == 4)]
  for(i in 1:(Nvars-1))
  {
    betas<-cbind(betas,values(hashStat)[which((1:(lHash * linx)) %% linx == (4+i))])
  }
  betas<-cbind(betas,values(hashStat)[which((1:(lHash * linx)) %% linx == (0))])


  system.time({

    res<-mySearch$forecast.matrix.na(link.g = g,covariates = te,betas = betas,mliks.in = mliks)$forecast

  })

  summary(res)

  length(res)
  res<-as.integer(res>=0.5)
  length(which(res>=0.5))
  length(which(res<0.5))
  length(res)

  (1-sum(abs(res-te$LABL),na.rm = T)/length(res))*100


  #FNR
  ps<-which(te$LABL==1)
  sum(abs(res[ps]-te$LABL[ps]))/(sum(abs(res[ps]-te$LABL[ps]))+length(ps))*100

  #FPR
  ns<-which(te$LABL==0)
  sum(abs(res[ns]-te$LABL[ns]))/(sum(abs(res[ns]-te$LABL[ns]))+length(ns))*100


  # compare e.g. with the following example with such patterns:



  #define your working directory, where the data files are stored
  workdir<-""

  #prepare the test set data
  simx <- read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/asteroid%20data/Recognize/NEAs.txt"),sep = ",",header = T,fill=TRUE)
  simy <-  read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/asteroid%20data/Recognize/NotNeas8%2B.txt"),sep = ",",header = T,fill=TRUE)
  simx$neo<-1
  simy$neo<-0
  #data.example <- as.data.frame(t(cbind(t(simy[sample.int(size = 400,n = 6621,replace = T),]),t(simx[sample.int(size = 600,n = 14099,replace = T),]))),stringsAsFactors = F)
  #data.example$epoch<-factor(data.example$epoch,labels = c(0,1))


  data.example <- as.data.frame(t(cbind(t(simy),t(simx))),stringsAsFactors = F)

  transform<-colnames(data.example)[-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25)]

  for(i in 1:length(transform))
  {
    print(i)
    data.example[[transform[i]]]<-as.numeric(as.character(data.example[[transform[i]]]))
  }

  data.example$esuar<-data.example$eccentricity^2
  data.example$asuar<-data.example$absolute_magnitude^2
  data.example$rsuar<-data.example$semi_major_axis^2
  data.example$rcube<-data.example$semi_major_axis^3
  data.example$anoml<-data.example$mean_anomaly*data.example$semi_major_axis
  data.example$anoms<-data.example$mean_anomaly*data.example$semi_major_axis^2
  data.example$anomv<-data.example$mean_anomaly*data.example$semi_major_axis^3
  data.example$perihell<-data.example$argument_of_perihelion*data.example$semi_major_axis
  data.example$perihels<-data.example$argument_of_perihelion*data.example$semi_major_axis^2
  data.example$perihelv<-data.example$argument_of_perihelion*data.example$semi_major_axis^3
  data.example$longitudel<-data.example$longitude_of_the_ascending.node*data.example$semi_major_axis
  data.example$longitudes<-data.example$longitude_of_the_ascending.node*data.example$semi_major_axis^2
  data.example$longitudev<-data.example$longitude_of_the_ascending.node*data.example$semi_major_axis^3


  data.example1<-data.example

  #prepare the training set data
  simx <- read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/asteroid%20data/Teach/NeoPHA.txt"),sep = ",",header = T,fill=TRUE)
  simy <-  read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/asteroid%20data/Teach/NotNeo-Type7.txt"),sep = ",",header = T,fill=TRUE)
  simx$neo<-1
  simy$neo<-0

  data.example <- as.data.frame(t(cbind(t(simy),t(simx))),stringsAsFactors = F)

  transform<-colnames(data.example)[-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25)]

  for(i in 1:length(transform))
  {
    print(i)
    data.example[[transform[i]]]<-as.numeric(as.character(data.example[[transform[i]]]))
  }

  data.example$esuar<-data.example$eccentricity^2
  data.example$asuar<-data.example$absolute_magnitude^2
  data.example$rsuar<-data.example$semi_major_axis^2
  data.example$rcube<-data.example$semi_major_axis^3
  data.example$anoml<-data.example$mean_anomaly*data.example$semi_major_axis
  data.example$anoms<-data.example$mean_anomaly*data.example$semi_major_axis^2
  data.example$anomv<-data.example$mean_anomaly*data.example$semi_major_axis^3
  data.example$perihell<-data.example$argument_of_perihelion*data.example$semi_major_axis
  data.example$perihels<-data.example$argument_of_perihelion*data.example$semi_major_axis^2
  data.example$perihelv<-data.example$argument_of_perihelion*data.example$semi_major_axis^3
  data.example$longitudel<-data.example$longitude_of_the_ascending.node*data.example$semi_major_axis
  data.example$longitudes<-data.example$longitude_of_the_ascending.node*data.example$semi_major_axis^2
  data.example$longitudev<-data.example$longitude_of_the_ascending.node*data.example$semi_major_axis^3

  #define the covariates and theobservations
  fparam.example <- colnames(data.example)[-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)]
  fobserved.example <- colnames(data.example)[1]

  # run MJMCMC
  system.time({

    formula1 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)],collapse = "+")))

    res = runemjmcmc(formula = formula1,data = data.example,estimator =estimate.bas.glm,estimator.args =  list(data = data.example,prior = aic.prior(),family = binomial(), logn = log(64)),recalc_margin = 50, save.beta = T,interact = F,relations = c("","sin","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=2,mutation_rate = 500, max.tree.size = 200000, Nvars.max = 30,p.allow.replace=0.001,p.allow.tree=0.001,p.nor=0.3,p.and = 0.7),n.models = 10000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 100,advanced.param = list(
      max.N.glob=as.integer(10),
      min.N.glob=as.integer(5),
      max.N=as.integer(3),
      min.N=as.integer(1),
      printable = F))
  })


  ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
  ppp$p.post

  mySearch$g.results[,]
  mySearch$fparam

  g<-function(x)
  {
    return((x = 1/(1+exp(-x))))
  }


  Nvars<-mySearch$Nvars.max
  linx <-mySearch$Nvars.max+4
  lHash<-length(hashStat)
  mliks <- values(hashStat)[which((1:(lHash * linx)) %% linx == 1)]
  betas <- values(hashStat)[which((1:(lHash * linx)) %% linx == 4)]
  for(i in 1:(Nvars-1))
  {
    betas<-cbind(betas,values(hashStat)[which((1:(lHash * linx)) %% linx == (4+i))])
  }
  betas<-cbind(betas,values(hashStat)[which((1:(lHash * linx)) %% linx == (0))])


  system.time({

    res<-mySearch$forecast.matrix.na(link.g = g,covariates = (data.example1[1:20720,-c(1,2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)]),betas = betas,mliks.in = mliks)$forecast

  })

  summary(res)

  length(res)
  res<-as.integer(res>=0.5)
  length(which(res>=0.5))
  length(which(res<0.5))
  length(res)
  length(which(data.example1$neo==1))

  (1-sum(abs(res-data.example1$neo),na.rm = T)/20720)*100


  #FNR
  ps<-which(data.example1$neo==1)
  sum(abs(res[ps]-data.example1$neo[ps]))/(sum(abs(res[ps]-data.example1$neo[ps]))+length(ps))*100

  #FPR
  ns<-which(data.example1$neo==0)
  sum(abs(res[ns]-data.example1$neo[ns]))/(sum(abs(res[ns]-data.example1$neo[ns]))+length(ns))*100

  # these happens when the patterns are present!
