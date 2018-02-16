library(qtl)
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
        train[[paste0(var,".",lev)]]<-0
        test[[paste0(var,".",lev)]]<-0
        train[[paste0(var,".",lev)]][which(train[[var]]==levels[lev])]<-1
        test[[paste0(var,".",lev)]][which(test[[var]]==levels[lev])]<-1
      }
      train[[var]] <- NULL
      test[[var]]  <- NULL
      
    }
  }
  
  # Return cleaned data to global environment
  tr <<- train
  te <<- test
}



phen1<-read.cross(dir = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/",file="journal.pone.0004318.s0021.csv",format = "csv")

#phen2<-read.cross(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/QTL%20logic%20regression/s003.csv"),sep = ";",header = T)

phen1 <- fill.geno(phen1)
write.cross(cross = phen1,format = "csv",filestem = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/s002TDf")


phen1<-fread("/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/QTL logic regression/s002TDf.csv",sep = ",",stringsAsFactors = F)[-c(1,2,3),]
names(phen1)<-c("Y",names(phen1)[-length(names(phen1))])
tr<-phen1[-which(phen1$Y=="-"),]
te<-phen1[which(phen1$Y=="-"),]
te$Y<-"-1"
tr$Y<-as.numeric(tr$Y)
te$Y<-as.numeric(te$Y)

#create dummy variables for all of the factors


process_tr_and_te(train = tr,test =te)

tr<-as.data.frame(tr)
te<-as.data.frame(te)
rm(phen1)
gc()


system.time({
  
  formula1 = as.formula(paste(colnames(tr)[1],"~ 1 +",paste0(colnames(tr)[-1],collapse = "+")))
  
  res = runemjmcmc(formula = formula1,data = X1,estimator = estimate.logic.lm,estimator.args =  list(data = data.example,n = 100, m = 50),recalc_margin = 250, save.beta = F,interact = T,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=1,mutation_rate = 300,last.mutation = 5000, max.tree.size = 4, Nvars.max =15,p.allow.replace=0.9,p.allow.tree=0.2,p.nor=0,p.and = 0.9),n.models = 20000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,print.freq = 100,advanced.param = list(
    max.N.glob=as.integer(10),
    min.N.glob=as.integer(5),
    max.N=as.integer(3),
    min.N=as.integer(1),
    printable = F))
})
