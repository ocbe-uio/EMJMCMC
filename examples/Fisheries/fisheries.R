#Code to run the bycatch model for the positive bycatch with negative binomial likelihood conditioned on positive observations.
#
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package2.r")
library("INLA")

#Read the data
d = readRDS("dataOVT19942006.Rda")
#Use a small sample of the data for testing (to save time) 
tmp = runif(dim(d)[1])
d = d[tmp>0.9,]


#Use only the positive observations of the species selected
d = d[d$cod>0,]

#Read the range paramters previousli estimated
rangeParameters = readRDS("rangeFrame.Rda")

#Variables describuing how strict we truncate the spatio-temporal precision matrix
truncateC01 = 0.01
truncateC02 = 0.1


#Number of observations
n = dim(d)[1]

#Define a variable need for the spatio-temporal interaction
d$idGeneric0 = 1:n

#Define some variables needed for the spatial grid
cutoff = 20
maxEdge = c(100,200)


#Define the spatial grid
points = matrix(0,n,2)
points[,1] = d$UTMX
points[,2] = d$UTMY
boundary = INLA::inla.nonconvex.hull(points)
boundary2 = INLA::inla.nonconvex.hull(points,convex = -0.35)
mesh = INLA::inla.mesh.2d(loc = points, boundary = list(boundary,boundary2),max.edge = maxEdge,cutoff = cutoff)

#Create the spde-object
spde = INLA::inla.spde2.matern(mesh = mesh , alpha = 2)

#Define the A-matrix needed for inference with the spatial grid
A.est = INLA::inla.spde.make.A(mesh = mesh,
                               loc = points)

#Define field indices
field.indices = INLA::inla.spde.make.index("field", n.spde = spde$n.spde)



#Extract the range parameters in the spatio-temporal interaction
#rangeParameters = readRange()
thetaA = rangeParameters$rCodNB[1]
thetaT= rangeParameters$rCodNB[2]
#Define some variables needed for the spatio-temporal interaction
avstand = matrix(0,n,2)
avstand[,1] = d$UTMX; avstand[,2] = d$UTMY;
dis = as.matrix(sp::spDists(avstand,longlat=FALSE))
timeDis = as.matrix(stats::dist(d[,"daysSince1994"]))
S = exp(-dis/thetaA - timeDis/thetaT)
cholesky = chol(S)
C0 = chol2inv(cholesky)
C0[C0 <truncateC01 & C0>-truncateC01& S<truncateC02] = 0



#Define the formula used by R-INLA
hyperGen0 = list(prec = list(initial = log(0.58)))


formula=cod~1


estimate.inla.Olav <- function(formula, d=d, reff = "+ offset(log(dist))+f(idGeneric0,model=\"generic0\",Cmatrix = C0,hyper = hyperGen0)+
  f(field,model = spde, initial = c(3.2,-5.2))+
  f(codYear,model = 'iid',initial = log(3.8400))")
{
  
  X <- cbind(idGeneric0 = 1:n,codYear = d$codYear,dist = d$dist,as.data.frame(model.matrix(object = formula,d = d)[,-1]))
  nnames<-paste0(array("V",length(X)-3),1:(length(X)-3))
  if(length(names(X))>4)
    names(X)[4:length(X)]<-nnames
  else nnames <- ""
  #Create the stack which is the structure used for inference with R-INLA
  stack = INLA::inla.stack(data = list(cod = d$cod),
                           A = list(A.est,1),
                           effects =
                             list(c(field.indices,list(Intercept = 1)),
                                  X))
  fla<-as.character(formula)
  formulaFish = as.formula(paste0(fla[2],"~-1+",paste(nnames,collapse  = "+"),reff))
  #Run R-INLA with good starting values
  timeBefore=Sys.time();
  out <- INLA::inla(formula = formulaFish,
                      data = INLA::inla.stack.data(stack,spde = spde),
                      family = "zeroinflatednbinomial0",
                      control.predictor = list(A=INLA::inla.stack.A(stack)),
                      control.compute=list(config=TRUE),
                      control.inla=list(int.strategy="eb",interpolator = "gaussian",optimiser = "domin"),
                      control.family = list(hyper = list(
                        "size" = list(initial = 0.78),
                        prob = list(
                          initial = -20, #Probability for 0 is here set to 0.
                          fixed = TRUE))))
  timeUsed = round((Sys.time()-timeBefore),1)
  
  
  print(paste("Det tok tid: ", timeUsed,sep = ""))
  coef<-out$summary.fixed$mode
  coef[1]<-coef[1]+out$summary.hyperpar$mode[1]/(1-out$summary.hyperpar$mode[2])
  rm(stack)
  gc()
  return(list(mlik = out$mlik[1],waic = -out$mlik[1] , dic = -out$mlik[1], summary.fixed =list(mean =coef)))
}

estimate.inla.Olav(formula = formula,d = d)


formula1 = as.formula(paste(colnames(d)[8],"~ 1 +",paste0(colnames(d)[-c(1,8,11,13,28,37)],collapse = "+")))
print(formula1)

cosi<-function(x)cos(x/180*pi)
sini<-function(x)sin(x/180*pi)
expi<-function(x)
{
  r<-exp(x)
  if(r==Inf)
    return(10000000)
  else
    return(r)
}

InvX<-function(x)
{
  if(x==0)
    return(10000000)
  else
    return(1/x)
  
}
#print(estimate.gamma.cpen (formula=formula1, data=train))
data.example<-d

res = runemjmcmc(formula = formula1,data = data.example,estimator = estimate.inla.Olav,estimator.args =  list(d = data.example),recalc_margin = 249, save.beta = F,interact = T,outgraphs=T,relations = c("","cosi","sigmoid","tanh","atan","sini"),relations.prob =c(0.9,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=2,mutation_rate = 250,last.mutation=5000, max.tree.size = 5, Nvars.max =15,p.allow.replace=0.9,p.allow.tree=0.2,p.nor=0.9,p.and = 0.9),n.models = 10000,unique =T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
  max.N.glob=as.integer(10),
  min.N.glob=as.integer(5),
  max.N=as.integer(3),
  min.N=as.integer(1),
  printable = T))

mySearch$fparam
ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
ppp$p.post
length(ppp$m.post)






