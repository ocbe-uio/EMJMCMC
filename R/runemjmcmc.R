runemjmcmc<-function(formula, data, secondary = vector(mode="character", length=0), latnames="",
                   estimator,estimator.args = "list",n.models,p.add.default = 1,p.add = 0.5, unique = F,save.beta=F, locstop.nd = F, latent="",max.cpu=4,max.cpu.glob=2,create.table=T, hash.length = 20, presearch=T, locstop =F ,pseudo.paral = F,interact = F,deep.method =1,relations = c("","sin","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.1,0.1,0.1,0.1,0.1,0.1),gen.prob = c(1,10,5,1,0),pool.cross = 0.9,p.epsilon = 0.0001, del.sigma = 0.5,pool.cor.prob = F, interact.param=list(allow_offsprings=2,mutation_rate = 100,last.mutation=2000, max.tree.size = 10000, Nvars.max = 100, p.allow.replace = 0.7,p.allow.tree=0.1,p.nor=0.3,p.and = 0.7), prand = 0.01,keep.origin = T, sup.large.n = 5000, recalc_margin = 2^10, create.hash=F,interact.order=1,burn.in=1, eps = 10^6, max.time = 120,max.it = 25000, print.freq = 100,outgraphs=F,advanced.param=NULL, distrib_of_neighbourhoods=t(array(data = c(7.6651604,16.773326,14.541629,12.839445,2.964227,13.048343,7.165434,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           0.9936905,15.942490,11.040131,3.200394,15.349051,5.466632,14.676458,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           1.5184551,9.285762,6.125034,3.627547,13.343413,2.923767,15.318774,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           14.5295380,1.521960,11.804457,5.070282,6.934380,10.578945,12.455602,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           6.0826035,2.453729,14.340435,14.863495,1.028312,12.685017,13.806295),dim = c(7,5))),  distrib_of_proposals = c(76.91870,71.25264,87.68184,60.55921,15812.39852))
{
# a function that creates an EMJMCMC2016 object with specified values of some parameters and default values of other parameters

#first create the object
assign("data.example",data, envir=globalenv())

variables <- simplify.formula(formula,names(data.example))
assign("fparam.example",variables$fparam, envir=globalenv())
assign("fobserved.example",variables$fobserved, envir=globalenv())

#for(i in 1:length(fparam.example))
#{
#  fparam.example[i]<<-paste("I(",variables$fparam[i],")",sep = "")
#}
fparam.tmp<- as.vector(sapply(FUN = paste,"I(",variables$fparam,")",sep="")[,1])
if(latnames[1]!="")
  fparam.example<<-c(fparam.tmp,latnames)
else
  fparam.example<<-fparam.tmp
#print(fparam.tmp)
#print(fparam.example)
assign("mySearch",methods::new(structure("EMJMCMC2016", package = "EMJMCMC")), envir=globalenv())
if(length(secondary)>0)
  mySearch$filtered <<- sapply(FUN = paste,"I(",secondary,")",sep="")
mySearch$estimator <<- estimator
mySearch$latnames <<- latnames
mySearch$estimator.args <<- estimator.args
mySearch$latent.formula <<- latent
mySearch$save.beta <<- save.beta
mySearch$prand<<-prand
mySearch$p.add<<-array(p.add,length(fparam.example))
mySearch$p.add.default<<-p.add.default
mySearch$recalc.margin <<- as.integer(recalc_margin)
mySearch$max.cpu <<- as.integer(max.cpu)
mySearch$locstop.nd <<- locstop.nd
mySearch$pool.cor.prob<<-pool.cor.prob
mySearch$sup.large.n<<-as.integer(sup.large.n)
mySearch$max.cpu.glob <<- as.integer(max.cpu.glob)
mySearch$deep.method <<- as.integer(deep.method)
if(interact)
{
  mySearch$allow_offsprings <<- as.integer(interact.param$allow_offsprings)
  mySearch$mutation_rate <<- as.integer(interact.param$mutation_rate)
  mySearch$Nvars.max <<- as.integer(interact.param$Nvars.max)
  mySearch$max.tree.size <<- as.integer(interact.param$max.tree.size)
  mySearch$p.allow.replace <<-  interact.param$p.allow.replace
  mySearch$p.allow.tree <<-  interact.param$p.allow.tree
  mySearch$p.epsilon <<-  p.epsilon
  mySearch$keep.origin<<- keep.origin
  mySearch$sigmas<<-relations
  mySearch$sigmas.prob<<-relations.prob
  mySearch$del.sigma<<-del.sigma
  mySearch$pool.cross<<-pool.cross
  mySearch$gen.prob<<-gen.prob
  mySearch$p.nor <<- interact.param$p.nor
  mySearch$p.and <<- interact.param$p.and
  mySearch$last.mutation <<- as.integer(interact.param$last.mutation)
}

if(!is.null(advanced.param))
{
  mySearch$max.N.glob<<-as.integer(advanced.param$max.N.glob)
  mySearch$min.N.glob<<-as.integer(advanced.param$min.N.glob)
  mySearch$max.N<<-as.integer(advanced.param$max.N)
  mySearch$min.N<<-as.integer(advanced.param$min.N)
  mySearch$printable.opt<<-advanced.param$printable
}

if(exists("hashStat"))
{
  hash::clear(hashStat)
  remove(hashStat,envir=globalenv())
}
if(exists("statistics1"))
{
  remove(statistics,envir=globalenv() )
  remove(statistics1,envir=globalenv())
}
if(exists("hash.keys1"))
{
  remove(hash.keys,envir=globalenv())
  remove(hash.keys1,envir=globalenv())
}

if(create.table)
{
  if(pseudo.paral)
    mySearch$parallelize <<- lapply
  #carry the search (training out)
  assign("statistics1",bigmemory::big.matrix(nrow = 2 ^min((length(fparam.example)),hash.length)+1, ncol =  16+length(fparam.example)*save.beta,init = NA, type = "double"), envir=globalenv())
  assign("statistics",bigmemory::describe(statistics1), envir=globalenv())

  mySearch$g.results[4,1]<<-0
  mySearch$g.results[4,2]<<-0
  mySearch$p.add <<- array(data = 0.5,dim = length(fparam.example))
  if((length(fparam.example))>20)
  {
    mySearch$hash.length<<-as.integer(hash.length)
    mySearch$double.hashing<<-T
    hash.keys1 <<- bigmemory::big.matrix(nrow = 2 ^(hash.length)+1, ncol = length(fparam.example),init = 0, type = "char")
    hash.keys <<- bigmemory::describe(hash.keys1)
  }

}else if(create.hash)
{

  assign("hashStat", hash::hash(), envir=globalenv())
  mySearch$parallelize <<- lapply
  mySearch$hash.length<<-as.integer(20)
  mySearch$double.hashing<<-F
}
# now as the object is created run the algorithm
initsol=stats::rbinom(n = length(fparam.example),size = 1,prob = 0.5)
if(unique)
  resm<-mySearch$modejumping_mcmc(list(varcur=initsol,locstop=locstop,presearch=presearch,statid=5, distrib_of_proposals =distrib_of_proposals,distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = eps, trit = n.models*100, trest = n.models, burnin = burn.in, max.time = max.time, maxit = max.it, print.freq = print.freq))
else
  resm<-mySearch$modejumping_mcmc(list(varcur=initsol,locstop=locstop,presearch=presearch,statid=5, distrib_of_proposals =distrib_of_proposals,distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = eps, trit =  n.models, trest = n.models*100, burnin = burn.in, max.time = max.time, maxit = max.it, print.freq = print.freq))
ppp<-1
print("MJMCMC is completed")
if(create.table)
{
  print("Post Proceed Results")
  ppp<-mySearch$post_proceed_results(statistics1 = statistics1)
  truth = ppp$p.post # make sure it is equal to Truth column from the article
  truth.m = ppp$m.post
  truth.prob = ppp$s.mass
  ordering = sort(ppp$p.post,index.return=T)
  print("pi truth")
  sprintf("%.10f",truth[ordering$ix])
  sprintf(fparam.example[ordering$ix])
}
else if(create.hash)
{
  print("Post Proceed Results")
  ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
  truth = ppp$p.post # make sure it is equal to Truth column from the article
  truth.m = ppp$m.post
  truth.prob = ppp$s.mass
  ordering = sort(ppp$p.post,index.return=T)
  print("pi truth")
  sprintf("%.10f",truth[ordering$ix])
  sprintf(fparam.example[ordering$ix])
}

if(outgraphs)
{
  graphics::par(mar = c(10,4,4,2) + 4.1)
  graphics::barplot(resm$bayes.results$p.post,density = 46,border="black",main = "Marginal Inclusion (RM)",ylab="Probability",names.arg = mySearch$fparam,las=2)
  graphics::barplot(resm$p.post,density = 46,border="black",main = "Marginal Inclusion (MC)",ylab="Probability",names.arg = mySearch$fparam,las=2)
}

return(ppp)
}
