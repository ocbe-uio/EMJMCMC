#inference
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")

#***********************IMPORTANT******************************************************
# if a multithreaded backend openBLAS for matrix multiplications
# is installed on your machine, please force it to use 1 thread explicitly
library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)
#***********************IMPORTANT******************************************************

X=read.csv("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/DBRM%20supplementaries/kepler%20and%20mass/exa1.csv")
data.example = as.data.frame(X)

#specify the initial formula
formula1 = as.formula(paste(colnames(X)[5],"~ 1 +",paste0(colnames(X)[-5],collapse = "+")))

#a set of nonlinearities that will be used in the DBRM model
sini=function(x)sin(x/180*pi)
expi=function(x)exp(-abs(x))
logi =function(x)log(abs(x)+1)
troot=function(x)abs(x)^(1/3)
to23=function(x)abs(x)^(2.3)
to35=function(x)abs(x)^(3.5)


#specify the estimator function returning p(Y|m)p(m), model selection criteria and the vector of the modes for the beta coefficients
estimate.gamma.cpen = function(formula, data,r = 1.0/223.0,logn=log(223.0),relat=c("to23","expi","logi","to35","sini","troot","sigmoid"))
{
  fparam=NULL
  fmla.proc=as.character(formula)[2:3]
  fobserved = fmla.proc[1]
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam =stri_split_fixed(str = fmla.proc[2],pattern = "+I",omit_empty = F)[[1]]
  sj=(stri_count_fixed(str = fparam, pattern = "*"))
  sj=sj+(stri_count_fixed(str = fparam, pattern = "+"))
  for(rel in relat)
    sj=sj+(stri_count_fixed(str = fparam, pattern = rel))
  sj=sj+1
  tryCatch(capture.output({
    out = glm(formula = formula,data = data, family = gaussian)
    mlik = (-(out$deviance -2*log(r)*sum(sj)))/2
    waic = -(out$deviance + 2*out$rank)
    dic =  -(out$deviance + logn*out$rank)
    summary.fixed =list(mean = coefficients(out))
    
  }, error = function(err) {
    print(err)
    mlik = -10000
    waic = -10000
    dic =  -10000
    summary.fixed =list(mean = array(0,dim=length(fparam)))
  }))
  return(list(mlik = mlik,waic = waic , dic = dic,summary.fixed =summary.fixed))
  
}
#define the number or cpus
M = 32
#define the size of the simulated samples
NM= 1000
#define \k_{max} + 1 from the paper
compmax = 32
#define treshold for preinclusion of the tree into the analysis
th=(10)^(-5)
#define a final treshold on the posterior marginal probability for reporting a tree
thf=0.05
#specify tuning parameters of the algorithm for exploring DBRM of interest
#notice that allow_offsprings=3 corresponds to the GMJMCMC runs and
#allow_offsprings=4 -to the RGMJMCMC runs
res1 = pinferunemjmcmc(n.cores = M, report.level =  0.2, num.mod.best = NM ,simplify = T,runemjmcmc.params = list(formula = formula1,data = data.example,estimator = estimate.gamma.cpen,estimator.args =  list(data = data.example),recalc_margin = 249, save.beta = F,interact = T,outgraphs=F,relations=c("to23","expi","logi","to35","sini","troot","sigmoid"),relations.prob =c(0.1,0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=3,mutation_rate = 250,last.mutation=10000, max.tree.size = 5, Nvars.max =15,p.allow.replace=0.9,p.allow.tree=0.01,p.nor=0.9,p.and = 0.9),n.models = 10000,unique =T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
  max.N.glob=as.integer(10),
  min.N.glob=as.integer(5),
  max.N=as.integer(3),
  min.N=as.integer(1),
  printable = F)))

print(res1$feat.stat)
