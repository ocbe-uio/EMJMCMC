#***********************IMPORTANT******************************************************
# if a multithreaded backend openBLAS for matrix multiplications
# is installed on your machine, please force it to use 1 thread explicitly
# library(RhpcBLASctl)
# blas_set_num_threads(1)
# omp_set_num_threads(1)
#***********************IMPORTANT******************************************************


data.example = read.csv("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/BGNLM/abalone%20age/abalone.data",header = F)
data.example$MS=as.integer(data.example$V1=="M")
data.example$FS=as.integer(data.example$V1=="F")
data.example$V1=data.example$V9
data.example$V9 = NULL

names(data.example) = c("Age","Length", "Diameter","Height","WholeWeight","ShuckedWeight","VisceraWeight","ShellWeight","Male","Femele")

set.seed(040590)
teid =  read.csv("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/BGNLM/abalone%20age/teid.csv",sep = ";")[,2]



test = data.example[teid,]
data.example = data.example[-teid,]

test_that("Input dataset is still roughly the same", {
  expect_equal(sum(test$Age), 9801)
})

#specify the initial formula
formula1 = as.formula(paste(colnames(test)[1],"~ 1 +",paste0(colnames(test)[-1],collapse = "+")))

#define the number or cpus
M = min(20, parallel::detectCores() - 1)
#define the size of the simulated samples
NM= 100
#define \k_{max} + 1 from the paper
compmax = 16
#define treshold for preinclusion of the tree into the analysis
th=(10)^(-5)
#define a final treshold on the posterior marginal probability for reporting a tree
thf=0.05
#specify tuning parameters of the algorithm for exploring DBRM of interest
#notice that allow_offsprings=3 corresponds to the GMJMCMC runs and
#
g = function(x) x
results=array(0,dim = c(2,100,5))

res1_seq <- suppressMessages(runemjmcmc(
  formula = formula1, data = data.example, estimator = estimate.gamma.cpen,
  estimator.args = list(data = data.example), recalc_margin = 249,
  save.beta = TRUE, interact = TRUE, outgraphs = FALSE,
  relations = c("to25", "expi", "logi", "to35", "troot", "sigmoid"),
  relations.prob = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
  interact.param = list(
    allow_offsprings = 3, mutation_rate = 250, last.mutation = 100,
    max.tree.size = 5, Nvars.max =15, p.allow.replace = 0.9,
    p.allow.tree = 0.01, p.nor = 0.9, p.and = 0.9
  ), n.models = 1000, unique = TRUE, max.cpu = 4, max.cpu.glob = 4,
  create.table = FALSE, create.hash = TRUE, pseudo.paral = TRUE,
  burn.in = 100, print.freq = 0L,
  advanced.param = list(
    max.N.glob = 10L,
    min.N.glob = 5L,
    max.N = 3L,
    min.N = 1L,
    printable = FALSE
  )
))

test_that("runemjmcmc outputs with correct elements", {
  expect_named(res1_seq, c("p.post", "m.post", "s.mass"))
  expect_length(res1_seq[["p.post"]], 15L)
  expect_gt(length(res1_seq[["m.post"]]), 1000L)
  expect_length(res1_seq[["s.mass"]], 1L)
  expect_lte(mean(res1_seq[["p.post"]]), 1)
  expect_lte(mean(res1_seq[["m.post"]]), 1)
  expect_equal(res1_seq[["s.mass"]], 0)
})

set.seed(2915224)

res1_par <- suppressMessages(pinferunemjmcmc(
  n.cores = M, report.level =  0.2, num.mod.best = NM, simplify = TRUE,
  predict = TRUE, test.data = as.data.frame(test), link.function = g,
  runemjmcmc.params = list(
    formula = formula1, data = data.example, estimator = estimate.gamma.cpen,
    estimator.args = list(data = data.example), recalc_margin = 249,
    save.beta = TRUE, interact = TRUE, outgraphs = FALSE,
    relations = c("to25", "expi", "logi", "to35", "troot", "sigmoid"),
    relations.prob = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
    interact.param = list(
      allow_offsprings = 3, mutation_rate = 250, last.mutation = 100,
      max.tree.size = 5, Nvars.max =15, p.allow.replace = 0.9,
      p.allow.tree = 0.01, p.nor = 0.9, p.and = 0.9
    ), n.models = 500, unique = TRUE, max.cpu = 1L, max.cpu.glob = 1L,
    create.table = FALSE, create.hash = TRUE, pseudo.paral = TRUE,
    burn.in = 100, print.freq = 0L,
    advanced.param = list(
      max.N.glob = 10L,
      min.N.glob = 5L,
      max.N = 3L,
      min.N = 1L,
      printable = FALSE
    )
  )
))

test_that("pinferunemjmcmc outputs with correct elements", {
  expect_named(
    res1_par,
    c("feat.stat", "predictions", "allposteriors", "threads.stats")
  )
  expect_gte(length(res1_par[["feat.stat"]]), 20L)
  expect_lte(length(res1_par[["feat.stat"]]), 30L)
  expect_equal(length(res1_par[["predictions"]]), 1000L)
  expect_equal(length(res1_par[["allposteriors"]]), 2L)
  expect_gte(length(res1_par[["threads.stats"]]), 1L)
  expect_lte(length(res1_par[["threads.stats"]]), 5L)
  expect_equal(mean(res1_par[["predictions"]]), 9.9, tolerance = 1e-1)
  expect_equal(
    res1_par[["threads.stats"]][[1]][["cterm"]], -6573, tolerance = 1e-1
  )
  expect_equal(
    res1_par[["threads.stats"]][[2]][["preds"]][1], 7.69, tolerance = 1e-1
  )
  expect_equal(
    res1_par[["threads.stats"]][[3]][["p.post"]][1], .977, tolerance = 1e-2
  )
  expect_equal(
    res1_par[["threads.stats"]][[4]][["post.populi"]], 0.0282, tolerance = 1e-3
  )
  expect_equal(
    res1_par[["threads.stats"]][[5]][["mliks"]][1], -7762, tolerance = 1e-1
  )
})

J <- seq_len(1L) # 1L to save time, but results are basically the same for 10L
for(j in J) {
  #specify the initial formula
  set.seed(j)

  res1 <- suppressMessages(
    pinferunemjmcmc(
      n.cores = M, report.level =  0.2, num.mod.best = NM, simplify = TRUE,
      predict = TRUE, test.data = as.data.frame(test), link.function = g,
      runemjmcmc.params = list(
        formula = formula1, data = data.example, estimator = estimate.gamma.cpen,
        estimator.args = list(data = data.example), recalc_margin = 249,
        save.beta = TRUE, interact = TRUE, outgraphs = FALSE,
        relations = c("to25", "expi", "logi", "to35", "troot", "sigmoid"),
        relations.prob = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
        interact.param = list(
          allow_offsprings = 3, mutation_rate = 250, last.mutation = 100,
          max.tree.size = 5, Nvars.max =15, p.allow.replace = 0.9,
          p.allow.tree = 0.01, p.nor = 0.9, p.and = 0.9
        ), n.models = 1000, unique = TRUE, max.cpu = 4, max.cpu.glob = 4,
        create.table = FALSE, create.hash = TRUE, pseudo.paral = TRUE,
        burn.in = 100, print.freq = 0L,
        advanced.param = list(
          max.N.glob = 10L,
          min.N.glob = 5L,
          max.N = 3L,
          min.N = 1L,
          printable = FALSE
        )
      )
    )
  )

  results[1,j,1]=  sqrt(mean((res1$threads.stats[[1]]$preds - test$Age)^2))
  results[1,j,2]=   sqrt(mean(abs(res1$threads.stats[[1]]$preds - test$Age)))
  results[1,j,3] =   cor(res1$threads.stats[[1]]$preds,test$Age)


  results[2,j,1]=  sqrt(mean((res1$predictions - test$Age)^2))
  results[2,j,2]=   sqrt(mean(abs(res1$predictions - test$Age)))
  results[2,j,3] =   cor(res1$predictions,test$Age)

  posteriorshell_1.csv <- res1$feat.stat # replaced write.csv()

  #print the run's metrics and clean the results
  resultsrun_1.csv <- results[, j, ] # replaces write.csv()

  sqrt_mean_preds <- sqrt(mean((res1$predictions - test$Age)^2))

  test_that(paste("Iteration", j, "results look ok"), {
    expect_named(
      res1,
      c("feat.stat", "predictions", "allposteriors", "threads.stats")
    )
    expect_equal(mean(res1[["predictions"]]), 9.9, tolerance = 1e-1)
    expect_equal(
      res1[["threads.stats"]][[1]][["cterm"]], -6573, tolerance = 1e-1
    )
    expect_equal(
      res1[["threads.stats"]][[2]][["preds"]][1], 6.48, tolerance = 1e-1
    )
    expect_equal(
      res1[["threads.stats"]][[3]][["p.post"]][1], .227, tolerance = 1e-2
    )
    expect_equal(
      res1[["threads.stats"]][[4]][["post.populi"]], 4.42e-26, tolerance = 1e-25
    )
    expect_equal(
      res1[["threads.stats"]][[5]][["mliks"]][1], -7761.802, tolerance = 1e-2
    )
    expect_equal(dim(results), c(2L, 100L, 5L))
    expect_equal(
      results[1, 1:10, 1:3][1, ],
      c(2.0479307, 1.2309200, 0.7598387),
      tolerance = 1e-6
    )
    expect_equal(
      results[2, 1:10, 1:3][1, ],
      c(2.0137489, 1.2230869, 0.7679209),
      tolerance = 1e-6
    )
    expect_equal(dim(posteriorshell_1.csv), c(10L, 2L))
    expect_equal(dim(resultsrun_1.csv), c(2L, 5L))
    expect_equal(mean(resultsrun_1.csv), .80, tolerance = 1e-1)
    expect_equal(sqrt_mean_preds, 2.01, tolerance = 1e-1)
  })
}

# for(j in 1:100)
# {

#   tmp = read.csv(paste0("resultsrun_",j,".csv"))

#   results[1,j,1]=  tmp[1,2]
#   results[1,j,2]=  tmp[1,3]
#   results[1,j,3] =   tmp[1,4]

#   if(tmp[2,2] == 0)
#     print(j)

#   results[2,j,1]=  tmp[2,2]
#   results[2,j,2]=   tmp[2,3]
#   results[2,j,3] =   tmp[2,4]

# }

# #make the joint summary of the runs, including min, max and medians of the performance metrics
# summary.results=array(data = NA,dim = c(2,15))

# for(i in 1:2){
# for(j in 1:5)
# {
#   summary.results[i,(j-1)*3+1]=min(results[i,,j])
#   summary.results[i,(j-1)*3+2]=median(results[i,,j])
#   summary.results[i,(j-1)*3+3]=max(results[i,,j])
# }
# }
# summary.results=as.data.frame(summary.results)


# #featgmj = hash()

# simplifyposteriors<-function(X,post,th=0.0001,thf=0.1,y = "Age")
# {
#   posteriors = (cbind(as.character(post[,1]),as.numeric(as.character(post[,2]))))

#   rhash<-hash()
#   for(i in 1:length(posteriors[,1]))
#   {


#     expr<-posteriors[i,1]
#     print(expr)
#     res<-model.matrix(data=X,object = as.formula(paste0(y,"~",expr)))
#     ress<-c(stri_flatten(round(sum(res[,2]),digits = 4),collapse = ""),stri_flatten(res[,2],collapse = ""),posteriors[i,2],expr)
#     if(!((ress[1] %in% keys(rhash))))
#       rhash[[ress[1]]]<-ress
#     else
#     {
#       if(ress[1] %in% keys(rhash))
#       {new
#         rhash[[ress[1]]][3]<- (as.numeric(rhash[[ress[1]]][3]) + as.numeric(ress[3]))
#         if(stri_length(rhash[[ress[1]]][4])>stri_length(expr))
#           rhash[[ress[1]]][4]<-expr
#       }
#     }

#   }
#   res<-as.data.frame(t(values(rhash)[c(3,4),]))
#   res$V1<-as.numeric(as.character(res$V1))
#   res<-res[which(res$V1>thf),]
#   res<-res[order(res$V1, decreasing = T),]
#   clear(rhash)
#   rm(rhash)
#   #res[which(res[,1]>1),1]<-1
#   colnames(res)<-c("posterior","tree")

#   row.names(res) = 1:length(res$posterior)
#   return(res)
# }




# featgmj = hash()

# for(j in 1:100)
# {
#   tmpf = read.csv(paste0("posteriorshell_",j,".csv"))
#   #tmp = simplifyposteriors(X = data.example,post =tmpf,y = "Age")
#   for(feat in as.character(tmpf$V1))
#   {
#     if(!has.key(hash = featgmj,key =  feat ))
#     {
#       featgmj[[feat]] = as.numeric(1)
#     } else{

#       featgmj[[feat]] =as.numeric(featgmj[[feat]]) + 1
#     }
#   }
# }

# tmp = simplifyposteriors(X = data.example,post =as.data.frame(cbind(keys(featgmj),as.numeric(values(featgmj)))),y = "Age")

# write.csv(x =tmp,row.names = T,file = "abalonefeat.csv")
# #print(paste0("end simulation ",j))
