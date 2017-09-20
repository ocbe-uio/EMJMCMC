snps<-unique(dataNeigbourhoodS2$causSNPid)


#absolutely best
bests<-NULL
for(snp in snps)
{
  ids<-which(dataNeigbourhoodS2$causSNPid==snp)
  best<-dataNeigbourhoodS2$SNPid[ids[length(ids)-1]]
  bests<-c(bests,best)
}

formula.best<- as.formula(paste0("Y~1+",paste(bests,collapse = "+")))

estimate.lm.MBIC2(formula.best,data.example)


for(i in 1:1000)
  print(estimate.lm.MBIC2(as.formula(paste0("Y~1+",paste(bests[sample.int(20,sample.int(20,1))],collapse = "+"))),data.example)$mlik)
  
# best for some subset corresponding to the runs detected<-results[[3]]$fparam with 90% of power
# and best model found -7977.09000 

snps<-unique(dataNeigbourhoodS2$causSNPid)

snpf<-results[[3]]$fparam

snpf<-stri_replace(str = snpf,fixed = "I(",replacement = "")
snpf<-stri_replace(str = snpf,fixed = ")",replacement = "")

#absolutely best for run 3 in our case 
bests<-NULL
for(snp in snps)
{
  ids<-which(dataNeigbourhoodS2$causSNPid==snp)
  
  
  
  best<-dataNeigbourhoodS2$SNPid[max(which(dataNeigbourhoodS2$SNPid[ids]%in%snpf))]
  bests<-c(bests,best)
}

bests<-bests[-which(is.na(bests))]

formula.best<- as.formula(paste0("Y~1+",paste(bests,collapse = "+")))

estimate.lm.MBIC2(formula.best,data.example)


for(i in 1:100)
  print(estimate.lm.MBIC2(as.formula(paste0("Y~1+",paste(bests[sample.int(18,sample.int(18,1))],collapse = "+"))),data.example)$mlik)

formula.best.mosgwa1<- (Y~1+
+rs4916071
+rs7515154
+rs594249
+rs1321108
+rs13294
+rs2841981
+rs12404596
+rs7522326
+rs4520444
+rs823130
+rs1342590
+rs1045287
+rs10802802
+rs10788764)

estimate.lm.MBIC2(formula.best.mosgwa1,data.example)

library(stringi)
j=0
fdr.tot=0
pow.tot=0
for(i in 1:22)
  {
      res1<-read.csv(paste0("postGMJSIM_",i,".csv"),header = T,stringsAsFactors = F)
      if(!is.null(res1$tree))
        {
            print(paste0("Converged Iteration ",i))
            detected<-res1$tree
            detected<-stri_replace(str = detected,fixed = "I(",replacement = "")
            detected<-stri_replace(str = detected,fixed = ")",replacement = "")
            
              detect.true.unique<-unique(dataNeigbourhoodS2$causSNPid[which(dataNeigbourhoodS2$SNPid %in% detected)])
              detect.true<-which(detected %in% dataNeigbourhoodS2$SNPid)
              
                detlen<-length(detect.true.unique)
                totlen<-length(detected)-length(detect.true)+length(detect.true.unique)
                
                  pow=detlen/20
                  print(pow)
                  fdr=(totlen-detlen)/totlen
                  print(fdr)
                  j=j+1
                  fdr.tot = fdr.tot + fdr
                  pow.tot = pow.tot + pow
                  
                  
                  }else
                      {
                          if(!is.null(res$tree))
                            {
                                print(paste0("Diverged Iteration ",i))
                                detected<-res1$X
                                detected<-stri_replace(str = detected,fixed = "I(",replacement = "")
                                detected<-stri_replace(str = detected,fixed = ")",replacement = "")
                                
                                  detect.true.unique<-unique(dataNeigbourhoodS2$causSNPid[which(dataNeigbourhoodS2$SNPid %in% detected)])
                                  detect.true<-which(detected %in% dataNeigbourhoodS2$SNPid)
                                  
                                    detlen<-length(detect.true.unique)
                                    totlen<-length(detected)-length(detect.true)+length(detect.true.unique)
                                    
                                      print(detlen/20)
                                    print((totlen-detlen)/totlen)
                                    
                                    }
                        }
        
        }
pow.tot/j
fdr.tot/j




estimate.lm.MBIC2 <- function(formula, data, n = 5402, m = 24602, c = 16,u=1700)
{
  size<-stri_count_fixed(str = as.character(formula)[3],pattern = "+")
  
  if(size>u)
  {
    return(list(mlik = (-50000 + rnorm(1,0,1) - size*log(m*m*n/c) + 2*log(factorial(size+1))),waic = 50000 + rnorm(1,0,1), dic =  50000+ rnorm(1,0,1),summary.fixed =list(mean = array(0,size+1))))
  }else{
    out <- lm(formula = formula,data = data)
    logmarglik <- (2*logLik(out) - out$rank*log(m*m*n/c) + 2*log(factorial(out$rank)))/2
    # use dic and aic as bic and aic correspondinly
    return(list(mlik = logmarglik,waic = AIC(out) , dic =  BIC(out),summary.fixed =list(mean = coef(out))))
  }
}

system.time({
  
estimate.lm.MBIC2(formula = formula1,data = data.example)

})


estimate.fastlm.MBIC2 <- function(formula, data, n = 5402, m = 24602, c = 16,u=1700)
{
  size<-stri_count_fixed(str = as.character(formula)[3],pattern = "+")
  
  if(size>u)
  {
    return(list(mlik = (-50000 + rnorm(1,0,1) - size*log(m*m*n/c) + 2*log(factorial(size+1))),waic = 50000 + rnorm(1,0,1), dic =  50000+ rnorm(1,0,1),summary.fixed =list(mean = array(0,size+1))))
  }else{
    out <- biglm(formula = formula,data = data)
    logmarglik <- (2*(-deviance(out)) - size*log(m*m*n/c) + 2*log(factorial(size)))/2
    # use dic and aic as bic and aic correspondinly
    return(list(mlik = logmarglik,waic = AIC(out) , dic =  AIC(out),summary.fixed =list(mean = coef(out))))
  }
}

system.time({
estimate.fastlm.MBIC2(formula = formula1,data = data.example)
})