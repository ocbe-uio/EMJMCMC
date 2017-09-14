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
