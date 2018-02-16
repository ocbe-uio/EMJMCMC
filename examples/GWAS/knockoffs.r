#install.packages("knockoff")
#install.packages("SNPknock")
#install.packages("doMC")


library(knockoff)
library(SNPknock)
library(plyr)
library(dplyr)
library(magrittr)
library(parallel)
library(stringi)
library(doMC)

setwd("/mn/sarpanitu/ansatte-u2/aliaksah/abeldata/gwassim/")



options("exppression" = 200000)
##########################################
#: Calculate True h^2 for each Scenario :#
##########################################
simPars <- 
  read.table("scripts/simulationParameters.txt", 
             header = TRUE, 
             stringsAsFactors = FALSE) %>%
  select(SNPId = SNP_causal, Pos = BP_causal, MAF = MAF_causal, S1, S2, S3, S4)

causSNPsS1 <- with(simPars, SNPId[S1 != 0])
causSNPsS2 <- with(simPars, SNPId[S2 != 0]) 
causSNPsS3 <- with(simPars, SNPId[S3 != 0]) 
causSNPsS4 <- with(simPars, SNPId[S4 != 0]) 

genoData <- read.table("data-geno/CHR1_NFBC.raw",
                       header = TRUE,
                       stringsAsFactors = FALSE)
names(genoData) <- gsub("_.$", "", names(genoData)) %>% gsub("\\.", "-", .)

genoData <- genoData[c("FID", "IID", simPars$SNPId)]

expectedData <- data.frame(as.matrix(genoData[-(1:2)]) %*% as.matrix(simPars[c('S1', 'S2', 'S3', 'S4')])) 

explainedVar <- apply(expectedData, 2, function(x) var(x)/(var(x)+1)) # True h^2


############################################
#: Calculate Neigbourhoods of causal SNPs :#
#: Definition: Neigbourhood 		  :#
#	       corr > .3 and +/- 1(2) MBP :#
############################################
distance_cutoff <- 1e6
corr_cutoff <- .3

physMap <- read.table("data-geno/CHR1_NFBC.bim",
                      header = FALSE,
                      stringsAsFactors = FALSE)[c(2, 4)] %>%
  select(SNPid = V2, pos = V4) %>%
  filter(!grepl('^cnv', SNPid))	   

genoData <- read.table("data-geno/CHR1_NFBC.raw", 
                       header = TRUE, 
                       stringsAsFactors = FALSE)
names(genoData) <- gsub("_.$", "", names(genoData)) %>% gsub("\\.", "-", .)	 

findNeigbour <- function(causSNP) {
  #causSNP <- "rs1498308"
  #print(causSNP)
  physMap %>%
    mutate(distance = abs(pos - pos[SNPid == causSNP])) %>%
    filter(distance <= distance_cutoff) %>%
    mutate(corr = sapply(SNPid, function(x) abs(cor(genoData[causSNP], genoData[x]))),
           causSNPid = causSNP) %>%
    filter(corr >= corr_cutoff)
}

#dataNeigbourhoodS1 <- plyr::ldply(causSNPsS1, findNeigbour)
dataNeigbourhoodS2 <- plyr::ldply(causSNPsS2, findNeigbour)
#dataNeigbourhoodS3 <- plyr::ldply(causSNPsS3, findNeigbour)
#dataNeigbourhoodS4 <- plyr::ldply(causSNPsS4, findNeigbour)

rm(genoData)
rm(expectedData)
rm(physMap)
rm(simPars)
rm(causSNPsS4)
rm(causSNPsS3)
rm(causSNPsS1)
rm(findNeigbour)
gc()


pheno<-read.csv(paste0("data_S2_nocausal_5402/pimass/data.recode.pheno_",1),header = F)
geno<-t(read.csv("data_S2_nocausal_5402/pimass/data.recode.mean.geno.txt",header = F,stringsAsFactors = F))
names<-geno[1,]
geno<-as.data.frame(geno[-c(1,2,3),])
names(geno)<-names
geno$Y<-pheno$V1
geno <- as.data.frame(mclapply(geno, as.numeric))


MM = 10
M = 31
size.init=1000
NM= 1000

compmax = 101
th<-(10)^(-5)
thf<-0.001

data.example<-geno
rm(geno)
gc()


j=1 
fdr.tot=0
pow.tot=0
fps=0
mis=0

for(j in 1:100)
{

    
    set.seed(j)
    
    pheno<-read.csv(paste0("data_S2_nocausal_5402/pimass/data.recode.pheno_",j),header = F)
    data.example$Y<-as.numeric(pheno$V1)
    rm(pheno)
    cors<-cor(data.example$Y,data.example[,1:(dim(data.example)[2]-1)])
    gc()
    cov.names<-names[which(abs(cors)>0.025)]
    detected<-cov.names
    detected<-stri_replace(str = detected,fixed = "I(",replacement = "")
    detected<-stri_replace(str = detected,fixed = ")",replacement = "")
    detect.true.unique<-unique(dataNeigbourhoodS2$causSNPid[which(dataNeigbourhoodS2$SNPid %in% detected)])
    detect.true<-which(detected %in% dataNeigbourhoodS2$SNPid)
    detlen<-length(detect.true.unique)
    totlen<-length(detected)-length(detect.true)+length(detect.true.unique)
    print(detlen/20)
    print((totlen-detlen)/totlen)
    
    #create genotypes matrix
    X = data.example[,which(abs(cors)>0.025)]
    #create knockoffs matrix
    Xk = X
    W = stat.glmnet_coefdiff(X, Xk, data.example$Y, family="gaussian", cores = 32)
    plot(W, pch=16, cex=1)
    
   
    
    t = knockoff.threshold(W, fdr=0.1, offset=0)
    detected = which(W >= t)
    
    colors = rep("gray",length(W))
    colors[detected] = "blue"
    plot(W, col=colors, pch=16, cex=1); abline(h=t, lty=2)
    
    detected = cov.names[detected]
    print(detected)
  
    
    detect.true.unique<-unique(dataNeigbourhoodS2$causSNPid[which(dataNeigbourhoodS2$SNPid %in% detected)])
    detect.true<-which(detected %in% dataNeigbourhoodS2$SNPid)
    
    detlen<-length(detect.true.unique)
    totlen<-length(detected)-length(detect.true)+length(detect.true.unique)
    
    pow=detlen/20
    print(pow)
    fdr=(totlen-detlen)/totlen
    print(fdr)
    fps=fps+totlen-detlen
    mis=mis+totlen-detlen + 20 - detlen
    fdr.tot = fdr.tot + fdr
    pow.tot = pow.tot + pow
    
    gc()
    

}
