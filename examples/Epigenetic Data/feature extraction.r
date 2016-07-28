
library(parallel)
require(stats)
#install.packages("stringi")
library("stringi")

getPlusStrand<-function(chromseq,start.p,stop.p)
{
  return(stri_sub(str = chromseq,from = start.p+1, to = stop.p+1))
}
getPlus3bp<-function(start.p)
{
  return(stri_sub(str = CH1Seq,from = start.p+1, to = start.p+3))
}
getMinus3bp<-function(start.p)
{
  res<-character(length = 3)
  for(i in 1:3)
  {
    cur<-stri_sub(str = CH1Seq,from = start.p-2+i, length = 1)
    res[3-i+1]<-ifelse(test =cur=="G","C",ifelse(cur=="C","G",ifelse(cur=="A","T",ifelse(cur=="T","A",cur))))
  }
  return(stri_flatten(res,collapse = ""))
}
getMinusStrand<-function(chromseq,start.p,stop.p)
{

  l.ca<-stop.p-start.p+1
  res<-character(length = l.ca)
  for(i in 1:l.ca)
  {
    cur<-stri_sub(str = CH1Seq,from = start.p-2+i, length = 1)
    res[l.ca-i+1]<-ifelse(test = cur=="G","C",ifelse(cur=="C","G",ifelse(cur=="A","T",ifelse(cur=="T","A",cur))))
  }
  return(stri_flatten(res,collapse = ""))
}

#read and preprocess the data
workdir<-"/mn/sarpanitu/ansatte-u2/aliaksah/workspace/epigenetic features/"


X <- read.table(paste(workdir,"GSM1085222_mC_calls_Col_0.tsv",sep = "",collapse = ""), header=T)

#read and preprocess the annotation
GGroup <- read.table(paste(workdir,"GeneGroups.csv",sep = "",collapse = ""), header=T,sep = '\t',blank.lines.skip = TRUE,fill = TRUE)

YY <- read.table(paste(workdir,"Genes.csv",sep = "",collapse = ""), header=T,sep = '\t')


Express<- read.table(paste(workdir,"GSM1086840_Col_0_expression2.tsv",sep = "",collapse = ""), header=T,sep = '\t')

# read the complete genome file
fileName <- paste(workdir,"GCF_000001735.3_TAIR10_genomic.fna",sep = "",collapse = "")
CH1Seqbuf<-paste(readLines(fileName), sep="\n",collapse = "")
CH1Seqbuf<-stri_split_fixed(str = CH1Seqbuf,pattern = ">NC")
length(CH1Seqbuf[[1]])

CH1Seq<-CH1Seqbuf[[1]][[5]]
CH1Seq<-substr(CH1Seq,63,nchar(CH1Seq))
CH1Seq<-toupper(CH1Seq)


# play around with the original genetic sequence and its translation between the strands
# first plus strand
substr(CH1Seq,76, 78)# 75 ok
getPlusStrand(CH1Seq,75,77)
getPlus3bp(75)
substr(CH1Seq,77, 79)# 76 ok
getPlusStrand(CH1Seq,76,78)
substr(CH1Seq,84, 86)# 83 ok
getPlusStrand(CH1Seq,83,85)
substr(CH1Seq,85, 87)# 84 ok
getPlusStrand(CH1Seq,84,86)
substr(CH1Seq,107, 109)# 106 ok
getPlusStrand(CH1Seq,106,108)
# now minus strand
substr(CH1Seq,32, 34)# 75 ok
getMinusStrand(CH1Seq,33,35)
getMinus3bp(33)
getMinusStrand(CH1Seq,1564692,1564694)

lprom_length<-1500
laftr_length<-500
Express$tracking_id<-as.character(Express$tracking_id)
YY$Locus.tag<-as.character(YY$Locus.tag)
nchrom <- 5
for(chromi in 1:nchrom)
{
  CH1Seq<-CH1Seqbuf[[1]][[1+chromi]]
  CH1Seq<-substr(CH1Seq,63,nchar(CH1Seq))
  CH1Seq<-toupper(CH1Seq)

  # seems to work well
  # proceed with extracting all the unknown C bases
  c.ids.pls<-(gregexpr("C", CH1Seq))[[1]]-1
  c.ids.min<-(gregexpr("G", CH1Seq))[[1]]-1
  c.ids.pls.pres.id<-which(X$chrom==chromi & X$strand=="+")
  c.ids.min.pres.id<-which(X$chrom==chromi & X$strand=="-")
  c.ids.pls.pres<-X$pos[c.ids.pls.pres.id]
  c.ids.min.pres<-X$pos[c.ids.min.pres.id]
  c.ids.data.pls<-which(c.ids.pls %in% c.ids.pls.pres)
  c.ids.data.min<-which(c.ids.min %in% c.ids.min.pres)
  X.ch1<-data.frame(as.array(sort(c(c.ids.pls,c.ids.min))))
  names(X.ch1)<-"pos"
  X.ch1$strand<-"+"
  min.ids<-which(X.ch1$pos %in% c.ids.min)
  pls.ids<-which(X.ch1$pos %in% c.ids.pls)
  X.ch1$strand[min.ids]<-"-"
  row.names(X.ch1)<-(1):(length(c.ids.pls)+length(c.ids.min))
  X.ch1$mc_class<-NA

  system.time(
    X.ch1$mc_class[pls.ids]<-getPlus3bp(start.p = c.ids.pls)
  )
  #the next transformation might take a bit of time (ca 5-25 min)
  system.time(
    X.ch1$mc_class[min.ids]<-mclapply(X= c.ids.min, FUN = getMinus3bp)
  )
  pres.ids<-sort(c(c.ids.pls.pres.id,c.ids.min.pres.id))
  c.min.ids<-which(X.ch1$pos %in% c.ids.pls.pres)
  c.pls.ids<-which(X.ch1$pos %in% c.ids.min.pres)
  c.pres.ids<-sort(c(c.pls.ids,c.min.ids))
  length(c.pres.ids)
  length(pres.ids)
  X.ch1$methylated_bases<-NA
  X.ch1$methylated_bases[c.pres.ids]<-X$methylated_bases[pres.ids]
  X.ch1$total_bases<-NA
  X.ch1$total_bases[c.pres.ids]<-X$total_bases[pres.ids]
  X.ch1$unmethylated_bases<-NA
  X.ch1$unmethylated_bases[c.pres.ids]<-X$total_bases[pres.ids]-X$methylated_bases[pres.ids]
  cc<-diff(X.ch1$pos,lag = 1)
  cc<-c(0,cc)
  X.ch1$base_dist <- cc
  dim(X.ch1)
  rm(cc)
  X.ch1$CHG <- as.integer((substr(as.character(X.ch1$mc_class), 3, 3)=="G" & substr(as.character(X.ch1$mc_class), 2, 2)!="G"))
  X.ch1$CG <- as.integer(substr(as.character(X.ch1$mc_class), 2, 2)=="G")
  X.ch1$CHH <- as.integer(substr(as.character(X.ch1$mc_class), 2, 2)!="G" & substr(as.character(X.ch1$mc_class), 3, 3)!="G")

  X.ch1$DT1 <- as.integer(X.ch1$base_dist==1)
  X.ch1$DT2 <- as.integer(X.ch1$base_dist==2)
  X.ch1$DT3 <- as.integer(X.ch1$base_dist==3)
  X.ch1$DT4 <- as.integer(X.ch1$base_dist==4)
  X.ch1$DT5 <- as.integer(X.ch1$base_dist==5)
  X.ch1$DT6_20 <- as.integer(X.ch1$base_dist>=6 & X.ch1$base_dist<20)
  X.ch1$DT20_inf <- as.integer(X.ch1$base_dist>=20)
  X.ch1$chrom<-chromi
  # the following procedure was not vectorized and thus might take a bit of time
  Y<-YY[which(YY$Replicon.Name==chromi),]
  for(i in 1:nrow(GGroup))
  {
    #print(i)
    if(!(as.character(GGroup$TYPE[[i]]) %in% colnames(X.ch1)))
      X.ch1[[as.character(GGroup$TYPE[[i]])]]<-0
    set<-which(gsub(" ", "", as.character(Y$Locus.tag), fixed = TRUE) == gsub(" ", "", as.character(GGroup$GENE_ID[[i]]), fixed = TRUE))
    set<-c(set,((which(toupper(gsub(" ", "", as.character(Y$Locus.tag), fixed = TRUE)) == gsub(" ", "", as.character(GGroup$GENE_ID[[i]]), fixed = TRUE)))))
    set<-unique(set)
    if((  length(set)>1 ))
    {

      print(paste(as.character(GGroup$GENE_ID[[i]])," is found with protein products:"))
      print(as.character(Y$Protein.product[set]))
      for(j in set)
      {
        #print(Y$Start[j])
        set4<-which(as.integer(X.ch1$pos)>=as.integer(Y$Start[j]) & as.integer(X.ch1$pos)<=as.integer(Y$Stop[j]))
        if(length(set4)>0)
        {
          print(paste(as.character(Y$Protein.product[set]),"found in epi data with ",length(set4)," base pairs"))
          X.ch1[[as.character(GGroup$TYPE[[i]])]][set4]<-1
        }
      }

    }
    else
    {
      print(paste(as.character(GGroup$GENE_ID[[i]])," is not found in EPI dataset"))
    }
  }
  X.ch1$none <- 1
  X.ch1$none <-X.ch1$none - (X.ch1$Mα + X.ch1$Mβ + X.ch1$Mβ + X.ch1$Mγ + X.ch1$Mδ + X.ch1$MIKC)
  Express.ch<-Express[which(Express$chrom==chromi),]
  # now mark all promoters, coding regions and transcriptions
  X.ch1$coding<-0
  X.ch1$promoter<-0
  X.ch1$aftron<-0
  #define length of the promoter regions
  lchrom<-length(X.ch1$chrom)
  X.ch1$express<-NA
  Y<-YY[which(YY$Replicon.Name==chromi & YY$Strand=="-"),]
  expr.id<- -1
  for(i in 1:nrow(Y))
  {
    if(i%%10==0)
      print(paste(i,"genes proceed"))
    l.expr.id<-expr.id
    expr.id<-which(Express.ch$tracking_id == Y$Locus.tag[i])
    if(length(expr.id)==0)
    {
      expr.id<- -1
      next
    }
    if(l.expr.id==expr.id)
      next
    set4<-NULL
    set4<-which((X.ch1$pos)>=(Express.ch$start[expr.id]) &(X.ch1$pos)<=(Express.ch$end[expr.id]))
    if(length(set4)>0)
    {

      ls4<-set4[length(set4)]
      fs4<-set4[1]


      set7<-which(X.ch1$pos>=(X.ch1$pos[fs4]-laftr_length) & X.ch1$pos<X.ch1$pos[fs4])
      set6<-which(X.ch1$pos<=(X.ch1$pos[ls4]+lprom_length) & X.ch1$pos>X.ch1$pos[ls4])

      X.ch1$coding[set4]<-1
      X.ch1$promoter[set6]<-1
      X.ch1$aftron[set7]<-1
      X.ch1$promoter[set4]<-0
      X.ch1$aftron[set4]<-0
      if(!is.null(expr.id))
      {
        X.ch1$express[set4]<-Express.ch$fpkm[expr.id]
        #X.ch1$express[set6]<-Express.ch$fpkm[expr.id]
        #X.ch1$express[set7]<-Express.ch$fpkm[expr.id]
      }

    }
  }
  Y<-YY[which(YY$Replicon.Name==chromi & YY$Strand=="+"),]
  expr.id<- -1
  for(i in 1:nrow(Y))
  {
    if(i%%10==0)
      print(paste(i,"genes proceed"))
    l.expr.id<-expr.id
    expr.id<-which(Express.ch$tracking_id == Y$Locus.tag[i])
    if(length(expr.id)==0)
    {
      expr.id<- -1
      next
    }
    if(l.expr.id==expr.id)
      next
    set4<-NULL
    set4<-which((X.ch1$pos)>=(Express.ch$start[expr.id]) &(X.ch1$pos)<=(Express.ch$end[expr.id]))
    if(length(set4)>0)
    {

      ls4<-set4[length(set4)]
      fs4<-set4[1]

      set7<-which(X.ch1$pos<=(X.ch1$pos[ls4]+laftr_length) & X.ch1$pos>X.ch1$pos[ls4])
      set6<-which(X.ch1$pos>=(X.ch1$pos[fs4]-lprom_length) & X.ch1$pos<X.ch1$pos[fs4])

      X.ch1$coding[set4]<-1
      X.ch1$promoter[set6]<-1
      X.ch1$aftron[set7]<-1
      X.ch1$promoter[set4]<-0
      X.ch1$aftron[set4]<-0
      if(!is.null(expr.id))
      {
        X.ch1$express[set4]<-Express.ch$fpkm[expr.id]
        #X.ch1$express[set6]<-Express.ch$fpkm[expr.id]
        #X.ch1$express[set7]<-Express.ch$fpkm[expr.id]
      }
    }
  }


  length(which(X.ch1$promoter==1))
  length(which(X.ch1$aftron==1))

  ss<-which(X.ch1$promoter + X.ch1$aftron == 2)
  length(ss)
  #X.ch1$promoter[ss]<-0
  #X.ch1$aftron[ss]<-0
  ss1<-which(X.ch1$promoter==1)
  ss2<-which(X.ch1$aftron==1)
  ss3<-which(X.ch1$coding==1)
  length(ss1)
  length(ss2)
  length(ss3)

  classes<-which(X.ch1$coding +  X.ch1$aftron + X.ch1$promoter > 0 & X.ch1$strand == "+")
  length(classes)
  classes<-which(X.ch1$coding +  X.ch1$aftron + X.ch1$promoter > 0 & X.ch1$strand == "-")
  length(classes)

  set5<-which(X.ch1$coding==0 & X.ch1$promoter==0 & X.ch1$aftron==0)
  print(paste(length(set5)," positions are not classified"))
  X.ch1$coding[set5]<-NA
  X.ch1$promoter[set5]<-NA
  X.ch1$aftron[set5]<-NA
  X.ch1$mc_class<-as.character(X.ch1$mc_class)
  write.csv(x = X.ch1, file = paste("EpigenNew",chromi,".csv",collapse = "",sep = ""))

}


View(X.ch1[212000:217000,])
