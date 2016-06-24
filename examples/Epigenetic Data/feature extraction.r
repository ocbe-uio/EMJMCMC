#read and preprocess the data
workdir<-"/mn/sarpanitu/ansatte-u2/aliaksah/workspace/epigenetic features/"


X <- read.table(paste(workdir,"GSM1085222_mC_calls_Col_0.tsv",sep = "",collapse = ""), header=T)

#read and preprocess the annotation
GGroup <- read.table(paste(workdir,"GeneGroups.csv",sep = "",collapse = ""), header=T,sep = '\t',blank.lines.skip = TRUE,fill = TRUE)

Y <- read.table(paste(workdir,"Genes.csv",sep = "",collapse = ""), header=T,sep = '\t')


Express<- read.table(paste(workdir,"GSM1086840_Col_0_expression2.tsv",sep = "",collapse = ""), header=T,sep = '\t')

fileName <- paste(workdir,"GCF_000001735.3_TAIR10_genomic.fna",sep = "",collapse = "")
getPlusStrand<-function(chromseq,start.p,stop.p)
{
  return(substr(chromseq,start.p+1, stop.p+1))
}
getPlus3bp<-function(start.p)
{
  return(substr(CH1Seq,start.p+1, start.p+3))
}
getMinus3bp<-function(start.p)
{
  ca<-strsplit(x = substr(CH1Seq,start.p-1, start.p+1),split = "")[[1]]
  l.ca<-length(ca)
  res<-character(length = l.ca)
  for(i in 1:l.ca)
  {
    res[l.ca-i+1]<-ifelse(test = ca[i]=="G","C",ifelse(ca[i]=="C","G",ifelse(ca[i]=="A","T","A")))
  }
  return(paste(res,collapse = ""))
}
getMinusStrand<-function(chromseq,start.p,stop.p)
{
  ca<-strsplit(x = substr(chromseq,start.p-1, stop.p-1),split = "")[[1]]
  l.ca<-length(ca)
  res<-character(length = l.ca)
  for(i in 1:l.ca)
  {
    res[l.ca-i+1]<-ifelse(test = ca[i]=="G","C",ifelse(ca[i]=="C","G",ifelse(ca[i]=="A","T","A")))
  }
  return(paste(res,collapse = ""))
}
CH1Seq<-paste(readLines(fileName), sep="\n",collapse = "")
CH1Seq<-substr(CH1Seq,66,nchar(CH1Seq))
CH1Seq<-toupper(CH1Seq)
# play around with the original genetic sequence
# first plus strand
substr(CH1Seq,76, 78)# 75 ok
getPlusStrand(CH1Seq,75,77)
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


getMinusStrand(CH1Seq,1564692,1564694)
# seems to work well
# proceed with extracting all the unknown C bases
c.ids.pls<-(gregexpr("C", CH1Seq))[[1]]-1
c.ids.min<-(gregexpr("G", CH1Seq))[[1]]-1
c.ids.pls.pres.id<-which(X$chrom==1 & X$strand=="+")
c.ids.min.pres.id<-which(X$chrom==1 & X$strand=="-")
c.ids.pls.pres<-X$pos[c.ids.pls.pres.id]
c.ids.min.pres<-X$pos[c.ids.min.pres.id]
c.ids.data.pls<-which(c.ids.pls %in% c.ids.pls.pres)
c.ids.data.min<-which(c.ids.min %in% c.ids.min.pres)
X.ch1<-data.frame(as.array(c(c.ids.pls,c.ids.min)))
names(X.ch1)<-"pos"
X.ch1$strand<-"+"
X.ch1$strand[(length(c.ids.pls)+1):(length(c.ids.pls)+length(c.ids.min))]<-"-"
row.names(X.ch1)<-(1):(length(c.ids.pls)+length(c.ids.min))
X.ch1$mc_class<-NA
X.ch1$mc_class[(length(c.ids.pls)+1):(length(c.ids.pls)+length(c.ids.min))]<-getMinusStrand(chromseq =  CH1Seq,start.p = c.ids.min,c.ids.min+2)
X.ch1$mc_class[(1):(length(c.ids.pls))]<-mclapply(X=c.ids.pls,FUN = getPlus3bp)

getPlusStrand(chromseq =  CH1Seq,start.p =c(1,3,5),stop.p =c(3,5,7))

dim(X)
cc<-diff(X$pos,lag = 1)
cc<-c(0,cc)
X$base_dist <- cc
X$base_dist[which(X$base_dist<0)]<-0
X$CHG <- as.integer((substr(as.character(X$mc_class), 3, 3)=="G" & substr(as.character(X$mc_class), 2, 2)!="G"))
X$CG <- as.integer(substr(as.character(X$mc_class), 2, 2)=="G")
X$CHH <- as.integer(substr(as.character(X$mc_class), 2, 2)!="G" & substr(as.character(X$mc_class), 3, 3)!="G")
X$DT1 <- as.integer(X$base_dist==1)
X$DT2 <- as.integer(X$base_dist==2)
X$DT3 <- as.integer(X$base_dist==3)
X$DT4 <- as.integer(X$base_dist==4)
X$DT5 <- as.integer(X$base_dist==5)
X$DT6_20 <- as.integer(X$base_dist>=6 & X$base_dist<20)
X$DT20_inf <- as.integer(X$base_dist>=20)
X$unmethylated_bases <-X$total_bases - X$methylated_bases
rm(cc)

for(i in 1:nrow(GGroup))
{
  #print(i)
  if(!(as.character(GGroup$TYPE[[i]]) %in% colnames(X)))
    X[[as.character(GGroup$TYPE[[i]])]]<-0
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
      set4<-which(as.integer(X$pos)>=as.integer(Y$Start[j]) & as.integer(X$pos)<=as.integer(Y$Stop[j]) & X$chrom==as.integer(Y$Replicon.Name[j]))
      if(length(set4)>0)
      {
        print(paste(as.character(Y$Protein.product[set]),"found in epi data with ",length(set4)," base pairs"))
        X[[as.character(GGroup$TYPE[[i]])]][set4]<-1
      }
    }

  }
  else
  {
    print(paste(as.character(GGroup$GENE_ID[[i]])," is not found in EPI dataset"))
  }
}
X$none <- 1
X$none <- X$none - (X$Mα + X$Mβ + X$Mβ + X$Mγ + X$Mδ + X$MIKC)

# now mark all promoters, coding regions and transcriptions
X$coding<-0
X$promoter<-0
X$aftron<-0
#define length of the promoter regions
lprom_length<-50
laftr_length<-50
for(i in 1:nrow(Y))
{
  if(i%%10==0)
    print(paste(i,"genes proceed"))
  set4<-which(as.integer(X$pos)>=as.integer(Y$Start[i]) & as.integer(X$pos)<=as.integer(Y$Stop[i]) & X$chrom==as.integer(Y$Replicon.Name[i]))
  set5<-which(X$chrom==as.integer(Y$Replicon.Name[i]))
  if(length(set4)>0)
  {

  ls4<-set4[length(set4)]
  fs4<-set4[1]
  set6<-(fs4-lprom_length):(fs4-1)
  set6<-intersect(set6, set5)
  set81<-NULL
  set82<-NULL
  if(i!=1){
    set8<-intersect(set6, set7)
    if(length(set8)>0){
      set81<-set8[1:round(length(set8)/2)]
      set82<-set8[(round(length(set8)/2)+1):length(set8)]
      if(length(set82)>0)
        X$aftron[set82]<-0
    }
  }
  set7<-(ls4+1):(ls4+laftr_length)
  set7<-intersect(set7, set5)
  #devide the sets in a clever way

    X$coding[set4]<-1
    X$promoter[set6]<-1
    if(length(set81)>0)
      X$promoter[set81]<-0
    X$aftron[set7]<-1
    X$promoter[set4]<-0
    X$aftron[set4]<-0
  }
}

write.csv(X, file = "EpigenNew.csv",row.names = F,col.names = T)

length(which(X$promoter==1))
length(which(X$aftron==1))

ss<-which(X$promoter + X$aftron == 2)
length(ss)
X$promoter[ss]<-0
X$aftron[ss]<-0
ss1<-which(X$promoter==1)
ss2<-which(X$aftron==1)
ss3<-which(X$coding==1)
length(ss1)
length(ss2)
length(ss3)

classes<-which(X$coding +  X$aftron + X$promoter > 0 & X$strand == "+")
length(classes)
classes<-which(X$coding +  X$aftron + X$promoter > 0 & X$strand == "-")
length(classes)



X$express<-0
for(i in 1:nrow(Express))
{
    if(i%%10==0)
      print(paste(i,"gene expressions proceed"))
    set4<-which(as.integer(X$pos)>=as.integer(Express$start[i]) & as.integer(X$pos)<= as.integer(Express$end[i]) & X$chrom==as.integer(Express$chrom[i]))
    if(length(set4)>0)
    {
      X$express[set4]<-Express$fpkm[i]
    }
}

data<-X[classes,]
View(data[which(data$promoter ==1),])
data$P <-data$methylated_bases/data$total_bases
data$pos <-1:dim(data)[1]


