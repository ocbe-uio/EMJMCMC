library(hash)
library(stringi)
library(qtl)


setwd("/nr/samba/user/ahu/Documents/tmp/BLRRES")
#read the output files of all of the MM simulations for the scenario of interest
temp = list.files(pattern="post2etaJQTL_*")
myfiles = lapply(FUN = read.csv,X = temp)


scenairo=5
#create a dummy matrix of binary covariates
X3= as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))
rhash=hash()
Z3 = X3
if(scenairo == 5)
{
  tps=array(data = 0, dim =c(1000,4))
  tps[,1]= X3$V4*X3$V17*X3$V30*X3$V10
  tps[,2]= X3$V7*X3$V20*X3$V12
  tps[,3]=(X3$V2*X3$V9)
  tps[,4]=(X3$V37)
}

#analyse the results of the simulations
N=length(myfiles)
alpha=0.5
clear(rhash)
fdr=0

WL = 0
for(i in 1:min(100,N))
{
  tpi=0
  fpi=0
  
  
  set.seed(i)
  
  map = sim.map(c(100,90,80,60,40), 50/5, include.x=FALSE, eq.spacing=TRUE)
  #plotMap(map)
  
  # The longer the chromosomal length the less correlated are markers
  # (at least theoretically)
  
  # Now simulate data from a backcross design using the sim.cross function
  n.ind = 1000   #sample size
  
  simbc = sim.cross(map, type="bc", n.ind=n.ind,
                    model=rbind(c(1,45,1), c(5,20,1), c(5,50,1)))
  
  # The relevant genotype data is in the structure geno 
  #str(simbc$geno)   #it is a bit complicated
  
  # Get an X matrix for each chromosome
  X.list = vector("list", 5)
  for (chr in 1:5){
    X.list[[chr]] = pull.geno(simbc, chr)
  }
  
  #Check the correlations between markers within the same chromosome
  lapply(X.list, cor)
  
  
  #Create one large X matrix which you can the use to make your own 
  # simulations of Y data with a logic regression model
  X = cbind(X.list[[1]],X.list[[2]],X.list[[3]],X.list[[4]],X.list[[5]])-1
  #permute elements of X
  X = X[,sample.int(n = 50,size = 50,replace = F)]
  
  map = unlist(map)
  lmap = stringi::stri_sub(labels(map),3)
  map = as.data.frame(t(map))
  names(map) = lmap
  
  map[1:10] = map[1:10] + 100
  map[11:20] = map[1:10] + 200
  map[21:30] = map[1:10] + 300
  map[31:40] = map[1:10] + 400
  map[41:50] = map[1:10] + 500
  
  X2=as.data.frame(X)
  names.old = names(X2)
  names(X2) = paste0("V",1:50)
  names.new = names(X2)
  
  X = X2
  


  names.all = as.data.frame(cbind(names.new,names.old))
  
  names.all$n1 = names.new
  names.all$n2 = names.new
  names.all$n3 = names.new
  names.all$n4 = names.new
  names.all$n5 = names.new
  names.all$n6 = names.new
  for(ncur in names.old[c(4,17,30,10,7,20,12,2,9,37)])
  {
    
    chr = as.integer(stringi::stri_sub(ncur,2,2))
    mar = as.integer(stringi::stri_sub(ncur,4,4))
    if(chr<4)
    {
      if(mar==10)
      {
        nab1 = paste0("D",chr,"M",9)
        nab2 = paste0("D",chr,"M",9)
        nab3 = paste0("D",chr,"M",9)
        nab4 = paste0("D",chr,"M",9)
        nab5 = paste0("D",chr,"M",9)
        nab6 = paste0("D",chr,"M",9)
      }else if(mar==1)
      {
        nab1 = paste0("D",chr,"M",2)
        nab2 = paste0("D",chr,"M",2)
        nab3 = paste0("D",chr,"M",2)
        nab4 = paste0("D",chr,"M",2)
        nab5 = paste0("D",chr,"M",2)
        nab6 = paste0("D",chr,"M",2)
      }else
      {
        nab1 = paste0("D",chr,"M",mar+1)
        nab2 = paste0("D",chr,"M",mar-1)
        nab3 = paste0("D",chr,"M",mar+1)
        nab4 = paste0("D",chr,"M",mar-1)
        nab5 = paste0("D",chr,"M",mar-1)
        nab6 = paste0("D",chr,"M",mar-1)
      }
      
    }else if(chr==4){
      
      if(mar==10)
      {
        nab1 = paste0("D",chr,"M",9)
        nab2 = paste0("D",chr,"M",8)
        nab3 = paste0("D",chr,"M",8)
        nab4 = paste0("D",chr,"M",9)
        nab5 = paste0("D",chr,"M",8)
        nab6 = paste0("D",chr,"M",9)
      }else if(mar==1)
      {
        nab1 = paste0("D",chr,"M",2)
        nab2 = paste0("D",chr,"M",3)
        nab3 = paste0("D",chr,"M",2)
        nab4 = paste0("D",chr,"M",3)
        nab5 = paste0("D",chr,"M",2)
        nab6 = paste0("D",chr,"M",3)
      }else if(mar==9)
      {
        nab1 = paste0("D",chr,"M",10)
        nab2 = paste0("D",chr,"M",10)
        nab3 = paste0("D",chr,"M",8)
        nab4 = paste0("D",chr,"M",7)
        nab5 = paste0("D",chr,"M",8)
        nab6 = paste0("D",chr,"M",7)
      }else if(mar==2)
      {
        nab1 = paste0("D",chr,"M",1)
        nab2 = paste0("D",chr,"M",1)
        nab3 = paste0("D",chr,"M",3)
        nab4 = paste0("D",chr,"M",4)
        nab5 = paste0("D",chr,"M",3)
        nab6 = paste0("D",chr,"M",4)
      }
      else
      {
        nab1 = paste0("D",chr,"M",mar+1)
        nab2 = paste0("D",chr,"M",mar-1)
        nab3 = paste0("D",chr,"M",mar+2)
        nab4 = paste0("D",chr,"M",mar-2)
        nab5 = paste0("D",chr,"M",mar+2)
        nab6 = paste0("D",chr,"M",mar-2)
      }
    }else {
      
      if(mar==10)
      {
        nab1 = paste0("D",chr,"M",9)
        nab2 = paste0("D",chr,"M",8)
        nab3 = paste0("D",chr,"M",8)
        nab4 = paste0("D",chr,"M",7)
        nab5 = paste0("D",chr,"M",7)
        nab6 = paste0("D",chr,"M",7)
      }else if(mar==1)
      {
        nab1 = paste0("D",chr,"M",2)
        nab2 = paste0("D",chr,"M",3)
        nab3 = paste0("D",chr,"M",2)
        nab4 = paste0("D",chr,"M",4)
        nab5 = paste0("D",chr,"M",4)
        nab6 = paste0("D",chr,"M",4)
      }else if(mar==9)
      {
        nab1 = paste0("D",chr,"M",10)
        nab2 = paste0("D",chr,"M",10)
        nab3 = paste0("D",chr,"M",8)
        nab4 = paste0("D",chr,"M",7)
        nab5 = paste0("D",chr,"M",6)
        nab6 = paste0("D",chr,"M",6)
      }else if(mar==2)
      {
        nab1 = paste0("D",chr,"M",1)
        nab2 = paste0("D",chr,"M",1)
        nab3 = paste0("D",chr,"M",3)
        nab4 = paste0("D",chr,"M",4)
        nab5 = paste0("D",chr,"M",5)
        nab6 = paste0("D",chr,"M",5)
      }
      else if(mar==3){
        nab1 = paste0("D",chr,"M",1)
        nab2 = paste0("D",chr,"M",2)
        nab3 = paste0("D",chr,"M",4)
        nab4 = paste0("D",chr,"M",5)
        nab5 = paste0("D",chr,"M",6)
        nab6 = paste0("D",chr,"M",6)
      }else if(mar==8){
        nab1 = paste0("D",chr,"M",10)
        nab2 = paste0("D",chr,"M",9)
        nab3 = paste0("D",chr,"M",7)
        nab4 = paste0("D",chr,"M",6)
        nab5 = paste0("D",chr,"M",5)
        nab6 = paste0("D",chr,"M",5)
      }else{
        nab1 = paste0("D",chr,"M",mar+1)
        nab2 = paste0("D",chr,"M",mar-1)
        nab3 = paste0("D",chr,"M",mar+2)
        nab4 = paste0("D",chr,"M",mar-2)
        nab5 = paste0("D",chr,"M",mar+3)
        nab6 = paste0("D",chr,"M",mar-3)
      }
      
    }
    
    names.all$n1[which(names.all$names.old==ncur)]=as.character(names.all$names.new[which(names.all$names.old==nab1)])
    names.all$n2[which(names.all$names.old==ncur)]=as.character(names.all$names.new[which(names.all$names.old==nab2)])
    names.all$n3[which(names.all$names.old==ncur)]=as.character(names.all$names.new[which(names.all$names.old==nab3)])
    names.all$n4[which(names.all$names.old==ncur)]=as.character(names.all$names.new[which(names.all$names.old==nab4)])
    #names(map[which(abs(map[ncur]-map)<=10)])
    
  }
  names.all$names.old=NULL
  X = Z3
  
  str.all = unique(stringi::stri_paste(as.matrix(names.all[c(37,2,9,7,20,12,4,17,30,10),]),sep = "  "))
  
  for(j in 1:length(myfiles[[i]]$posterior))
  {

    if(myfiles[[i]]$posterior[j]>=alpha)
    {
      expr=as.character(myfiles[[i]]$tree[j])
      expr.paresed = stringi::stri_replace_all_fixed(str = expr,pattern = "(","")
      expr.paresed = stringi::stri_replace_all_fixed(str = expr.paresed,pattern = ")","")
      expr.paresed = stringi::stri_replace_all_fixed(str = expr.paresed,pattern = "I","")
      expr.paresed = stringi::stri_replace_all_fixed(str = expr.paresed,pattern = "|","&")
      
      expr.paresed = stringi::stri_split_fixed(expr.paresed,"&")[[1]]
      #print(expr.paresed )
      WL= WL + sum(!expr.paresed %in% str.all)
      expr.old = expr
      count = stringi::stri_count_fixed(str = expr,pattern = "V")
      #print("before")
      #print(expr)
      if(count == 1 & sum(expr.paresed%in%"V37")==0)
      for(ncur in c(37))
      {
        if(!names.all$n1[ncur]%in%names.new[c(37)])
        {  
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n1[ncur],")"),paste0(names.all$names.new[ncur],")"))
        }
          #else
        #  expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n1[ncur],")"),paste0(names.all$names.new[ncur],")"))
        if(!names.all$n2[ncur]%in%names.new[c(37)])
        {
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n2[ncur],")"),paste0(names.all$names.new[ncur],")"))
        }
          #else
        #  expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n2[ncur],")"),paste0(names.all$names.new[ncur],")"))
        if(!names.all$n3[ncur]%in%names.new[c(37)])
        {
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n3[ncur],")"),paste0(names.all$names.new[ncur],")"))
        }  
        #else
        #  expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n3[ncur],")"),paste0(names.all$names.new[ncur],")"))
        if(!names.all$n4[ncur]%in%names.new[c(37)])
        {
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n4[ncur],")"),paste0(names.all$names.new[ncur],")"))
        }
        if(!names.all$n5[ncur]%in%names.new[c(37)])
        {
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n5[ncur],")"),paste0(names.all$names.new[ncur],")"))
        }
        if(!names.all$n6[ncur]%in%names.new[c(37)])
        {
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n6[ncur],")"),paste0(names.all$names.new[ncur],")"))
        }
        #else
        #  expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n4[ncur],")"),paste0(names.all$names.new[ncur],")"))
      }
      
      if(count == 2 & sum(expr.paresed%in%c("V2","V9"))<2)
      for(ncur in c(2,9))
      {
        if(!names.all$n1[ncur]%in%names.new[ncur])
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n1[ncur],")"),paste0(names.all$names.new[ncur],")"))
        if(!names.all$n2[ncur]%in%names.new[ncur])
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n2[ncur],")"),paste0(names.all$names.new[ncur],")"))
        if(!names.all$n3[ncur]%in%names.new[ncur])
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n3[ncur],")"),paste0(names.all$names.new[ncur],")"))
        if(!names.all$n4[ncur]%in%names.new[ncur])
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n4[ncur],")"),paste0(names.all$names.new[ncur],")"))
        if(!names.all$n5[ncur]%in%names.new[ncur])
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n5[ncur],")"),paste0(names.all$names.new[ncur],")"))
        if(!names.all$n6[ncur]%in%names.new[ncur])
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n6[ncur],")"),paste0(names.all$names.new[ncur],")"))
        #else
        #  expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n4[ncur],")"),paste0(names.all$names.new[ncur],")"))
      }
      
      if(count == 3 & sum(expr.paresed%in%c("V7","V20","V12"))<3)
      for(ncur in c(7,20,12))
      {
        if(!names.all$n1[ncur]%in%names.new[ncur])
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n1[ncur],")"),paste0(names.all$names.new[ncur],")"))
        if(!names.all$n2[ncur]%in%names.new[ncur])
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n2[ncur],")"),paste0(names.all$names.new[ncur],")"))
        if(!names.all$n3[ncur]%in%names.new[ncur])
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n3[ncur],")"),paste0(names.all$names.new[ncur],")"))
        if(!names.all$n4[ncur]%in%names.new[ncur])
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n4[ncur],")"),paste0(names.all$names.new[ncur],")"))
        if(!names.all$n5[ncur]%in%names.new[ncur])
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n5[ncur],")"),paste0(names.all$names.new[ncur],")"))
        if(!names.all$n6[ncur]%in%names.new[ncur])
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n6[ncur],")"),paste0(names.all$names.new[ncur],")"))
        #else
        #  expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n4[ncur],")"),paste0(names.all$names.new[ncur],")"))
      }
      
      if(count == 4 & sum(expr.paresed%in%c("V4","V17","V10","V30"))<3)
      for(ncur in c(4,17,30,10))
      {
        if(!names.all$n1[ncur]%in%names.new[ncur])
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n1[ncur],")"),paste0(names.all$names.new[ncur],")"))
        if(!names.all$n2[ncur]%in%names.new[ncur])
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n2[ncur],")"),paste0(names.all$names.new[ncur],")"))
        if(!names.all$n3[ncur]%in%names.new[ncur])
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n3[ncur],")"),paste0(names.all$names.new[ncur],")"))
        if(!names.all$n4[ncur]%in%names.new[ncur])
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n4[ncur],")"),paste0(names.all$names.new[ncur],")"))
        if(!names.all$n5[ncur]%in%names.new[ncur])
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n5[ncur],")"),paste0(names.all$names.new[ncur],")"))
        if(!names.all$n6[ncur]%in%names.new[ncur])
          expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n6[ncur],")"),paste0(names.all$names.new[ncur],")"))
        #else
        #  expr = stringi::stri_replace_first_fixed(expr,pattern = paste0(names.all$n4[ncur],")"),paste0(names.all$names.new[ncur],")"))
      }
      #print("after")
      #print(expr)
      
      if(expr.old!=expr)
      {
        print(paste0("In run ",i," Expression ",expr.old," was replaced with ", expr))
      }
      

      
      res=model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
      res[,1]=res[,1]-res[,2]
      TP=0
      for(u in 1:dim(tps)[2])
      {
        if(sum(tps[,u]-res[,2])==0)
        {
          TP=1
          break
        }
      }
      if(TP)
      {
        tpi=tpi+1
      }else{
        fpi=fpi+1
      }
      
      ress=c(stri_flatten(res[,1],collapse = ""),stri_flatten(res[,2],collapse = ""),1,expr)
      if(!(ress[1] %in% values(rhash)||(ress[2] %in% values(rhash))))
        rhash[[ress[1]]]=ress
      else
      {
        if(ress[1] %in% keys(rhash))
        {
          rhash[[ress[1]]][3]= (as.numeric(rhash[[ress[1]]][3])) + as.numeric(1)
          if(stri_length(rhash[[ress[1]]][4])>stri_length(expr))
            rhash[[ress[1]]][4]=expr
        }
        else
        {
          rhash[[ress[2]]][3]= (as.numeric(rhash[[ress[2]]][3])) + as.numeric(1)
          if(stri_length(rhash[[ress[2]]][4])>stri_length(expr))
            rhash[[ress[2]]][4]=expr
        }
      }
    }
  }
#  print(fpi)
#  print(tpi)
  if((fpi+tpi)>0)
    fdr=fdr+(fpi/(fpi+tpi))
}

#print the FDR out
print("FDR:")
print(fdr/100)
print("WL:")
print(WL)

#write the frequencies of all significant discoveries
res=as.data.frame(t(values(rhash)[c(4,3),]))
res$V1=as.numeric(as.character(res$V1))
res=res[order(res$V1, decreasing = T),]
colnames(res)=c("posterior","tree")
write.csv(x = t(values(rhash))[,c(3,4)],file = "found_corr.csv",row.names = F,col.names = F)
