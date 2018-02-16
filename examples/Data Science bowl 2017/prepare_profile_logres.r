#n/b these packages MUST be installed prior to analysis!
#install.packages("oro.dicom")
#install.packages("oro.nifti")
#source("http://bioconductor.org/biocLite.R")
#biocLite("EBImage")
#install.packages("stringi")
library("oro.dicom")
library("oro.nifti")
library("stringi")
require(EBImage)
library(parallel)

scansPath <- "scans/" 
picdir <- "results/scans/profile/"
resultLabelsPath<-"results/labels/profile/"
resultLabelsFile = paste0(resultLabelsPath, "data_profile.csv")
labelsPath <- "labels.csv"
patientLabelsPath <- "results/labels/profile/"

dir.create(picdir, recursive = TRUE)
dir.create(resultLabelsPath, recursive = TRUE)
dir.create(patientLabelsPath, recursive = TRUE)

patientDir <- sort(list.dirs(path = scansPath, full.names = TRUE, recursive = FALSE))
pat.names <- sort(list.dirs(path = scansPath, full.names = F, recursive = FALSE))
labels<-read.csv(file = labelsPath)
data.sample<-matrix(data = NA,nrow = length(patientDir),ncol = 113)


# vectorized call of picToX function
call.picToX<-function(vect)
{
  do.call(picToX,vect)
}
g<-function(x)
{
  return((x = 1/(1+exp(-x))))
}
#patient i proceeding function
proceed.head<-function(pat.id)
{
  print(pat.id)
  tryCatch({
  dicom <- readDICOM(patientDir[pat.id])
  dicomdf <- dicomTable(dicom$hdr)
  ps<-as.numeric(stri_split_fixed(dicomdf$`0028-0030-PixelSpacing`[1],pattern =  " ")[[1]][1])
  sl<-sort(as.numeric((dicomdf$`0020-1041-SliceLocation`)))
  sln<-as.numeric((dicomdf$`0020-1041-SliceLocation`))
  ss<-sl[2]-sl[1]
  
  if(is.na(ss))
  {
    poss<-stri_split_fixed(dicomdf$`0020-0032-ImagePositionPatient`,pattern =  " ")
    sln<-array(0,dim = length(poss))
    for(i in 1:length(poss))
    {
      sln[i]<-as.numeric(poss[[i]][3])
    }
    sl<-sort(sln)
    ss<-sl[2]-sl[1]
  }
  

  dicom2nifti <- dicom2nifti(dicom,rescale = T, reslice = T)
  hval = hist(dicom2nifti,breaks = 8)$density[1:7]
  slope = scl_slope(dicom2nifti)
  intercept = scl_inter(dicom2nifti)
  
  vect<-array(data = list(),dim = dim(dicom2nifti)[1])
  for(i in 1:dim(dicom2nifti)[1])
  {
    vect[[i]]$img = dicom2nifti[i,,]
    vect[[i]]$slope = as.numeric(slope)
    vect[[i]]$intercept = as.numeric(intercept)
    vect[[i]]$bone.th = 1000
    vect[[i]]$kern =makeBrush(5, shape='Gaussian',sigma = 0.3)
    vect[[i]]$size.x = dim(dicom2nifti)[2]
    vect[[i]]$size.y = dim(dicom2nifti)[3]
    vect[[i]]$cbrush = 5
    vect[[i]]$bin.th = -9999
    vect[[i]]$RU=0
    vect[[i]]$LU=0
    if(i < dim(dicom2nifti)[1]/2)
      vect[[i]]$RU=1
    else
      vect[[i]]$LU=1
    vect[[i]]$voxelspacing = (ps+ss)/2
    vect[[i]]$z=3
    #print( vect[[i]]$RU)
    #print( vect[[i]]$LU)
  }

  # vectorized (on cpu only) filtration of all of the 3d images for a given patient (map)
  
  
  res<-lapply(X = vect,FUN = call.picToX)
  
  stat<-labels[which(labels[,1]==pat.names[pat.id]),2]
  if(length(stat)==0)
    stat = "999"
  #analyze the results (reduce) and choose the slice with the largest area of the tumor
  # or average the results or choose other objective
  ocsv<-NULL
  rmax<-c(-Inf,-Inf,-Inf)
  imax<-1
  for(i in 1:length(res))
  {
    if(length(res[[i]])==0)
      next
    if(length(res[[i]]$X)==0)
      next
    
    ocsv<-cbind(ocsv,c(i,hval,g(res[[i]]$X[1]),g(res[[i]]$X[2]),g(res[[i]]$X[3]),sum(sl>sln[i])/length(sl),res[[i]]$X,ps,ps^2,1/ps,1/ps^2,ss,1/ss,ss*ps,1/ss/ps,stat))
    rcur<-res[[i]]$X[c(1,2,3,10)]
    if(rcur[1]>rmax[1] & rcur[2]>rmax[2] & rcur[3]>rmax[3] & rcur[4] < 2500)
    {
      rmax = rcur[-4]
      imax = i
    }
  }
  
  if(length(ocsv)!=0)
    write.csv(file = paste0(patientLabelsPath, "data_profile_",pat.id,"_",imax,".csv"),x = ocsv)
  
  if(length(res[[imax]])==0|length(res[[imax]]$lungs)==0|length(res[[imax]]$nodules)==0|length(res[[imax]]$lungs)==0)
  {
    warning(paste0("Patient ",pat.id,"not proceeded correctly"))
    rm(vect)
    rm(dicom)
    rm(dicomdf)
    return(NULL)
  }
  
  if(doPNG <- !dev.interactive())
  png(paste0(picdir,"nodules_p_",pat.id,"_",stat,".png"), width=512, height=512)
  image(res[[imax]]$nodules, useRaster = F,col = grey(0:64 / 64)) # should not suffer from anti-aliasing
  if(doPNG)
    dev.off()
  if(doPNG <- !dev.interactive())
    png(paste0(picdir,"lungs_p_",pat.id,"_",stat,".png"), width=512, height=512)
  image(res[[imax]]$lungs, useRaster = F,col = grey(0:64 / 64)) # should not suffer from anti-aliasing
  if(doPNG)
    dev.off()
  if(doPNG <- !dev.interactive())
    png(paste0(picdir,"img_p_",pat.id,"_",stat,".png"), width=512, height=512)
  image(vect[[imax]]$img, useRaster = F,col = grey(0:64 / 64)) # should not suffer from anti-aliasing
  if(doPNG)
    dev.off()
  
  rm(vect)
  rm(dicom)
  rm(dicomdf)
  
  return(list(X = c(g(res[[imax]]$X[1]),g(res[[imax]]$X[2]),g(res[[imax]]$X[3]),hval,sum(sl>sln[imax])/length(sl),res[[imax]]$X,ps,ps^2,1/ps,1/ps^2,ss,1/ss,ss*ps,1/ss/ps,stat), imax = imax,lungs = res[[imax]]$lungs,nodules = res[[imax]]$nodules))},
  error=function(e)
  {
    print(e)
    return(NULL)
  })
}

switch(Sys.info()[['sysname']],
       Windows= {data = lapply(FUN = proceed.head,X = 1:length(patientDir))},
       Linux  = {data = mclapply(FUN = proceed.head,X = 1:length(patientDir),mc.preschedule = F,mc.cores = 4,mc.cleanup = T)},
       Darwin = {data = lapply(FUN = proceed.head,X = 1:length(patientDir))})

for(i in 1:length(patientDir))
{
  if(length(data[[i]])==0|length(data[[i]]$X)==0)
  {
    warning(paste0("Patient ",i,"not proceeded correctly"))
    next
  }
  data.sample[i,]<-data[[i]]$X
  #print(data[[i]]$imax)
  
  
}


write.csv(file = resultLabelsFile,x = data.sample)




