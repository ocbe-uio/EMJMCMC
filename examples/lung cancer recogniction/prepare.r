#n/b these packages MUST be installed prior to analysis!
#install.packages("oro.dicom")
#install.packages("oro.nifti")
#source("http://bioconductor.org/biocLite.R")
#biocLite("EBImage")

library("oro.dicom")
library("oro.nifti")
library("geometry")
require(EBImage)
library(parallel)

patientDir <- sort(list.dirs(path = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/kaggle3/sample_images/sample_images", full.names = TRUE, recursive = FALSE))
#patientDir<-patientDir[1:6]
picdir <- "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/kaggle3/scans/"
pat.names <- sort(list.dirs(path = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/kaggle3/sample_images/sample_images", full.names = F, recursive = FALSE))
labels<-read.csv(file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/kaggle3/labels known.csv")
data.sample<-matrix(data = NA,nrow = length(patientDir),ncol = 85)


# vectorized call of picToX function
call.picToX<-function(vect)
{
  do.call(picToX,vect)
}

#patient i proceeding function
proceed<-function(pat.id)
{
  dicom <- readDICOM(patientDir[pat.id])
  dicomdf <- dicomTable(dicom$hdr)
  
  vect<-array(data = list(),dim = length(dicom$img))
  for(i in 1:length(dicom$img))
  {
    vect[[i]]$img = dicom$img[[i]]
    vect[[i]]$slope = as.numeric(dicom$hdr[[i]]$value[dicom$hdr[[i]]$name == "RescaleSlope"])
    vect[[i]]$intercept = as.numeric(dicom$hdr[[i]]$value[dicom$hdr[[i]]$name == "RescaleIntercept"])
    vect[[i]]$bone.th = 1000
    vect[[i]]$kern =makeBrush(5, shape='Gaussian',sigma = 0.3)
    vect[[i]]$size = as.numeric(dicomdf$`0028-0010-Rows`[i])
    vect[[i]]$cbrush = 5
  }
  
  # vectorized (on cpu only) filtration of all of the 3d images for a given patient (map)
  res<-lapply(X = vect,FUN = call.picToX)
  
  #analyze the results (reduce) and choose the slice with the largest area of the tumor
  # or average the results or choose other objective
  rmax<-0
  imax<-1
  for(i in 1:length(dicom$img))
  {
    if(is.null(res[[i]]))
      next
    if(is.null(res[[i]]$X))
      next
    rcur<-res[[i]]$X[1]
    if(rcur>rmax & rcur < 5000)
    {
      rmax = rcur
      imax = i
    }
  }
  
  stat<-labels[which(labels[,1]==pat.names[pat.id]),2]
  if(length(stat)==0)
    stat = "999"
  
  if(doPNG <- !dev.interactive())
  png(paste0(picdir,"nodules_",pat.id,"_",stat,".png"), width=vect[[imax]]$size, height=vect[[imax]]$size)
  image(res[[imax]]$nodules, useRaster = F,col = grey(0:64 / 64)) # should not suffer from anti-aliasing
  if(doPNG)
    dev.off()
  if(doPNG <- !dev.interactive())
    png(paste0(picdir,"lungs_",pat.id,"_",stat,".png"), width=vect[[imax]]$size, height=vect[[imax]]$size)
  image(res[[imax]]$lungs, useRaster = F,col = grey(0:64 / 64)) # should not suffer from anti-aliasing
  if(doPNG)
    dev.off()
  if(doPNG <- !dev.interactive())
    png(paste0(picdir,"img_",pat.id,"_",stat,".png"), width=vect[[imax]]$size, height=vect[[imax]]$size)
  image(vect[[imax]]$img, useRaster = F,col = grey(0:64 / 64)) # should not suffer from anti-aliasing
  if(doPNG)
    dev.off()
  
  rm(vect)
  rm(dicom)
  rm(dicomdf)
  
  return(list(X = c(res[[imax]]$X,stat), imax = imax,lungs = res[[imax]]$lungs,nodules = res[[imax]]$nodules))
}

data = mclapply(FUN = proceed,X = 1:length(patientDir),mc.preschedule = F,mc.cores = 4,mc.cleanup = T)

for(i in 1:length(patientDir))
{
  if(is.null(data[[i]]))
  {
    warning(paste0("Patient ",i,"not proceeded correctly"))
    next
  }
  data.sample[i,]<-data[[i]]$X
  
  
}

write.csv(file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/kaggle3/data.csv",x = data.sample)

