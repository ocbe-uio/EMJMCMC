picToX<-function(img,slope=1,intercept=1,size=512,bone.th = 700,kern = makeBrush(5, shape='Gaussian',sigma = 0.3),cbrush = 5)
{
  ret = T
  features = NULL
  nodules = NULL
  lungs = NULL
  tryCatch({
  #begin convert to HU scale
  #convert empty to water
  img[img == -2000] = 0
  # Convert to Hounsfield units (HU)
  img <- img * slope + intercept
  #convert bones to water
  img[img > bone.th] = -1000
  #end convert to HU scale
  
  #begin potential nodules detection
  #normalize the image
  curHU<-normalize(img)
  #find OSTU treshold
  th<-otsu(curHU, range = c(0, 1), levels = 256)
  #filter image w.r.t. OSTU treshold
  opn <-(curHU > th)
  #apply opening morphological operation
  opn<-opening(x=opn,kern = kern)
  
  #find contour of the lungs
  oc <- ocontour(opn)
  #delete everything outside the contour
  pids<-1:size
  rgn<-range(oc[[1]][,1])
  opn[which(pids<rgn[1]|pids>rgn[2]),]<-1
  rgn<-range(oc[[1]][,2])
  opn[,which(pids<rgn[1]|pids>rgn[2])]<-1
  for(i in 1:dim(oc[[1]])[1])
  {
    if(oc[[1]][i,1]<=size/2)
      opn[1:min(oc[[1]][i,1]+cbrush,size),oc[[1]][i,2]]<-1
    else
      opn[max(oc[[1]][i,1]-cbrush,1):size,oc[[1]][i,2]]<-1
    if(oc[[1]][i,2]<=size/2)
      opn[oc[[1]][i,1],1:min(oc[[1]][i,2]+cbrush,size)]<-1
    else
      opn[oc[[1]][i,1],max(1,oc[[1]][i,2]-cbrush):size]<-1
  }
  
  filled <- floodFill(x = opn, pt = c(1,1), col=0, tolerance=0)
  #end potential nodules detection
  
  #begin feature extraction
  
  if(is.null(filled)|is.null(img)|is.null(curHU))
    ret = F
  else
    features<-as.array(c(computeFeatures.shape(x = filled,properties = F),
                  computeFeatures.basic(x = filled, ref = img),
                  computeFeatures.basic(x = filled, ref = curHU),
                  computeFeatures.moment(x = filled, ref = img),
                  computeFeatures.moment(x = filled, ref = curHU),
                  computeFeatures.haralick(x = filled,ref = img),
                  computeFeatures.haralick(x = filled,ref = curHU)
  ))
  
  #end feature extraction
  rm(oc)
  rm(pids)
  
  }, 
  error=function(e)
  {
    print(e)
    ret = F
  },
  finally={
    if(ret)
      return(list(X=features, nodules = filled,lungs = opn))
    else
      return(NULL)
  }
  )
}

