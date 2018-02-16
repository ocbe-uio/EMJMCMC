
#segment and convert pictures to data 
picToX<-function(img,slope=1,intercept=-1024,size.x=512,size.y=512,mz = 1,proj=1,bone.th = 700,bin.th = -9999,kern = makeBrush(5, shape='Gaussian',sigma = 0.3),cbrush = 5,
                 voxelspacing = 0.028, RU=1,LU=1,z=1,useCont = F)
{
  ret = T
  features = NULL
  nodules = NULL
  lungs = NULL
  tryCatch({
    #begin convert to HU scale
    #convert empty to air
    if(min(img)!=-1024)
    {
      img[img == -2000] = 0
      # Convert to Hounsfield units (HU)
      img <- img * slope + intercept
      #convert bones to air
      #end convert to HU scale
    }
    img[img > bone.th] = -1000
    #begin potential nodules detection
    #normalize the image
    curHU<-normalize(img)
    #find OSTU treshold
    
    #filter image w.r.t. OSTU treshold
    if(bin.th == -9999)
    { 
      th<-otsu(curHU, range = c(0, 1), levels = 256)
      opn <-(curHU > th)
    }else{
      opn <-img > bin.th
    }
    
    #apply opening morphological operation
    opn<-opening(x=opn,kern = kern)
    sopn<-sum(opn)
    if(sopn==0|sopn==size.x*size.y)
    {
      ret = F
      return(NULL)
    }
    #delete the regions out of the lungs
    if(!useCont){
      opn<- floodFill(x = opn, pt = c(1,1), col=1, tolerance=0)
      opn<- floodFill(x = opn, pt = c(size.x,1), col=1, tolerance=0)
      opn<- floodFill(x = opn, pt = c(size.x/2,1), col=1, tolerance=0)
      opn<- floodFill(x = opn, pt = c(1,size.y), col=1, tolerance=0)
      opn<- floodFill(x = opn, pt = c(1,size.y/2), col=1, tolerance=0)
      opn<- floodFill(x = opn, pt = c(size.x,size.y), col=1, tolerance=0)
      opn<- floodFill(x = opn, pt = c(size.x/2,size.y), col=1, tolerance=0)
      opn<- floodFill(x = opn, pt = c(size.x,size.y/2), col=1, tolerance=0)
      
    }else{
      # #find contour of the lungs
      oc <- ocontour(opn)
      #delete everything outside the contour
      pids.x<-1:size.x
      pids.y<-1:size.y
      rgn<-range(oc[[1]][,1])
      opn[which(pids.x<rgn[1]|pids.x>rgn[2]),]<-1
      rgn<-range(oc[[1]][,2])
      opn[,which(pids.y<rgn[1]|pids.y>rgn[2])]<-1
      for(i in 1:dim(oc[[1]])[1])
      {
        if(oc[[1]][i,1]<=size.x/2)
          opn[1:min(oc[[1]][i,1]+cbrush,size.x),oc[[1]][i,2]]<-1
        else
          opn[max(oc[[1]][i,1]-cbrush,1):size.x,oc[[1]][i,2]]<-1
        if(oc[[1]][i,2]<=size.y/2)
          opn[oc[[1]][i,1],1:min(oc[[1]][i,2]+cbrush,size.y)]<-1
        else
          opn[oc[[1]][i,1],max(1,oc[[1]][i,2]-cbrush):size.y]<-1
      }
    }
    
    sopn<-sum(opn)
    if(sopn==0|sopn==size.x*size.y)
    {
      ret = F
      return(NULL)
    }
    
    if(sopn/size.x/size.y>0.99)
      opn<-1-opn
    
    filled <- floodFill(x = opn, pt = c(1,1), col=0, tolerance=0)
    #end potential nodules detection
    #image(filled)
    #begin feature extraction
    sf<-sum(filled)
    if(is.null(filled)|is.null(img)|is.null(curHU)|sf== (size.x)*(size.y) | sf==0)
    {
      ret = F
    }else{
      features<-as.array(c(computeFeatures.shape(x = filled,properties = F),
                           computeFeatures.moment(x = filled, ref = img),
                           computeFeatures.moment(x = filled, ref = curHU),
                           computeFeatures.basic(x = filled, ref = img),
                           computeFeatures.basic(x = filled, ref = curHU),
                           computeFeatures.haralick(x = filled,ref = img),
                           computeFeatures.haralick(x = filled,ref = curHU)
      ))
    }
    
    mx<-features[7]
    my<-features[8]
    
    #print(mx)
    #print(my)
    
    if(length(mx)==0|length(my)==0)
    {
      ref = F
      return(NULL)
    }
    #print(mx)
    #print(my)
    
    NL<-0
    if(proj==1)
    {
      
      if((my<size.y/2-10)&(LU))
      {
        NL<-1
      }else{
        if((my>size.y/2+10)&(RU))
          NL<-1
      }
    }else{
      if(proj==2)
      {
        if((mx<size.x/2-10)&my>size.y/2)
        {
          NL<-1
        }else{
          if((mx>size.x/2+10)&my>size.y/3)
            NL<-1
        }
      }else
      {
        NL<-z*(RU*(my>size.y/3)|RU*(my>size.y/2))
      }
    }
    
    NS1<-(((features[3]*voxelspacing/10)^(-0.5))-1.58113883)
    m1.max<- (-6.5929+0.5806 - 5.86116*NS1 + 0.64*NL)
    
    NS2<-(((features[5]*voxelspacing/10)^(-0.5))-1.58113883)
    m2.max<- (-6.5929+0.5806 - 5.86116*NS1 + 0.64*NL)
    
    NS3<-(((features[6]*voxelspacing/10)^(-0.5))-1.58113883)
    m3.max<- (-6.5929+0.5806 - 5.86116*NS3 + 0.64*NL)
    
    #print(c(NS1,NS2,NS3))
    #print(c(features[3],features[5],features[6]))
    #print(c(1/(1+exp(-m1.max)),1/(1+exp(-m2.max)),1/(1+exp(-m3.max))))
    
    features<-c(m1.max,m2.max,m3.max,NS1,NS2,NS3,NL,RU,LU,features)
    
    
    #end feature extraction
    if(useCont){
      rm(oc)
      rm(pids.x)
      rm(pids.y)
    }
    rm(sf)
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

