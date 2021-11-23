simplifyposteriors<-function(X,posteriors,th=0.0001,thf=0.2, resp)
{
posteriors<-posteriors[-which(posteriors[,2]<th),]
rhash <- hash::hash()
for(i in 1:length(posteriors[,1]))
{
  expr<-posteriors[i,1]
  #print(expr)
  res<-stats::model.matrix(data=X,object = stats::as.formula(paste0(resp,"~",expr)))
  ress<-c( stringi::stri_flatten(round(res[,2],digits = 4),collapse = ""), stringi::stri_flatten(res[,2],collapse = ""),posteriors[i,2],expr)
  if(!((ress[1] %in% hash::values(rhash))))
    rhash[[ress[1]]]<-ress
  else
  {
    if(ress[1] %in% hash::keys(rhash))
    {
      rhash[[ress[1]]][3]<- (as.numeric(rhash[[ress[1]]][3]) + as.numeric(ress[3]))
      if(stringi::stri_length(rhash[[ress[1]]][4])>stringi::stri_length(expr))
        rhash[[ress[1]]][4]<-expr
    }
  }

}
res<-as.data.frame(t(hash::values(rhash)[c(3,4),]))
res$V1<-as.numeric(as.character(res$V1))
res<-res[which(res$V1>thf),]
res<-res[order(res$V1, decreasing = T),]
hash::clear(rhash)
rm(rhash)
res[which(res[,1]>1),1]<-1
colnames(res)<-c("posterior","tree")
return(res)
}
