simplifyposteriors.infer=function(X,posteriors,th=0.0000001,thf=0.5, resp)
{
#define the function simplifying logical expressions at the end of the search
  todel = which(posteriors[,2]<th)
  if(length(todel)>0)
    posteriors=posteriors[-todel,]
  rhash=hash::hash()
  for(i in 1:length(posteriors[,1]))
  {
    expr=posteriors[i,1]
    res=stats::model.matrix(data=X,object = stats::as.formula(paste0(resp,"~",expr)))
    res[,1]=res[,1]-res[,2]
    ress=c( stringi::stri_flatten(res[,1],collapse = ""), stringi::stri_flatten(res[,2],collapse = ""),posteriors[i,2],expr)
    ress[1] =  stringi::stri_sub(ress[1],from = 1,to = 9999)
    if(!(ress[1] %in% hash::values(rhash)||(ress[2] %in% hash::values(rhash))))
      rhash[[ress[1]]]=ress
    else
    {
      if(ress[1] %in% hash::keys(rhash))
      {
        rhash[[ress[1]]][3]= (as.numeric(rhash[[ress[1]]][3]) + as.numeric(ress[3]))
        if(stringi::stri_length(rhash[[ress[1]]][4])>stringi::stri_length(expr))
          rhash[[ress[1]]][4]=expr
      }
      else
      {
        rhash[[ress[2]]][3]= (as.numeric(rhash[[ress[2]]][3]) + as.numeric(ress[3]))
        if(stringi::stri_length(rhash[[ress[2]]][4])>stringi::stri_length(expr))
          rhash[[ress[2]]][4]=expr
      }
    }

  }
  res=as.data.frame(t(hash::values(rhash)[c(3,4),]))
  res$V1=as.numeric(as.character(res$V1))
  tokeep = which(res$V1>thf)
  if(length(tokeep)>0)
  {
    res=res[tokeep,]
  }else
    warning(paste0("No features with posteriors above ",thf,". Returning everything"))
  res=res[order(res$V1, decreasing = T),]
  hash::clear(rhash)
  rm(rhash)
  tokeep = which(res[,1]>1)
  if(length(tokeep)>0)
    res[tokeep,1]=1
  colnames(res)=c("posterior","feature")
  rownames(res) = NULL
  return(res)
}
