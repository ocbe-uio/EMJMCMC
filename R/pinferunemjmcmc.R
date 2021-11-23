pinferunemjmcmc = function(n.cores = 4, mcgmj = mcgmjpse, report.level =  0.5, simplify = F, num.mod.best = 1000, predict = F,  test.data = 1, link.function = function(z)z, runemjmcmc.params)
{

  if(predict)
  {
    runemjmcmc.params$save.beta = T

    if(length(test.data)==0)
    {
      warning("Test data is not provided. No predictions will be made!")
    }
  }

  params = list(runemjmcmc.params)[rep(1,n.cores)]
  for(i in 1:n.cores)
  {
    params[[i]]$test = test.data
    params[[i]]$link = link.function
    params[[i]]$predict = predict
    params[[i]]$NM = num.mod.best
    params[[i]]$cpu=i
  }
  M = n.cores
  #results = runpar.infer(params[[1]])
  results=mcgmj(X = params,FUN = runpar.infer,mc.cores = n.cores)
  #print(results)
  #clean up
  gc()
  #prepare the data structures for final analysis of the runs
  compmax = runemjmcmc.params$interact.param$Nvars.max + 1
  resa=array(data = 0,dim = c(compmax,M*3))
  post.popul = array(0,M)
  max.popul = array(0,M)
  nulls=NULL
  not.null=1
  #check which threads had non-zero exit status
  for(k in 1:M)
  {
    if(length(results[[k]])<=1||length(results[[k]]$cterm)==0||length(results[[k]]$p.post)!=runemjmcmc.params$interact.param$Nvars.max)
    {
      nulls=c(nulls,k)
      #warning(paste0("Thread ",k,"did not converge or was killed by OS!"))
      next
    }
    else
    {
      not.null = k
    }

  }

  if(length(nulls) == M)
  {
    warning("All threads did not converge or gave an error! Returning stats from the threads only!")
    return(list(feat.stat = NULL,predictions = NULL,allposteriors = NULL, threads.stats = results))
  }


  #for all of the successful runs collect the results into the corresponding data structures
  for(k in 1:M)
  {
    if(k %in% nulls)
    {
      results[[k]]=results[[not.null]]
    }
    max.popul[k]=results[[k]]$cterm
    post.popul[k]=results[[k]]$post.populi
    resa[,k*3-2]=c(results[[k]]$fparam,"Post.Gen.Max")
    resa[,k*3-1]=c(results[[k]]$p.post,results[[k]]$cterm)
    resa[,k*3]=rep(post.popul[k],length(results[[k]]$p.post)+1)

  }
  #renormalize estimates of the marginal inclusion probabilities
  #based on all of the runs
  ml.max=max(max.popul)
  post.popul=post.popul*exp(-ml.max+max.popul)
  p.gen.post=post.popul/sum(post.popul)

  #perform BMA of the redictions across the runs
  pred = NULL
  if(predict){
    pred = results[[1]]$preds*p.gen.post[1]
    if(M > 1) {

      for(i in 2:M)
      {

        pred=pred+results[[i]]$preds*p.gen.post[i]

      }
    }
  }
  hfinal=hash::hash()
  for(ii in 1:M)
  {
    resa[,ii*3]=p.gen.post[ii]*as.numeric(resa[,ii*3-1])
    resa[length(resa[,ii*3]),ii*3]=p.gen.post[ii]
    if(p.gen.post[ii]>0)
    {
      for(jj in 1:(length(resa[,ii*3])-1))
      {
        if(resa[jj,ii*3]>0)
        {
          if(as.integer(hash::has.key(hash = hfinal,key =resa[jj,ii*3-2]))==0)
            hfinal[[resa[jj,ii*3-2]]]=as.numeric(resa[jj,ii*3])
          else
            hfinal[[resa[jj,ii*3-2]]]=hfinal[[resa[jj,ii*3-2]]]+as.numeric(resa[jj,ii*3])
        }

      }
    }
  }

  posteriors=hash::values(hfinal)
  hash::clear(hfinal)
  #delete the unused further variables
  rm(hfinal)
  rm(resa)
  rm(post.popul)
  rm(max.popul)
  #simplify the found trees and their posteriors
  posteriors=as.data.frame(posteriors)
  posteriors=data.frame(X=row.names(posteriors),x=posteriors$posteriors)
  posteriors$X=as.character(posteriors$X)
  res1 = NULL
  if(simplify){


    res1=simplifyposteriors.infer(X = runemjmcmc.params$data,posteriors = posteriors, thf = report.level,resp = as.character(runemjmcmc.params$formula)[2])
    rownames(res1) = NULL
    res1$feature = as.character(res1$feature)
  }
  posteriors=posteriors[order(posteriors$x, decreasing = T),]
  colnames(posteriors)=c("feature","posterior")
  rm(params)
  gc()
  return(list(feat.stat = cbind(res1$feature,res1$posterior),predictions = pred,allposteriors = posteriors, threads.stats = results))

}
