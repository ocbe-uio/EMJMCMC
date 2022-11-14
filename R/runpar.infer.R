runpar.infer=function(vect)
{
#define a function performing the map step for a given thread
  ret = NULL
  tryCatch({
    set.seed(as.integer(vect$cpu))
    do.call(runemjmcmc, vect[1:(length(vect)-5)])
    vals=hash::values(hashStat)
    fparam=mySearch$fparam
    cterm=max(vals[1,],na.rm = T)
    ppp=mySearch$post_proceed_results_hash(hashStat = hashStat)
    post.populi=sum(exp(hash::values(hashStat)[1,][1:vect$NM]-cterm),na.rm = T)

    betas = NULL
    mliks = NULL

    if(vect$save.beta){
      #get the modes of beta coefficients for the explored models
      Nvars=mySearch$Nvars
      linx =mySearch$Nvars+4
      lHash=length(hashStat)
      mliks = hash::values(hashStat)[which((1:(lHash * linx)) %% linx == 1)]
      betas = hash::values(hashStat)[which((1:(lHash * linx)) %% linx == 4)]
      for(i in 1:(Nvars-1))
      {
        betas=cbind(betas,hash::values(hashStat)[which((1:(lHash * linx)) %% linx == (4+i))])
      }
      betas=cbind(betas,hash::values(hashStat)[which((1:(lHash * linx)) %% linx == (0))])
    }

    preds = NULL
    if(vect$predict)
    {
      preds=mySearch$forecast.matrix.na.fast(link.g = vect$link, covariates = (vect$test),betas = betas,mliks.in = mliks)$forecast
    }


    ret = list(post.populi = post.populi, p.post =  ppp$p.post, cterm = cterm, preds = preds, fparam = fparam, betas = betas, mliks = mliks )
    if(length(cterm)==0){
      vect$cpu=as.integer(vect$cpu)+as.integer(stats::runif(1,1,10000))
      if(vect$cpu<50000)
        ret = runpar.infer(vect)
      else
        ret = NULL
    }

  },error = function(err){
    print(paste0("error in thread",  vect[length(vect)]))
    print(err)
    vect$cpu=as.integer(vect$cpu)+as.integer(stats::runif(1,1,10000))
    if(vect$cpu<50000)
      ret = runpar.infer(vect)
    else
      ret =err
  },finally = {
    return(ret)
  })
}
