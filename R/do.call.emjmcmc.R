do.call.emjmcmc<-function(vect)
{

set.seed(as.integer(vect$cpu))
do.call(runemjmcmc, vect[1:vect$simlen])
vals<-hash::values(hashStat)
fparam<-mySearch$fparam
cterm<-max(vals[1,],na.rm = T)
ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
post.populi<-sum(exp(hash::values(hashStat)[1,][1:vect$NM]-cterm),na.rm = T)
hash::clear(hashStat)
#rm(hashStat)
rm(vals)
gc()
return(list(post.populi = post.populi, p.post =  ppp$p.post, cterm = cterm, fparam = fparam))
}
