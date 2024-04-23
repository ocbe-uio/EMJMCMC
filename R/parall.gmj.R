#' @title A function to run parallel chains of (R)(G)MJMCMC algorithms
#' @param X a vector of lists of parameters of runemjmcmc as well as
#' several additional fields that must come after runemjmcmc parameters
#' such as:
#' \describe{
#'  \item{vect$simlen}{the number of parameters of runemjmcmc in vect}
#'  \item{vect$cpu}{the cpu id for to set the unique seed}
#'  \item{vect$NM}{the number of unique best models from runemjmcmc to
#'    base the output report upon}
#' }
#' @param M a number of cpus to be used (can only be equal to 1 on
#' Windows OS currently, up to a maximal number of cores can be used on
#' linux based systems)
#' @param preschedule if pseudoscheduling should be used for the jobs if
#' their number exceeds M (if TRUE) otherwise the jobs are performed
#' sequentially w.r.t. their order
#' @return a vector of lists of
#' \describe{
#'  \item{post.populi}{the total mass (sum of the marginal likelihoods times
#'    the priors of the visited models) from the addressed run of runemjmcmc}
#'  \item{p.post}{posterior probabilities of the covariates approximated by the
#'    addressed run of runemjmcmc}
#'  \item{cterm}{the best value of marginal likelihood times the prior from
#'    the addressed run of runemjmcmc}
#'  \item{fparam}{the final set of covariates returned by the addressed
#'    run of runemjmcmc}
#' }
#' @example /inst/examples/parall.gmj_example.R
#' @seealso runemjmcmc parall.gmj
#' @keywords methods models
#' @importFrom stringi stri_locate_all
#' @export
parall.gmj <- function(X, M = 16, preschedule = FALSE) {
  parallel::mclapply(
    X              = X,
    FUN            = do.call.emjmcmc,
    mc.preschedule = preschedule,
    mc.cores       = M,
    mc.cleanup     = TRUE
  )
}

#' @title A help function used by parall.gmj to run parallel chains of (R)(G)MJMCMC algorithms
#' @param vect a vector of parameters of runemjmcmc as well as several
#' additional fields that must come after runemjmcmc parameters such as:
#' \describe{
#'  \item{vect$simlen}{the number of parameters of runemjmcmc in vect}
#'  \item{vect$cpu}{the cpu id for to set the unique seed}
#'  \item{vect$NM}{the number of unique best models from runemjmcmc to
#'    base the output report upon}
#' }
#' @return a list of
#' \describe{
#'  \item{post.populi}{the total mass (sum of the marginal likelihoods times
#'    the priors of the visited models) from the addressed run of runemjmcmc}
#'  \item{p.post}{posterior probabilities of the covariates approximated by the
#'    addressed run of runemjmcmc}
#'  \item{cterm}{the best value of marginal likelihood times the prior from
#'    the addressed run of runemjmcmc}
#'  \item{fparam}{the final set of covariates returned by the addressed
#'    run of runemjmcmc}
#' }
#' @seealso runemjmcmc, parall.gmj
#' @keywords  methods models
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
rm(vals)
return(list(post.populi = post.populi, p.post =  ppp$p.post, cterm = cterm, fparam = fparam))
}
