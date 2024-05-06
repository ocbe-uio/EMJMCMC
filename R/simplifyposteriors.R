#' @title A function that ads up posteriors for the same expression written
#' in different character form in different parallel runs of the algorithm
#' (mainly for Logic Regression and Deep Regression contexts)
#' @param X a data.frame containing the data on the covariates
#' @param posteriors a data.frame with expressions in the first column and their posteriors in the second column from all of the runs
#' @param th initial filtering before summary threshold
#' @param thf threshold for final filtering after summary
#' @param resp the response to be addressed
#' @return res, a data.frame with the summarized across runs expressions and
#' their posteriors
#' @seealso runemjmcmc
#' @keywords methods models
simplifyposteriors<-function(X,posteriors,th=0.0001,thf=0.2, resp)
{
posteriors<-posteriors[which(posteriors[, 2] >= th), ]
rhash <- hash::hash()
for(i in 1:length(posteriors[,1])) {
  expr<-posteriors[i,1]
  res<-stats::model.matrix(data=X,object = stats::as.formula(paste0(resp,"~",expr)))
  ress <- c(
    stringi::stri_flatten(round(sum(res[, 2]), digits = 4), collapse = ""),
    stringi::stri_flatten(res[, 2], collapse = ""),
    posteriors[i, 2],
    expr
  )
  if(!((ress[1] %in% hash::values(rhash)))) {
    rhash[[ress[1]]]<-ress
  } else {
    if(ress[1] %in% hash::keys(rhash)) {
      rhash[[ress[1]]][3]<- (as.numeric(rhash[[ress[1]]][3]) + as.numeric(ress[3]))
      if(stringi::stri_length(rhash[[ress[1]]][4])>stringi::stri_length(expr)) {
        rhash[[ress[1]]][4]<-expr
      }
    }
  }

}
res<-as.data.frame(t(hash::values(rhash)[c(3,4),]))
res$V1<-as.numeric(as.character(res$V1))
res<-res[which(res$V1>thf),]
res<-res[order(res$V1, decreasing = TRUE),]
hash::clear(rhash)
rm(rhash)
res[which(res[,1]>1),1]<-1
colnames(res)<-c("posterior","tree")
return(res)
}
