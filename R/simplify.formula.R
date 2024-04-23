#' @title A function parsing the formula into the vectors of character arrays
#' of responses and covariates
#' @param fmla an R formula object
#' @param names all column names from the data.frame to be used with the formula
#' @return a list of
#' \describe{
#'  \item{fobserved}{a vector of character arrays corresponding to the observations}
#'  \item{fparam}{a vector of character arrays corresponding to the covariates}
#' }
#' @seealso formula data.frame
#' @example inst/examples/simplify.formula_example.R
#' @keywords methods models
#' @export
simplify.formula<-function(fmla,names)
{
fmla.proc<-as.character(fmla)[2:3]
fobserved <- fmla.proc[1]
fmla.proc[2]<-stringi::stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
fmla.proc[2]<-stringi::stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
fparam <- names[which(names %in% stringi::stri_split_fixed(str = fmla.proc[2],pattern = "+",omit_empty = F)[[1]] )]
return(list(fparam = fparam,fobserved = fobserved))
}
