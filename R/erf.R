#' @title erf activation function
#' @param x a real number
#' @return \code{erf(x)}, erf value
#' @examples erf(10)
#' @keywords methods models
#' @export
erf <- function(x) {
  return(2 * stats::pnorm(x * sqrt(2)) - 1)
}
