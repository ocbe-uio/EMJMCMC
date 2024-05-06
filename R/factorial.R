#' @title Truncated factorial to avoid stack overflow for huge values
#' @description truncated factorial to avoid stack overflow for huge values
#' @param x a non-negative integer number
#' @return \code{truncfactorial(x)}, truncated factorial as min(x!,171!)
#' @examples truncfactorial(10)
#' @export
#' @keywords methods models
truncfactorial <- function(x) ifelse(x <= 170, gamma(abs(x) + 1), gamma(171))
