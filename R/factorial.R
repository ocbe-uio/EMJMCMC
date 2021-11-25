#' @title Truncated factorial to avoid stack overflow for huge values
#' @description truncated factorial to avoid stack overflow for huge values
#' @param x a non-negative integer number
#' @return \code{factorial(x)}, truncated facctorial as min(x!,171!)
#' @examples factorial(10)
#' @export
#' @keywords methods models
factorial<-function(x) ifelse(x<=170,gamma(abs(x) + 1),gamma(171))
