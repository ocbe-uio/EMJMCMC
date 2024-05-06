#' @title sigmoid activation function
#' @param x a real number
#' @return sigmoid value
#' @examples sigmoid(10)
#' @keywords methods models
#' @export
sigmoid<- function(x) {
 return(1/(1+(exp(-x))))
}
