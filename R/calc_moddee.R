calc_moddee <- function(mliks) {
  if (is.null(unlist(mliks))) warning("mliks is NULL")
  which(unlist(mliks) == max(unlist(mliks), na.rm = TRUE))[1]
}
