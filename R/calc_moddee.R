calc_moddee <- function(mliks) {
  if (length(mliks) == 0) {
    return(0)
  }
  which(unlist(mliks) == max(unlist(mliks), na.rm = TRUE))[1]
}
