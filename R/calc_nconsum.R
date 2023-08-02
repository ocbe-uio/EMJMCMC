calc_nconsum <- function(mliks, moddee, xyz) {
  if (is.null(unlist(mliks))) warning("mliks is NULL")
  sum(exp(-mliks[moddee] + mliks[xyz]), na.rm = TRUE)
}
