calc_nconsum <- function(mliks, moddee, xyz) {
  sum(exp(-mliks[moddee] + mliks[xyz]), na.rm = TRUE)
}
