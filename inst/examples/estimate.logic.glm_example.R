X1 <- as.data.frame(
  array(data = rbinom(n = 50 * 1000, size = 1, prob = 0.3), dim = c(1000, 50))
)
Y1 <- -0.7 + 1 * ((1 - X1$V1) * (X1$V4)) + 1 * (X1$V8 * X1$V11) + 1 * (X1$V5 * X1$V9)
X1$Y1 <- round(1.0 / (1.0 + exp(-Y1)))

formula1 <- as.formula(
  paste(colnames(X1)[51], "~ 1 +", paste0(colnames(X1)[-c(51)], collapse = "+"))
)

estimate.logic.glm(
  formula = formula1, data = X1, family = binomial(), n = 1000, m = 50
)
