EMJMCMC2016$methods(
  forecast.matrix.na.fast = function(covariates, link.g, betas, mliks.in) {
    formula.cur <- stats::as.formula(paste(fparam[1], "/2 ~", paste0(fparam, collapse = "+")))
    na.ids <- matrix(data = 0, nrow = dim(covariates)[1], ncol = dim(covariates)[2])
    na.ids[which(is.na(covariates))] <- 1
    na.bc <- colSums(na.ids)
    ids.betas <- NULL
    res.na <- array(NA, dim(covariates)[1])
    for (i in which(na.bc > 0))
    {
      if (((names(covariates)[i] %in% fparam) || (paste0("I(", names(covariates)[i], ")", collapse = "") %in% fparam))) {
        ids.betas <- c(ids.betas, which(fparam == names(covariates)[i] | fparam == paste0("I(", names(covariates)[i], ")", collapse = "")))
        next
      } else {
        na.ids[which(is.na(covariates[, i])), i] <- 0
        na.bc[i] <- 0
      }
    }
    na.br <- rowSums(na.ids)
    lv.br <- sort(unique(na.br))
    if (lv.br[1] == 0) {
      ids <- which(na.br == 0)
      mliks <- mliks.in
      xyz <- which(unlist(mliks) != -10000)
      moddee <- which(unlist(mliks) == max(unlist(mliks), na.rm = TRUE))[1]
      zyx <- array(data = NA, dim = lHash)
      nconsum <- sum(exp(-mliks[moddee] + mliks[xyz]), na.rm = TRUE)
      betas1 <- betas
      betas1[which(is.na(betas1))] <- 0
      if (nconsum > 0) {
        zyx[xyz] <- exp(mliks[xyz] - mliks[moddee]) / nconsum
      } else {
        diff <- 0 - mliks[moddee]
        mliks <- mliks + diff
        nconsum <- sum(exp(-mliks[moddee] + mliks[xyz]), na.rm = TRUE)
        zyx[xyz] <- exp(mliks[xyz] - mliks[moddee]) / nconsum
      }

      covariates1 <- covariates[ids, ]
      res <- t(zyx) %*% g(betas1 %*% t(stats::model.matrix(object = formula.cur, data = covariates1)))
      res.na[ids] <- res
      rm(mliks)
      rm(res)
      rm(zyx)
      rm(xyz)
      rm(covariates1)
      rm(betas1)
    }
    ids <- (which(is.na(res.na)))
    for (iii in which(na.bc > 0))
    {
      if (sum(na.ids[ids, iii]) > 0) {
        var.del <- (which(fparam == names(covariates)[iii] | fparam == paste0("I(", names(covariates)[iii], ")", collapse = "")))
        if (length(var.del > 0)) {
          w.ids <- union(w.ids, which(!is.na(betas[, (var.del + 1)])))
        }
      }
    }
    mliks <- mliks.in[-w.ids]
    if (length(mliks) == 0) {
      warning("not enough models for bagging in prediction. please train the model longer!")
      return(-1)
    }
    xyz <- which(unlist(mliks) != -10000)
    moddee <- which(unlist(mliks) == max(unlist(mliks), na.rm = TRUE))[1]
    zyx <- array(data = NA, dim = length(mliks))
    nconsum <- sum(exp(-mliks[moddee] + mliks[xyz]), na.rm = TRUE)
    betas1 <- betas[-w.ids, ]
    betas1[which(is.na(betas1))] <- 0
    if (nconsum > 0) {
      zyx[xyz] <- exp(mliks[xyz] - mliks[moddee]) / nconsum
    } else {
      diff <- 0 - mliks[moddee]
      mliks <- mliks + diff
      nconsum <- sum(exp(-mliks[moddee] + mliks[xyz]), na.rm = TRUE)
      zyx[xyz] <- exp(mliks[xyz] - mliks[moddee]) / nconsum
    }
    covariates1 <- as.matrix(covariates[ids, ])
    covariates1[which(is.na(covariates1))] <- 0
    covariates1 <- as.data.frame(covariates1)
    res <- t(zyx) %*% g(betas1 %*% t(stats::model.matrix(object = formula.cur, data = covariates1)))
    res.na[ids] <- res
    rm(mliks)
    rm(zyx)
    rm(xyz)
    rm(res)
    rm(covariates1)


    return(list(forecast = res.na))
  }
)
