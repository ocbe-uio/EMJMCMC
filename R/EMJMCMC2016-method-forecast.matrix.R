EMJMCMC2016$methods(
  forecast.matrix = function(covariates, link.g) {
    res <- NULL
    ncases <- dim(covariates)[1]
    formula.cur <- stats::as.formula(paste(fparam[1], "/2 ~", paste0(fparam, collapse = "+")))
    nvars <- Nvars
    if (save.beta) {
      if (exists("statistics1")) {
        ids <- which(!is.na(statistics1[, 15]))
        lids <- length(ids)
        statistics1[ids, (17:(nvars + 16))][which(is.na(statistics1[ids, (17:(nvars + 16))]))] <- 0
        res <- t(statistics1[ids, 15]) %*% g(statistics1[ids, (16:(nvars + 16))] %*% t(stats::model.matrix(object = formula.cur, data = covariates)))
      } else if (exists("hashStat")) {
        if (allow_offsprings == 0) {
          if (fparam[1] == "Const") {
            linx <- Nvars + 3
          } else {
            linx <- Nvars + 1 + 3
          }
        } else {
          if (fparam[1] == "Const") {
            linx <- Nvars.max + 3
          } else {
            linx <- Nvars.max + 1 + 3
          }
        }

        lHash <- length(hashStat)
        mliks <- hash::values(hashStat)[which((1:(lHash * linx)) %% linx == 1)]
        betas <- hash::values(hashStat)[which((1:(lHash * linx)) %% linx == 4)]
        for (i in 1:(Nvars - 1))
        {
          betas <- cbind(betas, hash::values(hashStat)[which((1:(lHash * linx)) %% linx == (4 + i))])
        }
        betas <- cbind(betas, hash::values(hashStat)[which((1:(lHash * linx)) %% linx == (0))])
        betas[which(is.na(betas))] <- 0
        xyz <- which(unlist(mliks) != -10000)
        g.results[4, 2] <<- lHash
        moddee <- which(unlist(mliks) == max(unlist(mliks), na.rm = TRUE))[1]
        zyx <- array(data = NA, dim = lHash)
        nconsum <- sum(exp(-mliks[moddee] + mliks[xyz]), na.rm = TRUE)

        if (nconsum > 0) {
          zyx[xyz] <- exp(mliks[xyz] - mliks[moddee]) / nconsum
        } else {
          diff <- 0 - mliks[moddee]
          mliks <- mliks + diff
          nconsum <- sum(exp(-mliks[moddee] + mliks[xyz]), na.rm = TRUE)
          zyx[xyz] <- exp(mliks[xyz] - mliks[moddee]) / nconsum
        }


        res <- t(zyx) %*% g(betas %*% t(stats::model.matrix(object = formula.cur, data = covariates)))
      }
    } else {
      warning("No betas were saved. Prediction is impossible. Please change the search parameters and run the search again.")
    }

    return(list(forecast = res))
  }
)
