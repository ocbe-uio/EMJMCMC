EMJMCMC2016$methods(
  post_proceed_results_hash = function(hashStat) {
    if (save.beta) {
      if (allow_offsprings == 0 || Nvars > Nvars.max) {
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
    } else {
      linx <- 3
    }

    lHash <- length(hashStat)
    mliks <- hash::values(hashStat)[which((1:(lHash * linx)) %% linx == 1)]
    xyz <- which(unlist(mliks) != -10000)
    g.results[4, 2] <<- lHash
    moddee <- calc_moddee(mliks)
    zyx <- array(data = NA, dim = lHash)
    nconsum <- calc_nconsum(mliks, moddee, xyz)

    if (nconsum > 0) {
      zyx[xyz] <- exp(mliks[xyz] - mliks[moddee]) / nconsum
    } else {
      diff <- 0 - mliks[moddee]
      mliks <- mliks + diff
      nconsum <- calc_nconsum(mliks, moddee, xyz)
      zyx[xyz] <- exp(mliks[xyz] - mliks[moddee]) / nconsum
    }


    keysarr <- as.array(hash::keys(hashStat))
    p.post <- array(data = 0, dim = Nvars)
    for (i in 1:lHash)
    {
      if (is.na(zyx[i])) {
        del(x = keysarr[i], hash = hashStat)
        next
      }
      # vec<-dectobit(strtoi(keysarr[i], base = 0L)-1) # we will have to write some function that would handle laaaaargeee integers!
      varcur <- as.integer(strsplit(keysarr[i], split = "")[[1]])
      if (length(which(is.na(varcur))) > 0) {
        del(x = keysarr[i], hash = hashStat)
        next
      }
      if (length(varcur) > Nvars) {
        varcur <- varcur[1:Nvars]
      }
      p.post <- (p.post + varcur * zyx[i])
    }

    if (!exists("p.post") || is.null(p.post) || sum(p.post, na.rm = T) == 0 || sum(p.post, na.rm = T) > Nvars) {
      p.post <- array(data = 0.5, dim = Nvars)
    }
    # hash::values(hashStat)[which((1:(lHash * linx)) %%linx == 4)]<-zyx

    return(list(p.post = p.post, m.post = zyx, s.mass = sum(exp(mliks), na.rm = TRUE)))
  }
)
