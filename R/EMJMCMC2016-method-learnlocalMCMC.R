EMJMCMC2016$methods(
  # local mcmc procedure
  learnlocalMCMC = function(model) {
    M <- M.mcmc
    mlikcur <- model$mlikcur
    waiccur <- model$waiccur
    varcand <- model$varcur
    varcur <- model$varcur
    varcurb <- model$varcur
    varglob <- model$varcur

    modglob <- NULL
    fm <- NULL
    fmb <- NULL

    # estimate large jump in a reverse move
    if (model$reverse || is.infinite(mlikcur)) {
      vectbg <- buildmodel(max.cpu = 1, varcur.old = varcur, statid = model$statid, switch.type = 8, min.N = min.N, max.N = max.N)
      if (!is.null(vectbg[[1]]$formula)) {
        bgmod <- lapply(X = vectbg, FUN = .self$fitmodel)
        waiccur <- bgmod[[1]]$waic
        mlikcur <- bgmod[[1]]$mlik
      } else if (exists("statistics1")) {
        iidd <- bittodec(varcur) + 1
        waiccur <- statistics1[iidd, 2]
        mlikcur <- statistics1[iidd, 1]
      } else if (exists("hashStat")) {
        iidd <- paste(varcur, collapse = "")
        waiccur <- hash::values(hashStat[iidd])[2]
        mlikcur <- hash::values(hashStat[iidd])[1]
      }
    }


    if (printable.opt) print(paste("Begin with ", mlikcur))

    mlikglob <- mlikcur
    mlikcand <- mlikcur
    waiccand <- waiccur
    waicglob <- waiccur
    waiccur <- waiccur

    for (m in 1:M)
    {
      withRestarts(tryCatch({

        # statistics <- bigmemory::describe(statistics)
        vect <- buildmodel(varcur.old = varcur, statid = model$statid, max.cpu = max.cpu, switch.type = switch.type, min.N = min.N, max.N = max.N)
        cluster <- TRUE
        flag1 <- 0

        for (mod_id in 1:max.cpu)
        {
          if (is.null(vect[[mod_id]]$formula)) {
            flag1 <- flag1 + 1
          }
        }

        # flag1<-sum(is.null(vect[[]]$formula))

        if (flag1 == max.cpu) {
          cluster <- FALSE
          if (printable.opt) print("!!!!MTMCMC models already estimated!!!!")
        } else {
          res.par <- parallelize(X = vect, FUN = .self$fitmodel)
        }
        p.select.y <- array(data = 0, dim = max.cpu)
        for (mod_id in 1:max.cpu)
        {
          if (cluster) {
            fm <- res.par[[mod_id]]

            if (is.null(fm) && (is.na(res.par[[mod_id]]$waic))) {
              varcand <- varcurb
              if (printable.opt) print("locMTMCMC Model Fit Error!?")
              next
            }
          }

          varcand <- vect[[mod_id]]$varcur

          if (cluster) {
            waiccand <- res.par[[mod_id]]$waic
            mlikcand <- res.par[[mod_id]]$mlik
          } else if (exists("statistics1")) {
            iidd <- bittodec(varcand) + 1
            waiccand <- statistics1[iidd, 2]
            mlikcand <- statistics1[iidd, 1]
          } else if (exists("hashStat")) {
            iidd <- paste(varcand, collapse = "")
            waiccand <- hash::values(hashStat[iidd])[2]
            mlikcand <- hash::values(hashStat[iidd])[1]
          }

          if ((mlikcand > mlikglob)) # update the parameter of interest
            {
              if (printable.opt) print(paste("locMTMCMC update waic.glob = ", waiccand))
              if (printable.opt) print(paste("locMTMCMC update waic.glob.mlik = ", mlikglob))
              mlikglob <- mlikcand
              waicglob <- waiccand
              varglob <- varcand
              if (cluster) {
                modglob <- fm
              }
            }


          g1 <- waiccur

          if (waiccur == Inf) {
            g1 <- 1
          }

          p.select.y[mod_id] <- (mlikcand + vect[[mod_id]]$log.mod.switchback.prob + log(lambda(c = cc, alpha = aa, g1 = -g1, g2 = -waiccand, g.domain.pos = FALSE))) # correct for different criteria later

          if (is.na(p.select.y[mod_id])) {
            p.select.y[mod_id] <- 0
          }
          if (is.infinite(p.select.y[mod_id]) || p.select.y[mod_id] > 100000000) {
            # if(printable.opt)print(paste("very large log.w.y detected ",p.select.y[mod_id]))
            p.select.y[mod_id] <- 100000000
          }
        }

        max.p.select.y <- max(p.select.y)
        p.select.y <- p.select.y - max.p.select.y

        if (printable.opt) print(paste("max log.w.y is ", max.p.select.y, "normilized log.w.n.y is ", paste(p.select.y, collapse = ", ")))


        ID <- sample(x = max.cpu, size = 1, prob = exp(p.select.y))

        if (printable.opt) print(paste("cand ", ID, " selected"))

        varcand <- vect[[ID]]$varcur

        if (cluster) {
          waiccand <- res.par[[ID]]$waic
          mlikcand <- res.par[[ID]]$mlik
        } else if (exists("statistics1")) {
          iidd <- bittodec(varcand) + 1
          waiccand <- statistics1[iidd, 2]
          mlikcand <- statistics1[iidd, 1]
        } else if (exists("hashStat")) {
          iidd <- paste(varcand, collapse = "")
          waiccand <- hash::values(hashStat[iidd])[2]
          mlikcand <- hash::values(hashStat[iidd])[1]
        }

        # p.Q.cand<- p.select.y[ID]/sum(p.select.y)

        if (printable.opt) print("do reverse step")

        p.select.z <- array(data = 0.01, dim = max.cpu)


        if (max.cpu != 1) {
          vect1 <- buildmodel(max.cpu = max.cpu - 1, varcur.old = varcand, statid = model$statid, switch.type = switch.type, min.N = min.N, max.N = max.N)

          cluster <- TRUE

          flag1 <- 0

          for (mod_id in 1:(max.cpu - 1))
          {
            if (is.null(vect1[[mod_id]]$formula)) {
              flag1 <- flag1 + 1
            }
          }

          if (flag1 == (max.cpu - 1)) {
            cluster <- FALSE
            if (printable.opt) print("!!!!MTMCMC reverse models already estimated!!!!")
          } else {
            res.par.back <- parallelize(X = vect1, FUN = .self$fitmodel)
          }

          for (mod_id in 1:(max.cpu - 1))
          {
            if (cluster) {
              if (is.null(fm) && (is.na(res.par.back[[mod_id]]$waic))) {
                if (printable.opt) print("locMTMCMC Model Fit Error!?")
                next
              }
            }

            varcand.b <- vect1[[mod_id]]$varcur

            if (cluster) {
              waiccand.b <- res.par.back[[mod_id]]$waic
              mlikcand.b <- res.par.back[[mod_id]]$mlik
            } else if (exists("statistics1")) {
              iidd <- bittodec(varcand.b) + 1
              waiccand.b <- statistics1[iidd, 2]
              mlikcand.b <- statistics1[iidd, 1]
            } else if (exists("hashStat")) {
              iidd <- paste(varcand.b, collapse = "")
              waiccand.b <- hash::values(hashStat[iidd])[2]
              mlikcand.b <- hash::values(hashStat[iidd])[1]
            }

            if ((mlikcand.b > mlikglob)) {
              if (printable.opt) print(paste("locMTMCMC update waic.glob = ", waiccand.b))
              if (printable.opt) print(paste("locMTMCMC update waic.glob.mlik = ", mlikcand.b))
              mlikglob <- mlikcand.b
              waicglob <- waiccand.b
              varglob <- varcand.b
              if (cluster) {
                modglob <- fm
              }
            }

            g1 <- waiccand

            if (waiccand == Inf) {
              g1 <- 1
            }

            p.select.z[mod_id] <- (mlikcand.b + vect1[[mod_id]]$log.mod.switchback.prob + (lambda(c = cc, alpha = aa, g1 = -g1, g2 = -waiccand.b, g.domain.pos = FALSE))) # correct for different criteria later

            if (is.na(p.select.z[mod_id])) {
              p.select.z[mod_id] <- 0
            }
            if (is.infinite(p.select.z[mod_id]) || p.select.z[mod_id] > 100000000) {
              # if(printable.opt)print(paste("very large log.w.y detected ",p.select.z[mod_id]))
              p.select.z[mod_id] <- 100000000
            }
          }
        }

        if (waiccur == Inf) {
          g1 <- 1
        }
        p.select.z[max.cpu] <- (mlikcur + vect[[ID]]$log.mod.switch.prob + (lambda(c = cc, alpha = aa, g1 = -g1, g2 = -waiccand, g.domain.pos = FALSE)))

        if (is.na(p.select.z[mod_id])) {
          p.select.z[mod_id] <- 0
        }
        if (is.infinite(p.select.z[mod_id]) || p.select.z[mod_id] > 100000000) {
          # if(printable.opt)print(paste("very large log.w.y detected ",p.select.z[mod_id]))
          p.select.z[mod_id] <- 100000000
        }

        max.p.select.z <- max(p.select.z)
        p.select.z <- p.select.z - max.p.select.z

        if (printable.opt) print(paste("max log.w.z is ", max.p.select.z, "normilized log.w.n.z is ", paste(p.select.z, collapse = ", ")))

        if (log(stats::runif(n = 1, min = 0, max = 1)) < (log(sum(exp(p.select.y))) - log(sum(exp(p.select.z)))) + max.p.select.y - max.p.select.z) {
          mlikcur <- mlikcand
          if (printable.opt) print(paste("locMTMCMC update ratcur = ", mlikcand))
          if (printable.opt) print(paste("locMTMCMC accept move with ", waiccand))
          varcur <- varcand
          waiccur <- waiccand
        }
      }), abort = function() {
        fm <- fmb
        closeAllConnections()
        options(error = traceback)
        onerr <- TRUE
      })
    }

    # if(printable.opt)print("FINISH LOCAL MTMCMC")

    # !#if(printable.opt)print(points)

    if (model$reverse == FALSE) {
      if (is.null(varcur)) {
        # if(printable.opt)print("No moves acceoted in the procedure")
        varcur <- varcand
        waiccur <- waiccand
        varcur <- varcand
      }

      vect <- buildmodel(max.cpu = 1, varcur.old = varcur, statid = model$statid, switch.type = type.randomize, min.N = min.N.randomize, max.N = max.N.randomize)


      varcur <- vect[[1]]$varcur
      # if(printable.opt)print(varcur)

      cluster <- TRUE


      if (is.null(vect[[1]]$formula)) {
        cluster <- FALSE
        if (printable.opt) print("!!!!MTMCMC reverse models already estimated!!!!")
      } else {
        mod <- fitmodel(vect[[1]])
      }

      if (cluster) {
        waiccur <- mod$waic
        mlikcur <- mod$mlik
      } else if (exists("statistics1")) {
        iidd <- bittodec(varcur) + 1
        waiccur <- statistics1[iidd, 2]
        mlikcur <- statistics1[iidd, 1]
      } else if (exists("hashStat")) {
        iidd <- paste(varcur, collapse = "")
        waiccur <- hash::values(hashStat[iidd])[2]
        mlikcur <- hash::values(hashStat[iidd])[1]
      }
      # incorporate what happens for the backward optimization

      model.prob <- vect[[1]]$log.mod.switch.prob
      model.prob.fix <- vect[[1]]$log.mod.switchback.prob
    } else # incorporate what happens for the reverse move
    {
      if (is.null(varcur)) {
        # if(printable.opt)print("No moves acceoted in the reverse procedure")
        varcur <- varcand
        waiccur <- waiccand
      }
      model.probs <- calculate.move.logprobabilities(varold = varcur, varnew = model$varold, switch.type = type.randomize, min.N = min.N.randomize, max.N = max.N.randomize)
      model.prob <- model.probs$log.switch.forw.prob
      model.prob.fix <- model.probs$log.switch.back.prob
    }

    if (is.null(varcur)) {
      # if(printable.opt)print("NO VARCUR OBTAINED")
      varcur <- i
      waiccur <- Inf
      varcur <- model$varold
    }

    return(list(varcur = varcur, waiccur = waiccur, mlikcur = mlikcur, log.prob.cur = model.prob, log.prob.fix = model.prob.fix, varglob = varglob, waicglob = waicglob, mlikglob = mlikglob, modglob = modglob))
  }
)
