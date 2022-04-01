EMJMCMC2016$methods(
  # local greedy optimization
  learnlocalND = function(model) {

    # Step.nd
    varcand <- model$varcur
    varglob <- model$varcur
    varcurb <- model$varcur
    mlikcand <- model$mlikcur
    waiccand <- model$waiccur
    modglob <- NULL
    fm <- NULL
    fmb <- NULL
    opt.achieved <- FALSE


    # estimate large jump in a reverse move
    if (model$reverse || is.infinite(mlikcand)) {
      vectbg <- buildmodel(max.cpu = 1, varcur.old = varcand, statid = model$statid, switch.type = 8, min.N = min.N, max.N = max.N)
      if (!is.null(vectbg[[1]]$formula)) {
        bgmod <- lapply(X = vectbg, FUN = .self$fitmodel)
        waiccand <- bgmod[[1]]$waic
        mlikcand <- bgmod[[1]]$mlik
      } else if (exists("statistics1")) {
        iidd <- bittodec(varcand) + 1
        waiccand <- statistics1[iidd, 2]
        mlikcand <- statistics1[iidd, 1]
      } else if (exists("hashStat")) {
        iidd <- paste(varcand, collapse = "")
        waiccand <- hash::values(hashStat[iidd])[2]
        mlikcand <- hash::values(hashStat[iidd])[1]
      }
    }


    if (printable.opt) print(paste("Begin with ", mlikcand))

    mlikglob <- mlikcand
    mlikcand <- mlikcand
    waiccand <- waiccand
    waicglob <- waiccand
    waiccur <- waiccand

    buf.M.nd <- M.nd
    if (model$switch.type == 5) {
      buf.M.nd <- Nvars - sum(varcurb)
    } else if (model$switch.type == 6) {
      buf.M.nd <- sum(varcurb)
    }
    if (buf.M.nd == 0) {
      buf.M.nd <- 1
    }
    if (M.nd < buf.M.nd) {
      buf.M.nd <- M.nd
    }

    for (iterat in 1:buf.M.nd)
    {
      withRestarts(tryCatch({
        # statistics <- bigmemory::describe(statistics)
        mmax.cpu <- max.cpu
        if (model$switch.type == 5) {
          mmax.cpu <- Nvars - sum(varcurb)
        } else if (model$switch.type == 6) {
          mmax.cpu <- sum(varcurb)
        }
        if (mmax.cpu == 0) {
          mmax.cpu <- 1
        }
        vect <- buildmodel(max.cpu = mmax.cpu, varcur.old = varcurb, statid = model$statid, switch.type = model$switch.type, min.N = min.N, max.N = max.N)

        cluster <- TRUE

        flag1 <- 0

        for (mod_id in 1:mmax.cpu)
        {
          if (is.null(vect[[mod_id]]$formula)) {
            flag1 <- flag1 + 1
          }
        }

        if (flag1 == mmax.cpu) {
          cluster <- FALSE
          if (printable.opt) print("!!!!Greedy Models already estimated!!!!")
        } else {
          res.par <- parallelize(X = vect, FUN = .self$fitmodel)
        }

        for (mod_id in 1:mmax.cpu)
        {
          varcand1 <- vect[[mod_id]]$varcur
          if (cluster) {
            waiccand1 <- res.par[[mod_id]]$waic
            mlikcand1 <- res.par[[mod_id]]$mlik
          } else if (exists("statistics1")) {
            iidd <- bittodec(varcand1) + 1
            waiccand1 <- statistics1[iidd, 2]
            mlikcand1 <- statistics1[iidd, 1]
          } else if (exists("hashStat")) {
            iidd <- paste(varcand1, collapse = "")
            waiccand1 <- hash::values(hashStat[iidd])[2]
            mlikcand1 <- hash::values(hashStat[iidd])[1]
          }

          if (objective == 0) {
            objcand <- waiccand1
            objcur <- waiccand
            objglob <- waicglob
          } else {
            objcand <- -mlikcand1
            objcur <- -mlikcand
            objglob <- -mlikglob
          }

          if (objcand < objcur || mod_id == 1) {
            varcand <- varcand1
            waiccand <- waiccand1
            mlikcand <- mlikcand1
            if (printable.opt) print(paste("GREEDY update local optima with ", objcand))
            if (cluster) {
              fm <- res.par[[mod_id]]
            }
          }
        }
        # if(printable.opt)print(waiccand)

        if (objective == 0) {
          objcand <- waiccand1
          objcur <- waiccand
          objglob <- waicglob
        } else {
          objcand <- -mlikcand1
          objcur <- -mlikcand
          objglob <- -mlikglob
        }

        if (objcur < objglob) {
          waicglob <- waiccand
          varglob <- varcand
          varcurb <- varcand
          mlikglob <- mlikcand
          if (cluster) {
            modglob <- fm
          }

          if (printable.opt) print(paste("GREEDY update global optima with ", objcur))
        }
      }), abort = function() {
        opt.achieved <- TRUE
        fm <- fmb
        closeAllConnections()
        options(error = traceback)
        onerr <- TRUE
      })

      if (objcur != objglob) {
        if (locstop.nd) {
          break
        } else {
          varcurb <- varcand
        }
      }
    }

    # !#if(printable.opt)print(points)

    if (model$reverse == FALSE) {
      vect <- buildmodel(max.cpu = 1, varcur.old = varcand, statid = model$statid, switch.type = type.randomize, min.N = min.N.randomize, max.N = max.N.randomize)

      varcur <- vect[[1]]$varcur
      # if(printable.opt)print(varcur)

      cluster <- TRUE



      if (is.null(vect[[1]]$formula)) {
        cluster <- FALSE
        if (printable.opt) print("!!!!Greedy reverse model already estimated!!!!")
      } else {
        mod <- lapply(X = vect, FUN = fitmodel)
      }

      if (cluster) {
        waiccur <- mod[[1]]$waic
        mlikcur <- mod[[1]]$mlik
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
      model.probs <- calculate.move.logprobabilities(switch.type = type.randomize, varold = varcand, varnew = model$varold, min.N = min.N.randomize, max.N = max.N.randomize)
      model.prob <- model.probs$log.switch.forw.prob
      model.prob.fix <- model.probs$log.switch.back.prob
      varcur <- varcand
      waiccur <- waiccand
      mlikcur <- mlikcand
    }

    return(list(varcur = varcur, waiccur = waiccur, mlikcur = mlikcur, log.prob.cur = model.prob, log.prob.fix = model.prob.fix, varglob = varglob, waicglob = waicglob, mlikglob = mlikglob, modglob = modglob))
  }
)
