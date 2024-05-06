EMJMCMC2016$methods(
  # local simulated annealing optimization
  learnlocalSA = function(model) {
    t.min <- SA.param$t.min
    t <- SA.param$t.init
    dt <- SA.param$dt
    M <- SA.param$M
    varcand <- model$varcur
    varcurb <- model$varcur
    varcur <- model$varcur
    varglob <- model$varcur
    varglob <- NULL
    modglob <- NULL
    probcur <- 1
    probrevcur <- 1
    first.prob <- 1
    fm <- NULL
    fmb <- NULL

    mlikcur <- model$mlikcur
    waiccur <- model$waiccur
    # estimate large jump in a reverse move
    # estimate large jump in a reverse move
    if ((model$reverse && !model$sa2) || is.infinite(mlikcur)) {
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

    while (t > t.min) {
      if (printable.opt) print(paste("anneal to ", t))
      t.new <- t * exp(-dt)
      if (model$reverse == TRUE && t.new <= t.min) {
        M <- M - 1
      }
      for (m in 1:M)
      {
        withRestarts(tryCatch({
          mmax.cpu <- max.cpu
          if (model$switch.type == 5) {
            mmax.cpu <- Nvars - sum(varcur)
          } else if (model$switch.type == 6) {
            mmax.cpu <- sum(varcur)
          }
          if (mmax.cpu == 0) {
            mmax.cpu <- 1
          }

          vect <- buildmodel(max.cpu = mmax.cpu, varcur.old = varcur, statid = model$statid, switch.type = model$switch.type, min.N = min.N, max.N = max.N)
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
            if (printable.opt) print("!!!!SA Models already estimated!!!!")
          } else {
            res.par <- parallelize(X = vect, FUN = .self$fitmodel)
          }
          for (mod_id in 1:mmax.cpu)
          {
            if (cluster) {
              fmb <- fm
              fm <- res.par[[mod_id]]
              waiccand <- Inf
              if (is.null(fm) && (is.na(res.par[[mod_id]]$waic))) {
                varcand <- varcurb
                if (printable.opt) print("SA Model Fit Error!?")
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
              iidd <- stri_paste(varcand, collapse = "")
              waiccand <- hash::values(hashStat[iidd])[2]
              mlikcand <- hash::values(hashStat[iidd])[1]
            }

            if (objective == 0) {
              objcand <- waiccand
              objcur <- waiccur
              objglob <- waicglob
            } else {
              objcand <- -mlikcand
              objcur <- -mlikcur
              objglob <- -mlikglob
            }

            if (t == SA.param$t.init && mod_id == 1 && m == 2) {
              delta <- objcand - objcur
              first.prob <- vect[[mod_id]]$log.mod.switchback.prob + log(punif(q = exp(x = delta / t), min = 0, max = 1))
            }

            if (objcand < objcur) {
              if (printable.opt) print(paste("SA accept move with ", objcand))
              waiccur <- waiccand
              varcur <- varcand
              mlikcur <- mlikcand
              if (objcand < objglob) {
                waicglob <- waiccand
                varglob <- varcand
                mlikglob <- mlikcand
                if (cluster) {
                  modglob <- fm
                }

                if (printable.opt) print(paste("SA update global optima with", objcand))
              }
            } else {
              delta <- objcand - objcur
              if (stats::runif(n = 1, min = 0, max = 1) <= exp(x = -delta / t)) {
                model.probs <- calculate.move.logprobabilities(varold = varcur, varnew = varcand, switch.type = model$switch.type, min.N = min.N, max.N = max.N)
                probcur <- model.probs$log.switch.forw.prob
                probrevcur <- model.probs$log.switch.back.prob
                waiccur <- waiccand
                varcur <- varcand
                mlikcur <- mlikcand
                if (printable.opt) print(paste("SA accept move with ", objcand))
              }
            }
          }
        }), abort = function() {
          varcur <- varcurb
          fm <- fmb
          closeAllConnections()
          withr::local_options(error = traceback)
          onerr <- TRUE
        })
      }
      t <- t.new
    }
    t <- t / exp(-dt)

    if (model$reverse == FALSE) {
      model.prob <- log(punif(q = exp(x = -delta / t), min = 0, max = 1)) + probcur # log(P(Mk,Mk-1))
      model.prob.fix <- log(punif(q = exp(x = delta / t), min = 0, max = 1)) + probrevcur # log(P(Mk-1,Mk))

      if (model$sa2 == TRUE) {
        model.prob.fix <- model.prob.fix + first.prob # correcting for the term for local improvements of type 3.
      }
    } else # incorporate what happens for the reverse move
    {
      if (is.null(varcur)) {
        if (printable.opt) print("No moves accepted in the reverse procedure")
        varcur <- model$varcur
        objcur <- model$objcur
      }

      delta <- objcur - model$objold


      model.probs <- calculate.move.logprobabilities(varold = varcur, varnew = model$varold, switch.type = model$switch.type, min.N = min.N, max.N = max.N)
      model.prob <- punif(q = exp(x = -delta / t), min = 0, max = 1, log.p = TRUE) + model.probs$log.switch.forw.prob
      model.prob.fix <- punif(q = exp(x = delta / t), min = 0, max = 1, log.p = TRUE) + model.probs$log.switch.back.prob

      if (model.prob == -Inf) {
        model.prob <- -100000000
      }
      if (model.prob.fix == -Inf) {
        model.prob.fix <- -100000000
      }
    }

    return(list(varcur = varcur, waiccur = waiccur, mlikcur = mlikcur, log.prob.cur = model.prob, log.prob.fix = model.prob.fix, varglob = varglob, waicglob = waicglob, mlikglob = mlikglob, modglob = modglob))
  }
)
