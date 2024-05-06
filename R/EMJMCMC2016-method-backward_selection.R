EMJMCMC2016$methods(
  # backward selection procedure
  backward_selection = function(model) {
    # if(printable.opt)print("begin backward selection procedure")
    varcand <- model$varcur
    varcurb <- model$varcur
    varglob <- varcand
    mlikglob <- model$mlikcur
    mlikcur <- model$mlikcur
    waiccur <- model$waiccur
    waicglob <- model$waiccur
    varglob <- NULL
    modglob <- NULL
    waiccurb <- model$waiccur

    fm <- NULL
    fmb <- NULL

    ub <- bittodec(array(1, length(varcurb)))
    layer <- length(which(varcurb == 1))

    if (is.infinite(mlikcur)) {
      vectbg <- buildmodel(max.cpu = 1, varcur.old = varcurb, statid = model$statid, switch.type = 8, min.N = min.N, max.N = max.N)
      if (!is.null(vectbg[[1]]$formula)) {
        bgmod <- lapply(X = vectbg, FUN = .self$fitmodel)
        waiccur <- bgmod[[1]]$waic
        mlikcur <- bgmod[[1]]$mlik
      } else if (exists("statistics1")) {
        iidd <- bittodec(varcand) + 1
        waiccur <- statistics1[iidd, 2]
        mlikcur <- statistics1[iidd, 1]
      } else if (exists("hashStat")) {
        iidd <- paste(varcand, collapse = "")
        waiccur <- hash::values(hashStat[iidd])[2]
        mlikcur <- hash::values(hashStat[iidd])[1]
      }
      if (!is.na(mlikcur) && !is.na(waiccur)) {
        mlikglob <- mlikcur
        mlikcur <- mlikcur
        waiccand <- waiccur
        waicglob <- waiccur
        waiccur <- waiccur
        waiccurb <- waiccur
      }
    }


    while (layer > 0) {
      withRestarts(tryCatch({
        if (printable.opt) print(paste("backward proceed with layer", layer))
        if (printable.opt) print(paste("current backward solution is", as.character(varcand)))
        vect <- buildmodel(max.cpu = layer, varcur.old = varcurb, statid = model$statid, switch.type = 6, min.N = min.N, max.N = max.N)

        if (printable.opt) print(paste("finish backward preparing models at layer", layer))

        cluster <- TRUE
        flag1 <- 0
        for (mod_id in 1:layer)
        {
          if (is.null(vect[[mod_id]]$formula)) {
            flag1 <- flag1 + 1
          }
        }
        if (flag1 == layer) {
          cluster <- FALSE
          if (printable.opt) print("!!!!backward Models already estimated!!!!")
        } else {
          res.par <- parallelize(X = vect, FUN = .self$fitmodel)
        }
        if (printable.opt) print(paste("end backward optimizing at layer", layer))

        for (mod_id in 1:layer)
        {
          if (cluster) {
            fmb <- fm
            fm <- res.par[[mod_id]]
            waiccand <- Inf
            if (is.null(fm) && (is.na(res.par[[mod_id]]$waic))) {
              varcand <- varcurb
              if (printable.opt) print("backward Model Fit Error!?")
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

          if (objective == 0) {
            objcand <- waiccand
            objcur <- waiccur
            objglob <- waicglob
          } else {
            objcand <- -mlikcand
            objcur <- -mlikcur
            objglob <- -mlikglob
          }


          if (objcand <= objcur || mod_id == 1) {
            if (printable.opt) print(paste("backward accept with ", objcand))
            objcur <- objcand
            varcurb <- varcand
            waiccur <- waiccand
            varcur <- varcand
            mlikcur <- mlikcand

            if (objcur < objglob) {
              objglob <- objcur
              waicglob <- waiccand
              varglob <- varcand
              mlikglob <- mlikcand
              if (!is.null(fm)) {
                modglob <- fm
              }

              if (printable.opt) print(paste("backward global optima with ", objcand))
            }
          }
        }
        if (objcur != objglob) {
          if (model$locstop) {
            break
          } else {
            varcurb <- varcur
          }
        }
      }), abort = function() {
        varcur <- varcurb
        fm <- fmb
        closeAllConnections()
        withr::local_options(error = traceback)
        onerr <- TRUE
      })


      layer <- layer - 1
    }



    model.prob <- 1


    model.prob.fix <- 1


    return(list(varcur = varglob, waiccur = waicglob, mlikcur = mlikglob, log.prob.cur = model.prob, log.prob.fix = model.prob.fix, varglob = varglob, waicglob = waicglob, mlikglob = mlikglob, modglob = modglob))
  }
)
