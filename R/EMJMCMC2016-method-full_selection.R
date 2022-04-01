
EMJMCMC2016$methods(
  # full selection procedure
  full_selection = function(model) {
    if (printable.opt) print(paste("begin full selection procedure!", "Careful, ", 2^Nvars, " models have to be estimated"))
    if (Nvars > 30) {
      if (printable.opt) print("Finishing the procedure might well take forever!")
    }

    varcand <- array(0, Nvars)
    varcurb <- varcand
    varglob <- varcand
    varglob <- NULL
    modglob <- NULL
    mlikglob <- model$mlikcur
    mlikcur <- model$mlikcur
    waiccand <- model$waiccur
    waicglob <- model$waiccur
    waiccur <- model$waiccur
    waiccurb <- model$waiccur


    fm <- NULL
    fmb <- NULL

    ubs <- as.integer(bittodec(array(1, Nvars)) + 1)

    ub <- model$ub

    totit <- as.integer(ubs / ub) + 1

    if (model$totalit < totit) {
      totit <- model$totalit
    }

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
      if (!is.na(waiccur) && !is.na(mlikcur)) {
        mlikglob <- mlikcur
        mlikcand <- mlikcur
        waiccand <- waiccur
        waicglob <- waiccur
        waiccur <- waiccur
        waiccurb <- waiccur
      }
    }


    for (i in 1:totit)
    {
      if (ub * i > ubs) {
        ub <- ubs - ub * (i - 1) - 1
        if (printable.opt) print(paste("last ", ub, " iterations to complete"))
        varcurb <- varcand
      }
      withRestarts(tryCatch({
        vect <- buildmodel(max.cpu = ub, varcur.old = varcurb, statid = model$statid, switch.type = 7, shift.cpu = model$ub * (i - 1), min.N = min.N, max.N = max.N)
        if (printable.opt) print(paste("proceed with full ecumeration"))
        if (printable.opt) print(paste("current solution is", as.character(varcand)))
        cluster <- TRUE
        flag1 <- 0
        for (mod_id in 1:ub)
        {
          if (is.null(vect[[mod_id]]$formula)) {
            flag1 <- flag1 + 1
          }
        }
        if (flag1 == ub) {
          cluster <- FALSE
          if (printable.opt) print("!!!!full models already estimated!!!!")
        } else {
          res.par <- parallelize(X = vect, FUN = .self$fitmodel)
        }

        if (printable.opt) print(paste("end optimizing full ecumeration"))

        for (mod_id in 1:ub)
        {
          if (cluster) {
            fmb <- fm
            fm <- res.par[[mod_id]]
            waiccand <- Inf
            if (is.null(fm) && (is.na(res.par[[mod_id]]$waic))) {
              varcand <- varcurb
              if (printable.opt) print("full Model Fit Error!?")
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


          if (objcand <= objcur) {
            if (printable.opt) print(paste("full accept with ", objcand))
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

              if (printable.opt) print(paste("full global optima with ", objcand))
            }
          }
          varcurb <- varcand
        }

        # waiccurb<-waiccur
      }), abort = function() {
        varcur <- varcurb
        fm <- fmb
        closeAllConnections()
        options(error = traceback)
        onerr <- TRUE
      })
    }


    #                                                     if(length(which(varcur == 1))==(Nvars-1))
    #                                                     {
    #                                                       vectbg<-buildmodel(max.cpu = 1,varcur.old = varcurb,statid = model$statid,switch.type=8, min.N = min.N,max.N = max.N)
    #                                                       if(!is.null(vectbg[[1]]$formula))
    #                                                       {
    #                                                         bgmod <- lapply(X = vectbg,FUN = .self$fitmodel)
    #                                                         waiccand<-bgmod[[1]]$waic
    #                                                         mlikcand<-bgmod[[1]]$mlik
    #                                                       }
    #                                                       else
    #                                                       {
    #                                                         iidd<-bittodec(varcur)+1
    #                                                         waiccand<-statistics1[iidd,2]
    #                                                         mlikcand<-statistics1[iidd,1]
    #                                                       }
    #                                                       if(objective==0)
    #                                                       {
    #                                                         objcand<-waiccand
    #                                                         objcur<-waiccur
    #                                                         objglob<-waicglob
    #                                                       }else
    #                                                       {
    #                                                         objcand<- -mlikcand
    #                                                         objcur<-  -mlikcur
    #                                                         objglob<- -mlikglob
    #                                                       }
    #
    #
    #                                                       if(objcand<=objcur)
    #                                                       {
    #                                                         if(printable.opt)print(paste("full accept with ", objcand))
    #                                                         objcur<-objcand
    #                                                         varcurb<-varcand
    #                                                         waiccur<-waiccand
    #                                                         varcur<-varcand
    #                                                         mlikcur<-mlikcand
    #
    #                                                         if(objcur<objglob)
    #                                                         {
    #                                                           objglob<-objcur
    #                                                           waicglob<-waiccand
    #                                                           varglob<-varcand
    #                                                           mlikglob<-mlikcand
    #                                                           if(!is.null(fm))
    #                                                             modglob<-fm
    #
    #                                                           if(printable.opt)print(paste("full global optima with ", objcand))
    #                                                         }
    #                                                       }
    #                                                     }
    #
    model.prob <- 1


    model.prob.fix <- 1


    return(list(varcur = varglob, waiccur = waicglob, mlikcur = mlikglob, log.prob.cur = model.prob, log.prob.fix = model.prob.fix, varglob = varglob, waicglob = waicglob, mlikglob = mlikglob, modglob = modglob))
  }
)
