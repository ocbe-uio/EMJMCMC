EMJMCMC2016$methods(
  # fit the given model by means of the specified estimator
  fitmodel = function(model) {
    if (!is.null(model)) {
      fm <- NULL
      id <- bittodec(model$varcur)

      if (is.null(id)) {
        id <- 0
      }
      id <- id + 1

      if (exists("statistics1")) {
        if (is.na(statistics1[id, 1])) {
          onerr <- FALSE
          statistics1[id, c(2, 3)] <- 100000
          statistics1[id, 1] <- -100000
          statistics1[id, 4:14] <- 0
          utils::capture.output({
            withRestarts(tryCatch(utils::capture.output({
              fm <- do.call(estimator, c(estimator.args, model$formula))
            })), abort = function() {
              onerr <- TRUE
              fm <- NULL
            })
          }) # fit the model, get local improvements



          if (!is.null(fm)) {
            statistics1[id, 2] <- fm$waic[[1]]
            statistics1[id, 1] <- fm$mlik[[1]]
            statistics1[id, 3] <- fm$dic[[1]]
            if (save.beta) {
              if (fparam[1] == "Const") {
                inxx <- which(model$varcur == 1)
                if (length(inxx) == length(fm$summary.fixed$mean)) {
                  statistics1[id, 15 + inxx] <- fm$summary.fixed$mean
                }
              } else {
                inxx <- c(0, which(model$varcur == 1))
                if (length(inxx) == length(fm$summary.fixed$mean)) {
                  statistics1[id, 16 + inxx] <- fm$summary.fixed$mean
                }
              }
            }

            if (fm$waic[[1]] < g.results[2, 1] && !is.na(fm$waic[[1]])) {
              g.results[2, 1] <<- fm$waic[[1]]
              g.results[2, 2] <<- as.integer(id)
            }
            if (fm$mlik[[1]] > g.results[1, 1] && !is.na(fm$mlik[[1]])) {
              g.results[1, 1] <<- fm$mlik[[1]]
              g.results[1, 2] <<- as.integer(id)
            }

            if (fm$dic[[1]] < g.results[3, 1] && !is.na(fm$dic[[1]])) {
              g.results[3, 1] <<- fm$dic[[1]]
              g.results[3, 2] <<- as.integer(id)
            }

            g.results[4, 2] <<- g.results[4, 2] + 1
            if (g.results[4, 2] %% recalc.margin == 0) {
              p.add <<- as.array(post_proceed_results(statistics1)$p.post)
            }
          }
        }
        if (model$statid != -1) {
          statistics1[id, model$statid + 1] <- statistics1[id, model$statid + 1] + 1
        }
        g.results[4, 1] <<- g.results[4, 1] + 1
        return(list(mlik = statistics1[id, 1], waic = statistics1[id, 2], dic = statistics1[id, 3]))
      } else if (exists("statistics")) {
        if (is.na(statistics[id, 1])) {
          # if(printable.opt)print("Invoked from EMJMCMC SUB environment")
          onerr <- FALSE
          statistics[id, c(2, 3)] <- 100000
          statistics[id, 1] <- -100000
          statistics[id, 4:14] <- 0

          utils::capture.output({
            withRestarts(tryCatch(utils::capture.output({
              fm <- do.call(estimator, c(estimator.args, model$formula))
            })), abort = function() {
              onerr <- TRUE
              fm <- NULL
            })
          }) # fit the modal, get local improvements



          if (!is.null(fm)) {
            statistics[id, 2] <- fm$waic[[1]]
            statistics[id, 1] <- fm$mlik[[1]]
            statistics[id, 3] <- fm$dic[[1]]
            if (save.beta) {
              if (fparam[1] == "Const") {
                inxx <- which(model$varcur == 1)
                if (length(inxx) == length(fm$summary.fixed$mean)) {
                  statistics[id, 15 + inxx] <- fm$summary.fixed$mean
                }
              } else {
                inxx <- c(0, which(model$varcur == 1))
                if (length(inxx) == length(fm$summary.fixed$mean)) {
                  statistics[id, 16 + inxx] <- fm$summary.fixed$mean
                }
              }
            }

            if (fm$waic[[1]] < g.results[2, 1] && !is.na(fm$waic[[1]])) {
              g.results[2, 1] <<- fm$waic[[1]]
              g.results[2, 2] <<- as.integer(id)
            }
            if (fm$mlik[[1]] > g.results[1, 1] && !is.na(fm$mlik[[1]])) {
              g.results[1, 1] <<- fm$mlik[[1]]
              g.results[1, 2] <<- as.integer(id)
            }

            if (fm$dic[[1]] < g.results[3, 1] && !is.na(fm$dic[[1]])) {
              g.results[3, 1] <<- fm$dic[[1]]
              g.results[3, 2] <<- as.integer(id)
            }

            g.results[4, 2] <<- g.results[4, 2] + 1
            if (g.results[4, 2] %% recalc.margin == 0) {
              proceeeded <- post_proceed_results(statistics)
              p.add <<- as.array(proceeeded$p.post)
              # g.results[4, 2] <<-
            }
          }
        }
        if (model$statid != -1) {
          statistics[id, model$statid + 1] <- statistics[id, model$statid + 1] + 1
        }
        g.results[4, 1] <<- g.results[4, 1] + 1
        return(list(mlik = statistics[id, 1], waic = statistics[id, 2], dic = statistics[id, 3]))
      } else if (exists("hashStat")) {
        if (Nvars.max > Nvars) {
          idd <- as.character(paste(c(model$varcur, array(0, Nvars.max - Nvars)), collapse = ""))
        } else {
          idd <- as.character(paste(c(model$varcur), collapse = ""))
        }

        if (!hash::has.key(key = idd, hash = hashStat)) {
          onerr <- FALSE
          utils::capture.output({
            withRestarts(tryCatch(utils::capture.output({
              fm <- do.call(estimator, c(estimator.args, model$formula))
            })), abort = function() {
              onerr <- TRUE
              fm <- NULL
            })
          }) # fit the modal, get local improvements

          if (!save.beta) {
            hashBuf <- array(data = NA, dim = 3)
          } else {
            if (allow_offsprings == 0 || Nvars.max < Nvars) {
              if (fparam[1] == "Const") {
                linx <- Nvars
                inxx <- which(model$varcur == 1)
              } else {
                linx <- Nvars + 1
                inxx <- c(0, which(model$varcur == 1))
              }
            } else {
              if (fparam[1] == "Const") {
                linx <- Nvars.max
                inxx <- which(model$varcur == 1)
              } else {
                linx <- Nvars.max + 1
                inxx <- c(0, which(model$varcur == 1))
              }
            }

            hashBuf <- array(data = NA, dim = 3 + linx)
          }
          hashBuf[c(2, 3)] <- 100000
          hashBuf[1] <- -100000


          if (!is.null(fm)) {
            hashBuf[2] <- fm$waic[[1]]
            hashBuf[1] <- fm$mlik[[1]]
            hashBuf[3] <- fm$dic[[1]]

            if (save.beta) {
              if (fparam[1] == "Const") {
                if (length(inxx) == length(fm$summary.fixed$mean)) {
                  hashBuf[3 + inxx] <- fm$summary.fixed$mean
                }
              } else {
                if (length(inxx) == length(fm$summary.fixed$mean)) {
                  hashBuf[4 + inxx] <- fm$summary.fixed$mean
                }
              }
            }

            hashStat[idd] <- hashBuf
            #                                                          if(id>1)
            #                                                          {
            #                                                            inxx<-which(model$varcur==1)
            #                                                            if(length(inxx)==length(fm$summary.fixed$mean))
            #                                                              statistics[id,14+inxx]<-fm$summary.fixed$mean
            #                                                          }
            if (fm$waic[[1]] < g.results[2, 1] && !is.na(fm$waic[[1]])) {
              g.results[2, 1] <<- fm$waic[[1]]
              g.results[2, 2] <<- (id)
            }
            if (fm$mlik[[1]] > g.results[1, 1] && !is.na(fm$mlik[[1]])) {
              g.results[1, 1] <<- fm$mlik[[1]]
              g.results[1, 2] <<- (id)
            }

            if (fm$dic[[1]] < g.results[3, 1] && !is.na(fm$dic[[1]])) {
              g.results[3, 1] <<- fm$dic[[1]]
              g.results[3, 2] <<- (id)
            }

            g.results[4, 2] <<- g.results[4, 2] + 1
          }
        }
        g.results[4, 1] <<- g.results[4, 1] + 1
        if (hash::has.key(hash = hashStat, key = idd)) {
          hasRes <- hash::values(hashStat[idd])
        } else {
          hasRes <- c(-10000, 10000, 10000)
        }
        if (g.results[4, 2] %% recalc.margin == 0) {
          proceeeded <- post_proceed_results_hash(hashStat)
          p.add <<- as.array(proceeeded$p.post)
          # g.results[4, 2] <<-
        }
        return(list(mlik = hasRes[1], waic = hasRes[2], dic = hasRes[3]))
      } else {
        utils::capture.output({
          withRestarts(tryCatch(utils::capture.output({
            fm <- do.call(estimator, c(estimator.args, model$formula))
          })), abort = function() {
            onerr <- TRUE
            fm <- NULL
          })
        }) # fit the modal, get local improvements

        if (!is.null(fm)) {
          if (fm$waic[[1]] < g.results[2, 1] && !is.na(fm$waic[[1]])) {
            g.results[2, 1] <<- fm$waic[[1]]
            g.results[2, 2] <<- (id)
          }
          if (fm$mlik[[1]] > g.results[1, 1] && !is.na(fm$mlik[[1]])) {
            g.results[1, 1] <<- fm$mlik[[1]]
            g.results[1, 2] <<- (id)
          }

          if (fm$dic[[1]] < g.results[3, 1] && !is.na(fm$dic[[1]])) {
            g.results[3, 1] <<- fm$dic[[1]]
            g.results[3, 2] <<- (id)
          }
          g.results[4, 1] <<- g.results[4, 1] + 1
          g.results[4, 2] <<- g.results[4, 2] + 1
          return(list(mlik = fm$mlik[[1]], waic = fm$waic[[1]], dic = fm$dic[[1]]))
        } else {
          g.results[4, 1] <<- g.results[4, 1] + 1
          return(list(mlik = -Inf, waic = Inf, dic = Inf))
        }
      }
    }
    g.results[4, 1] <<- g.results[4, 1] + 1
    g.results[4, 2] <<- g.results[4, 2] + 1
    return(list(mlik = -Inf, waic = Inf, dic = Inf))
  }
)
