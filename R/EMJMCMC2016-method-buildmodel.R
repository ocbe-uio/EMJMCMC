EMJMCMC2016$methods(
  # build a new model to draw from the old one return selection probabilities
  buildmodel = function(varcur.old, statid, shift.cpu, max.cpu, switch.type, min.N, max.N, changeble.coord = NULL) {
    vect <- vector(length = max.cpu, mode = "list")
    shift <- 0
    for (cpu in 1:(max.cpu))
    {
      set.seed(stats::runif(1, 1, 10000), kind = NULL, normal.kind = NULL)
      if (!is.null(varcur.old) && length(varcur.old == Nvars)) {
        varcur.old[which(is.na(varcur.old))] <- 1
        varcur <- varcur.old
      } else {
        varcur <- stats::rbinom(n = (Nvars), size = 1, 0.5)
        varcur.old <- varcur
      }

      changeble <- FALSE
      if (!is.null(changeble.coord)) {
        changeble <- TRUE
      }
      if (switch.type == 1) # random size random N(x) # try avoiding when doing mcmc rather than optimization
        {
          log.mod.switch.prob <- ifelse(max.N < min.N, 0, log(1 / (max.N - min.N + 1)))
          KK <- floor(stats::runif(n = 1, min.N, max.N + 0.999999999))
          log.mod.switch.prob <- log.mod.switch.prob + KK * log(truncfactorial(Nvars - KK + 1) / truncfactorial(Nvars))
          log.mod.switchback.prob <- log.mod.switch.prob
          change.buf <- array(data = 0, dim = Nvars)
          if (changeble) {
            for (ttt in 1:KK)
            {
              iid <- floor(stats::runif(n = 1, min = 1, max = Nvars + 0.999999999))
              if (changeble.coord[iid] == 1) {
                next
              }
              change.buf[iid] <- 1
              varcur[iid] <- stats::rbinom(n = 1, size = 1, prob = round(p.add[iid], digits = 8))
              log.mod.switch.prob <- 0 # log.mod.switch.prob + log(dbinom(x = varcur[iid],size = 1,prob = p.add[iid])) # this is one of the pathes to get there in general
              log.mod.switchback.prob <- 0 # log.mod.switchback.prob + log(dbinom(x = 1-varcur[iid],size = 1,prob = p.add[iid])) # but this is one of the ways to get back only
            }
          } else {
            iid <- floor(stats::runif(n = KK, min = 1, max = Nvars + 0.999999999))
            change.buf[iid] <- 1
            varcur[iid] <- stats::rbinom(n = KK, size = 1, prob = round(p.add, digits = 8))
          }
        } else if (switch.type == 2) # fixed sized inverse N(x)
        {
          if (min.N != max.N) {
            min.N <<- max.N
          }
          log.mod.switch.prob <- ifelse(max.N < min.N, 0, log(1 / (max.N - min.N + 1)))
          KK <- max.N
          log.mod.switch.prob <- log.mod.switch.prob + KK * log(truncfactorial(Nvars - KK + 1) / truncfactorial(Nvars))
          log.mod.switchback.prob <- log.mod.switch.prob

          change.buf <- array(data = 0, dim = Nvars)
          if (changeble) {
            for (ttt in 1:KK)
            {
              iid <- floor(stats::runif(n = 1, min = 1, max = Nvars + 0.999999999))
              if (change.buf[iid] == 1 || changeble.coord[iid] == 1) {
                KK <- KK + 1
                next
              }
              change.buf[iid] <- 1
              varcur[iid] <- 1 - varcur[iid]
            }
          } else {
            iid <- floor(stats::runif(n = KK, min = 1, max = Nvars + 0.999999999))
            change.buf[iid] <- 1
            varcur[iid] <- 1 - varcur[iid]
          }
        } else if (switch.type == 3) # random sized inverse N(x)
        {
          log.mod.switch.prob <- ifelse(max.N < min.N, 0, log(1 / (max.N - min.N + 1)))
          KK <- floor(stats::runif(n = 1, min.N, max.N + 0.999999999))
          log.mod.switch.prob <- log.mod.switch.prob + KK * log(truncfactorial(Nvars - KK + 1) / truncfactorial(Nvars))
          log.mod.switchback.prob <- log.mod.switch.prob
          change.buf <- array(data = 0, dim = Nvars)
          if (changeble) {
            for (ttt in 1:KK)
            {
              iid <- floor(stats::runif(n = 1, min = 1, max = Nvars + 0.999999999))
              if (change.buf[iid] == 1 || changeble.coord[iid] == 1) {
                KK <- KK + 1
                next
              }
              change.buf[iid] <- 1
              varcur[iid] <- 1 - varcur[iid]
            }
          } else {
            iid <- floor(stats::runif(n = KK, min = 1, max = Nvars + 0.999999999))
            change.buf[iid] <- 1
            varcur[iid] <- 1 - varcur[iid]
          }
        } else if (switch.type == 4) # fixed N(x) for reverse from type 2 swaps
        {
          if (min.N != max.N) {
            min.N <<- max.N
          }
          log.mod.switch.prob <- ifelse(max.N < min.N, 0, log(1 / (max.N - min.N + 1)))
          KK <- max.N
          log.mod.switch.prob <- log.mod.switch.prob + KK * log(truncfactorial(Nvars - KK + 1) / truncfactorial(Nvars))
          log.mod.switchback.prob <- log.mod.switch.prob
          ids <- which(changeble.coord == 0)
          change.buf <- array(data = 0, dim = Nvars)
          if (changeble) {
            for (ttt in KK)
            {
              iid <- floor(stats::runif(n = 1, min = 1, max = Nvars + 0.999999999))
              if (change.buf[iid] == 1 || changeble.coord[iid] == 1) {
                KK <- KK + 1
                next
              }
              change.buf[iid] <- 1
              varcur[iid] <- 1 - varcur[iid]
            }
          } else {
            iid <- floor(stats::runif(n = KK, min = 1, max = Nvars + 0.999999999))
            change.buf[iid] <- 1
            varcur[iid] <- 1 - varcur[iid]
          }
        } else if (switch.type == 5) {
        change.buf <- array(data = 0, dim = (Nvars))
        log.mod.switch.prob <- 0
        log.mod.switchback.prob <- 0
        add <- 0
        if (cpu + shift > Nvars) {
          shift <- 1 - cpu
          varcur.old <- rep(0, times = Nvars)
        }
        if (varcur.old[cpu + shift] != 1) {
          varcur <- varcur.old
          varcur[cpu + shift] <- 1
        } else {
          if (cpu + shift < Nvars) {
            shift <- shift + 1
          }
          while (varcur.old[cpu + shift] == 1) {
            shift <- shift + 1
            if (cpu + shift >= Nvars) {
              shift <- (Nvars - cpu)
              break
            }
          }
          varcur <- varcur.old
          varcur[cpu + shift] <- 1
        }
      } else if (switch.type == 6) {
        change.buf <- array(data = 0, dim = (Nvars))
        log.mod.switch.prob <- 0
        log.mod.switchback.prob <- 0
        add <- 0
        if (cpu + shift > Nvars) {
          shift <- 1 - cpu
          varcur.old <- rep(1, times = Nvars)
        }

        if (varcur.old[cpu + shift] != 0) {
          varcur <- varcur.old
          varcur[cpu + shift] <- 0
        } else {
          if (cpu + shift < Nvars) {
            shift <- shift + 1
          }
          while (varcur.old[cpu + shift] == 0) {
            shift <- shift + 1
            if (cpu + shift >= Nvars) {
              shift <- (Nvars - cpu)
              break
            }
          }
          varcur <- varcur.old
          varcur[cpu + shift] <- 0
        }
      } else if (switch.type == 7) {
        change.buf <- array(data = 0, dim = (Nvars))
        log.mod.switch.prob <- 0
        log.mod.switchback.prob <- 0
        vec <- dectobit(cpu + shift.cpu)
        varcur <- c(array(0, dim = (Nvars - length(vec))), vec) # issues here
      } else if (switch.type == 8) {
        log.mod.switch.prob <- 0
        log.mod.switchback.prob <- 0
        varcur <- varcur.old
        change.buf <- array(data = 0, dim = Nvars)
        # if(printable.opt)print("type 8 invoked")
        # if(printable.opt)print(varcur)
      } else if (switch.type == 9) {
        log.mod.switch.prob <- 0
        log.mod.switchback.prob <- 0
        change.buf <- array(data = 1, dim = Nvars)
        changevar <- stats::rbinom(n = Nvars, size = 1, prob = prand)
        varcur <- varcur.old
        varcur[which(changevar == 1)] <- (1 - varcur[which(changevar == 1)])
        log.mod.switch.prob <- prand^length(which(changevar == 1))
      } else {
        log.mod.switch.prob <- 0
        log.mod.switchback.prob <- 0
        change.buf <- array(data = 1, dim = Nvars)
        varcur <- stats::rbinom(n = Nvars, size = 1, prob = round(p.add, digits = 8))
        # if(printable.opt)print("type 8 invoked")
        # if(printable.opt)print(varcur)
      }

      #     for(g in 1:max(isobsbinary))
      #     {
      #       if(length(varcur[which(isobsbinary == g && varcur %in% c(1,3))])==length(which(isobsbinary == g)))
      #         varcur[which(isobsbinary == g && varcur %in% c(1,3))[1]]=0
      #       if(length(varcur[which(isobsbinary == g && varcur %in% c(2,3))])==length(which(isobsbinary == g)))
      #         varcur[which(isobsbinary == g && varcur %in% c(2,3))[1]]=0
      #     }



      covobs <- if (fparam[1] == "Const") fparam[which(varcur[-1] == 1) + 1] else fparam[which(varcur == 1)]


      # obsconst<-2*as.integer((varcur[1] %in% c(1,3)) || length(covobs) == 0 || (length(covobs)==length(which(isobsbinary != 0))) ) -1

      obsconst <- ifelse(fparam[1] == "Const", 2 * as.integer((varcur[1])) - 1, 1)


      id <- bittodec(varcur)
      id <- id + 1



      if (ifelse(exists("statistics1"), is.na(statistics1[id, 1]), ifelse(exists("statistics"), is.na(statistics[id, 1, ]), ifelse(exists("hashStat"), !hash::has.key(hash = hashStat, key = paste(varcur, collapse = "")), TRUE)))) #||TRUE)
        {
          # formula <- NULL
          # formula <- stats::as.formula(ifelse(length(covobs)>0,(stri_join( stringi::stri_flatten(fobserved[1]), " ~ ",obsconst,"+",  stringi::stri_flatten(covobs, collapse=" + "), latent.formula)),(stri_join( stringi::stri_flatten(fobserved[1]), " ~ ",obsconst, latent.formula))))
          #
          # if(is.null(formula)){
          #   formula <- stats::as.formula(ifelse(length(covobs)>0,(stri_join( stringi::stri_flatten(fobserved[1]), " ~ ",obsconst,"+",  stringi::stri_flatten(covobs, collapse=" + "))),(stri_join( stringi::stri_flatten(fobserved[1]), " ~ ",obsconst))))
          #
          # }
          formula <- NULL
          utils::capture.output({
            withRestarts(tryCatch(utils::capture.output({
              formula <- stats::as.formula(paste(paste(fobserved[1]), " ~ ", obsconst, ifelse(length(covobs) > 0, " + ", ""), paste(covobs, collapse = " + "), latent.formula))
            })), abort = function() {
              onerr <- TRUE
              fm <- NULL
            })
          })
          if (is.null(formula)) {
            formula <- stats::as.formula(paste(paste(fobserved[1]), " ~ ", obsconst, ifelse(length(covobs) > 0, " + ", ""), paste(covobs, collapse = " + ")))
          }
        } else {
        formula <- NULL
      }

      vect[[cpu]] <- list(formula = formula, varcur = varcur, statid = statid, changed = change.buf, log.mod.switch.prob = log.mod.switch.prob, log.mod.switchback.prob = log.mod.switchback.prob)
    }

    return(vect)
  }
)
