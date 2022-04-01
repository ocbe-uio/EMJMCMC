EMJMCMC2016$methods(
  # forward backward random dance
  forw_backw_walk = function(model) {
    varcur <- stats::rbinom(n = Nvars, size = 1, prob = stats::runif(n = 1, min = 0, max = model$p1))
    mlikcur <- -Inf
    waiccur <- Inf
    for (i in 1:model$steps)
    {
      fff <- forward_selection(list(varcur = stats::rbinom(n = Nvars, size = 1, prob = stats::runif(n = 1, min = 0, max = model$p1)), mlikcur = -Inf, waiccur = Inf, locstop = FALSE, statid = -1))
      if (objective == 1) {
        if (fff$mlikglob > mlikcur) {
          mlikcur <- fff$mlikglob
          varcur <- fff$varglob
          waiccur <- fff$waicglob
        }
      } else {
        if (fff$waicglob < waiccur) {
          mlikcur <- fff$mlikglob
          varcur <- fff$varglob
          waiccur <- fff$waicglob
        }
      }
      set.seed(i * model$steps)
      bbb <- backward_selection(list(varcur = stats::rbinom(n = Nvars, size = 1, prob = stats::runif(n = 1, min = model$p1, max = 1)), mlikcur = -Inf, waiccur = Inf, locstop = FALSE, statid = -1))
      if (objective == 1) {
        if (bbb$mlikglob > mlikcur) {
          mlikcur <- bbb$mlikglob
          varcur <- bbb$varglob
          waiccur <- bbb$waicglob
        }
      } else {
        if (bbb$waicglob < waiccur) {
          mlikcur <- bbb$mlikglob
          varcur <- bbb$varglob
          waiccur <- bbb$waicglob
        }
      }
    }
    if (model$reverse == FALSE) {
      vect <- buildmodel(max.cpu = 1, varcur.old = varcur, statid = -1, switch.type = type.randomize, min.N = min.N.randomize, max.N = max.N.randomize)

      varcur <- vect[[1]]$varcur
      # if(printable.opt)print(varcur)

      cluster <- TRUE



      if (is.null(vect[[1]]$formula)) {
        cluster <- FALSE
        if (printable.opt) print("!!!!Back Forw reverse model already estimated!!!!")
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
      model.probs <- calculate.move.logprobabilities(switch.type = type.randomize, varold = varcur, varnew = model$varold, min.N = min.N.randomize, max.N = max.N.randomize)
      model.prob <- model.probs$log.switch.forw.prob
      model.prob.fix <- model.probs$log.switch.back.prob
    }

    return(list(varcur = varcur, waiccur = waiccur, mlikcur = mlikcur, log.prob.cur = model.prob, log.prob.fix = model.prob.fix, varglob = varcur, waicglob = waiccur, mlikglob = mlikcur))
  }
)
