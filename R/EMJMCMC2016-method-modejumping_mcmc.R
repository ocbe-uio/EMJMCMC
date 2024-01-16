EMJMCMC2016$methods(
  # global emjmcmc procedure for model selection (hyper heuristic logic in terms of COP)
  modejumping_mcmc = function(glob.model) {
    stm <- proc.time()
    if (printable.opt) print("Begin model selection EMJMCMC2016 procedure")
    set.seed(stats::runif(n = 1, min = 1, max = seed), kind = NULL, normal.kind = NULL)
    acc_moves <- 1
    accept_old <- 1
    distrib_of_proposals <- glob.model$distrib_of_proposals
    distrib_of_neighbourhoods <- glob.model$distrib_of_neighbourhoods
    # do the search and simulations accross the modes
    g.results[4, 1] <<- 0
    g.results[4, 2] <<- 0

    if (glob.model$presearch) {
      forward_selection(list(varcur = rep(0, length(fparam.example)), mlikcur = -Inf, waiccur = Inf, locstop = glob.model$locstop, statid = -1))
      backward_selection(list(varcur = rep(1, length(fparam.example)), mlikcur = -Inf, waiccur = Inf, locstop = glob.model$locstop, statid = -1))

      if (exists("statistics1") && recalc.margin < 2^Nvars) {
        p.add <<- as.array(post_proceed_results(statistics1)$p.post)
      } else if (exists("hashStat") && recalc.margin < 2^Nvars) {
        p.add <<- as.array(post_proceed_results_hash(hashStat)$p.post)
      }
    } else {
      p.add <<- array(data = 0.1, Nvars)
      p.post <- array(data = 0.1, Nvars)
      vec <- stats::rbinom(n = Nvars, size = 1, prob = 0.0000001) # generate an initial solution
      varcur <- c(array(0, dim = (Nvars - length(vec))), vec)
    }
    waiccur <- Inf
    waicglob <- Inf
    mlikcur <- -Inf
    mlikglob <- -Inf
    ratcur <- -Inf

    # set up initial parameters
    if (is.null(glob.model$varcur)) {
      if ((!is.na(g.results[1, 2])) && g.results[1, 2] > 0) {
        vec <- dectobit(g.results[1, 2] - 1)
        varcur <- c(array(0, dim = (Nvars - length(vec))), vec)
        waiccur <- g.results[2, 1]
        waicglob <- g.results[2, 1]
        mlikcur <- g.results[1, 1]
        mlikglob <- g.results[1, 1]
        ratcur <- g.results[1, 1]

        print(paste("initial solution is set with mlik of ", mlikcur))
      } else {
        vec <- stats::rbinom(n = Nvars, size = 1, prob = 0.5) # generate an initial solution
        varcur <- c(array(0, dim = (Nvars - length(vec))), vec)
      }
    } else if (length(glob.model$varcur[which(glob.model$varcur %in% c(0, 1))]) == Nvars) {
      varcur <- glob.model$varcur
    } else {
      if (printable.opt) print("Incorrect initial solution set be the user, a random one is generated")
      vec <- stats::rbinom(n = Nvars, size = 1, prob = 0.5) # generate an initial solution
      varcur <- c(array(0, dim = (Nvars - length(vec))), vec)
    }

    varcurb <- varcur
    varglob <- varcur
    modglob <- NULL

    p1 <- array(data = 0.0001, dim = Nvars)
    p2 <- array(data = 0.9999, dim = Nvars)
    j <- 0
    j.a <- 0
    p.post <- array(data = 1, dim = Nvars)
    waiccur <- Inf
    waicglob <- Inf
    mlikcur <- -Inf
    mlikglob <- -Inf
    ratcur <- -Inf
    fm <- NULL
    eps.emp <- normprob(p1, p2)
    max.cpu.buf <- max.cpu
    delta.time <- 0
    LocImprove <<- as.array(0)
    LocNeighbor <- 0
    max.cpu.buf <- max.cpu.glob

    while ((eps.emp >= glob.model$eps || j <= glob.model$maxit || j <= glob.model$burnin) && delta.time < glob.model$max.time && g.results[4, 1] <= glob.model$trit && g.results[4, 2] <= glob.model$trest) {
      p1 <- p.post / acc_moves
      set.seed(stats::runif(n = 1, min = 1, max = seed * 100), kind = NULL, normal.kind = NULL)
      LocImprove <<- as.array(sample(x = 5, size = 1, prob = distrib_of_proposals) - 1)
      LocNeighbor <- (sample(x = 7, size = 1, prob = distrib_of_neighbourhoods[LocImprove + 1, ]))
      switch.type.glob.buf <- LocNeighbor
      switch.type.buf <- LocNeighbor
      if (LocNeighbor == 7) {
        switch.type.glob.buf <- 9
        switch.type.buf <- 9
      }

      # if(printable.opt)print(LocImprove)
      j <- j + 1
      j.a <- j.a + 1
      if (glob.model$print.freq > 0 && j %% glob.model$print.freq == 0) {
        cat(
          formatC(j, width = 4L, format = "d"),
          "iterations completed up to now after",
          formatC(delta.time, digits = 6L, flag = "-", format = "f"),
          "cpu minutes",
          "best MLIK found",
          formatC(g.results[1, 1], digits = 3L, flag = "-", format = "f"),
          "current mlik found",
          formatC(mlikcur, digits = 3L, flag = "-", format = "f"),
          "current acceptance ratio",
          formatC(acc_moves / j.a, digits = 6L, flag = "-", format = "f"), "\n"
        )
      }
      if (j %% 100 == 0) {
        seed <<- stats::runif(n = 1, min = 0, max = 100000)
      }
      # the small part of the code to be upgraded at least slightly
      if (allow_offsprings %in% c(1, 2) && j %% mutation_rate == 0 && (j <= last.mutation || Nvars != Nvars.max)) {
        if (Nvars > Nvars.max || j == mutation_rate) {
          # do the stuff here
          if (j == mutation_rate) {
            fparam.pool <<- c(fparam.pool, filtered)
          }
          to.del <- which(p.add < p.allow.tree)
          if (length(to.del) == Nvars) {

            to.del <- to.del[-sample(x = Nvars, size = sample(x = Nvars - 1, size = 1), prob = p.add + p.epsilon)]
          }
          if (length(to.del) < Nvars - Nvars.max) {
            tdl.id <- order(p.add, decreasing = T)
            to.del <- to.del[-tdl.id[1:Nvars.max]]
          }
          if (glob.model$print.freq > 0L) {
            message("Data filtered! Insignificant variables deleted!")
          }
          if (length(to.del) > 0) {
            hash::clear(hashStat)
            hashStat <- hash::hash()
            fparam <<- fparam[-to.del]
            Nvars <<- length(fparam)
            Nvars.init <<- Nvars
            p.add <<- p.add[-to.del]
            p.post <- array(data = 1, dim = Nvars)
            varcurb <- varcurb[1:Nvars]
            varcand <- varcurb[1:Nvars]
            varglob <- varcurb[1:Nvars]
            p1 <- array(0, dim = (Nvars))
            p2 <- array(1, dim = (Nvars))
            acc_moves <- 1
            j.a <- 1
          }
        } else {
          if (Nvars >= Nvars.max) {
            idmut <- (which(p.add[(Nvars.init + 1):Nvars] <= p.allow.replace) + Nvars.init)
            lidmut <- length(idmut) # maximal number of covariates that can die out
            if (lidmut > 0) {
              p.del <- (lidmut - sum(p.add[idmut])) / lidmut
              lidmut <- stats::rbinom(n = 1, size = lidmut, prob = p.del)
            }
          } else {
            idmut <- (Nvars + 1):Nvars.max
            lidmut <- Nvars.max - Nvars
          }

          for (idel in 1:lidmut) {
            p.del <- 1 - (sum(p.add)) / Nvars
            mother <- ifelse(stats::runif(n = 1, min = 0, max = 1) <= p.del, fparam[which(rmultinom(n = 1, size = 1, prob = p.add / 2) == 1)], fparam[stats::runif(n = 1, min = 1, max = Nvars.init)])
            ltreem <- stringi::stri_length(mother)
            mother <- stringi::stri_sub(mother, from = 2, to = ltreem)

            if (allow_offsprings == 1) {
              sjm <- sum(stringi::stri_count_fixed(str = mother, pattern = c("&", "|")))
            } else {
              sjm <- sum(stringi::stri_count_fixed(str = mother, pattern = c("+", "*")))
            }

            if (sjm <= max.tree.size) {

              # p.del<-1-(sum(p.add))/Nvars
              father <- ifelse(stats::runif(n = 1, min = 0, max = 1) <= p.del, fparam.pool[stats::runif(n = 1, min = 1, max = length(fparam.pool))], fparam[which(rmultinom(n = 1, size = 1, prob = p.add / 2) == 1)])
              ltreef <- stringi::stri_length(father)
              father <- stringi::stri_sub(father, from = 2, to = ltreef)

              if (allow_offsprings == 1) {
                sjf <- sum(stringi::stri_count_fixed(str = father, pattern = c("&", "|")))
              } else {
                sjf <- sum(stringi::stri_count_fixed(str = father, pattern = c("+", "*")))
              }

              if (sjm + sjf + 1 <= max.tree.size) {
                if (allow_offsprings == 1) {
                  if (!grepl(father, mother, fixed = T) && !grepl(mother, father, fixed = T)) {
                    proposal <- stri_paste(paste(ifelse(stats::runif(n = 1, min = 0, max = 1) < p.nor, "I((1-", "I(("), mother, sep = ""), paste(ifelse(stats::runif(n = 1, min = 0, max = 1) < p.nor, "(1-", "("), father, "))", sep = ""), sep = ifelse(stats::runif(n = 1, min = 0, max = 1) < p.and, ")&", ")|"))
                  } else {
                    if (max(sjm, sjf) > 1) {
                      t.d <- sample(size = 1, x = (max(sjm, sjf) + 1))
                      if (sjm >= sjf) {
                        loc <- c(1, stri_locate_all(str = mother, regex = "\\&|\\||\\*|\\+")[[1]][, 1], stringi::stri_length(mother))
                        proposal <- stri_paste(stringi::stri_sub(mother, from = 1, to = loc[t.d] - 1), stringi::stri_sub(mother, from = (loc[t.d + 1] + (t.d == 1)), to = stringi::stri_length(mother)))
                      } else {
                        loc <- c(1, stri_locate_all(str = father, regex = "\\&|\\||\\*|\\+")[[1]][, 1], stringi::stri_length(father))
                        proposal <- stri_paste(stringi::stri_sub(father, from = 1, to = loc[t.d] - 1), stringi::stri_sub(father, from = (loc[t.d + 1] + (t.d == 1)), to = stringi::stri_length(father)))
                      }

                      diffs <- (stringi::stri_count_fixed(str = proposal, pattern = "(") - stringi::stri_count_fixed(str = proposal, pattern = ")"))
                      if (diffs > 0) {
                        proposal <- stri_paste(proposal, stri_paste(rep(")", diffs), collapse = ""), collapse = "")
                      }
                      if (diffs < 0) {
                        proposal <- stri_paste(stri_paste(rep("(", -diffs), collapse = ""), proposal, collapse = "")
                      }
                      proposal <- stri_paste("I", proposal)
                      # print(paste("&&&&&&",mother,"ssss",father))
                    } else {
                      proposal <- stri_paste("I", mother)
                    }
                  }
                  # proposal<-stri_paste("I",mother)
                } else {
                  proposal <- stri_paste(paste(ifelse(stats::runif(n = 1, min = 0, max = 1) < p.nor, "I(", "I(-"), mother, sep = ""), paste("(", father, "))", sep = ""), sep = ifelse(stats::runif(n = 1, min = 0, max = 1) < p.and, "*", "+"))
                  proposal <- stri_paste("I(", sigmas[sample(x = length(sigmas), size = 1, prob = sigmas.prob)], "(", proposal, "))", sep = "")
                  while ((proposal %in% fparam)) {
                    proposal <- stri_paste("I(", sigmas[sample(x = length(sigmas), size = 1, prob = sigmas.prob)], "(", proposal, "))", sep = "")
                  }
                }
              } else {
                t.d <- sample(size = 1, x = (max(sjm, sjf) + 1))
                if (sjm >= sjf) {
                  loc <- c(1, stri_locate_all(str = mother, regex = "\\&|\\||\\*|\\+")[[1]][, 1], stringi::stri_length(mother))
                  proposal <- stri_paste(stringi::stri_sub(mother, from = 1, to = loc[t.d] - 1), stringi::stri_sub(mother, from = (loc[t.d + 1] + (t.d == 1)), to = stringi::stri_length(mother)))
                } else {
                  loc <- c(1, stri_locate_all(str = father, regex = "\\&|\\||\\*|\\+")[[1]][, 1], stringi::stri_length(father))
                  proposal <- stri_paste(stringi::stri_sub(father, from = 1, to = loc[t.d] - 1), stringi::stri_sub(father, from = (loc[t.d + 1] + (t.d == 1)), to = stringi::stri_length(father)))
                }

                diffs <- (stringi::stri_count_fixed(str = proposal, pattern = "(") - stringi::stri_count_fixed(str = proposal, pattern = ")"))
                if (diffs > 0) {
                  proposal <- stri_paste(proposal, stri_paste(rep(")", diffs), collapse = ""), collapse = "")
                }
                if (diffs < 0) {
                  proposal <- stri_paste(stri_paste(rep("(", -diffs), collapse = ""), proposal, collapse = "")
                }
                proposal <- stri_paste("I", proposal)
                # print(paste("!!!!!",mother,"ssss",father))
              }
              # maybe check correlations here
              if ((!(proposal %in% fparam)) && Nvars < Nvars.max) {
                # if(cor())
                fparam <<- c(fparam, proposal)
                Nvars <<- as.integer(Nvars + 1)
                p.add <<- as.array(c(p.add, p.allow.replace))
                p.post <- as.array(c(p.post, 1))
                if (printable.opt) {
                  print(paste("mutation happended ", proposal, " tree  added"))
                }
              } else if (!(proposal %in% fparam)) {
                to.del <- (which(p.add[(Nvars.init + 1):Nvars] < p.allow.replace) + Nvars.init)
                lto.del <- length(x = to.del)
                if (lto.del > 0) {
                  id.replace <- to.del[round(stats::runif(n = 1, min = 1, max = lto.del))]
                  if (printable.opt) {
                    print(paste("mutation happended ", proposal, " tree  replaced ", fparam[id.replace]))
                  }
                  fparam[id.replace] <<- proposal
                  keysarr <- as.array(hash::keys(hashStat))
                  p.add[id.replace] <<- p.allow.replace
                  for (jjj in 1:length(keysarr))
                  {
                    if (length(keysarr > 0)) {
                      if (stringi::stri_sub(keysarr[jjj], from = id.replace, to = id.replace) == "1") {
                        del(x = keysarr[jjj], hash = hashStat)
                      }
                    }
                  }
                }
              }
            }
          }

          varcurb <- c(varcurb, array(1, dim = (Nvars - length(varcurb))))
          varcand <- c(varcand, array(1, dim = (Nvars - length(varcand))))
          varglob <- c(varglob, array(1, dim = (Nvars - length(varglob))))
          p.post <- array(1, dim = (Nvars))
          p1 <- c(p1, array(0, dim = (Nvars - length(p1))))
          p2 <- c(p1, array(1, dim = (Nvars - length(p1))))
          acc_moves <- 1
          j.a <- 1
        }
      } else if (allow_offsprings == 3 && j %% mutation_rate == 0 && (j <= last.mutation || Nvars != Nvars.max)) {
        # if(latnames[1]!="")

        # perform preliminary filtration here
        if (Nvars > Nvars.max || j == mutation_rate) {
          # pool.cor.prob = T
          # do the stuff here
          if (j == mutation_rate) {
            fparam.pool <<- unique(c(fparam.pool, filtered))
            if (!pool.cor.prob) {
              pool.probs <- array(data = 1 / length(fparam.pool), dim = length(fparam.pool))
            } else {
              fobserved.cleaned <- fobserved
              fobserved.cleaned <- stri_replace(str = fobserved.cleaned, fixed = "I(", replacement = "")
              fobserved.cleaned <- stri_replace(str = fobserved.cleaned, fixed = ")", replacement = "")
              fparam.pool.cleaned <- fparam.pool
              fparam.pool.cleaned <- stri_replace(str = fparam.pool.cleaned, fixed = "I(", replacement = "")
              fparam.pool.cleaned <- stri_replace(str = fparam.pool.cleaned, fixed = ")", replacement = "")
              pool.probs <- abs(cor(estimator.args$data[[fobserved.cleaned]], estimator.args$data[, which(fparam.pool.cleaned %in% names(estimator.args$data))])) + p.epsilon
              rm(fobserved.cleaned)
              rm(fparam.pool.cleaned)
            }
          }
          to.del <- which(p.add < p.allow.tree)
          if (length(to.del) == Nvars) {
            to.del <- to.del[-sample(x = Nvars, size = sample(x = Nvars - 1, size = 1), prob = p.add + p.epsilon)]
          }
          if (length(to.del) < Nvars - Nvars.max) {
            tdl.id <- order(p.add, decreasing = T)
            to.del <- to.del[-tdl.id[1:Nvars.max]]
          }
          if (glob.model$print.freq > 0L) {
            message("Data filtered! Insignificant variables deleted!")
          }
          if (length(to.del) > 0) {
            hash::clear(hashStat)
            hashStat <- hash::hash()
            fparam <<- fparam[-to.del]
            Nvars <<- length(fparam)
            Nvars.init <<- Nvars
            p.add <<- p.add[-to.del]
            p.post <- array(data = 1, dim = Nvars)
            varcurb <- varcurb[-to.del]
            varcand <- varcurb[-to.del]
            varglob <- varcurb[-to.del]
            p1 <- array(0, dim = (Nvars))
            p2 <- array(1, dim = (Nvars))
            acc_moves <- 1
            j.a <- 1
          }
        } else {
          if (Nvars >= Nvars.max) {
            if (keep.origin) {
              idmut <- (which(p.add[(Nvars.init + 1):Nvars] <= p.allow.replace) + Nvars.init)
              lidmut <- length(idmut) # maximal number of covariates that can die out
              if (lidmut > 0) {
                p.del <- (lidmut - sum(p.add[idmut])) / lidmut
                lidmut <- stats::rbinom(n = 1, size = lidmut, prob = p.del)
              }
            } else {
              idmut <- which(p.add <= p.allow.replace)
              lidmut <- length(idmut) # maximal number of covariates that can die out
              if (lidmut > 0) {
                p.del <- (lidmut - sum(p.add[idmut])) / lidmut
                lidmut <- stats::rbinom(n = 1, size = lidmut, prob = p.del)
              }
            }
            # if(lidmut>0)
            # {
            # p.del<-(lidmut - sum(p.add[idmut]))/lidmut
            # lidmut<-stats::rbinom(n = 1,size = lidmut,prob = p.del)
            # }
          } else {
            idmut <- (Nvars + 1):Nvars.max
            lidmut <- Nvars.max - Nvars
          }
          # now having chosen the candidates to be deleted we can propose new variables
          idel <- 1
          repsd <- 1
          while (idel <= lidmut && repsd <= lidmut * 100) {
            repsd <- repsd + 1
            # gen.prob<-c(1,1,1,1,1)#just uniform for now
            action.type <- sample(x = 5, size = 1, prob = gen.prob)


            if (action.type == 1) {
              # mutation (add a leave not in the search space)
              proposal <- fparam.pool[sample(x = length(fparam.pool), size = 1, prob = pool.probs)]
            } else if (action.type == 2) {
              # crossover type of a proposal
              # generate a mother
              # actvars<-which(varcurb==1)
              mother <- ifelse(stats::runif(n = 1, min = 0, max = 1) <= pool.cross, fparam[sample(x = Nvars, size = 1, prob = p.add + p.epsilon)], fparam.pool[sample(x = length(fparam.pool), size = 1, prob = pool.probs)])
              ltreem <- stringi::stri_length(mother)
              mother <- stringi::stri_sub(mother, from = 1, to = ltreem)
              # sjm<-sum(stringi::stri_count_fixed(str = mother, pattern = c("+","*")))
              # generate a father
              father <- ifelse(stats::runif(n = 1, min = 0, max = 1) <= pool.cross, fparam[sample(x = Nvars, size = 1, prob = p.add + p.epsilon)], fparam.pool[sample(x = length(fparam.pool), size = 1, prob = pool.probs)])
              ltreef <- stringi::stri_length(father)
              father <- stringi::stri_sub(father, from = 1, to = ltreef)
              # sjf<-sum(stringi::stri_count_fixed(str = father, pattern = c("+","*")))

              proposal <- stri_paste("I(", stri_paste(mother, father, sep = "*"), ")", sep = "")
            } else if (action.type == 3) {
              proposal <- ifelse(stats::runif(n = 1, min = 0, max = 1) <= pool.cross, fparam[sample(x = Nvars, size = 1, prob = p.add + p.epsilon)], fparam.pool[sample(x = length(fparam.pool), size = 1, prob = pool.probs)])
              proposal <- stri_paste("I(", sigmas[sample(x = length(sigmas), size = 1, prob = sigmas.prob)], "(", proposal, "))", sep = "")
            } else if (action.type == 4) {


              # select a subset for the projection
              spad <- sum(p.add)
              actvars <- which(stats::rbinom(n = length(fparam), size = 1, prob = p.add / spad + p.epsilon) == 1)


              if (length(actvars) <= 1) {
                proposal <- fparam.pool[sample(x = length(fparam.pool), size = 1)]
              } else {
                utils::capture.output(withRestarts(tryCatch(
                  {
                    # print(stats::as.formula(stri_paste(fobserved,"~ 1 +",paste0(fparam[actvars],collapse = "+"))))
                    # get the projection coefficients as the posterior mode of the fixed effects
                    if (deep.method == 1) {
                      bet.act <- do.call(.self$estimator, c(estimator.args, stats::as.formula(stri_paste(fobserved, "~ 1 +", paste0(fparam[actvars], collapse = "+")))))$summary.fixed$mean
                      nab <- which(is.na(bet.act))
                      if (length(nab) > 0) {
                        bet.act[nab] <- 0
                      }
                      bet.act <- round(bet.act, digits = 8)
                      bet.act <- stri_paste("m(", bet.act, ",", sep = "")
                      # make a projection
                      bet.act <- stri_paste(bet.act, c("1", fparam[actvars]), ")", sep = "")
                      bet.act <- stri_paste("I(", bet.act, ")", sep = "")
                      proposal <- stri_paste("I(", stri_paste(bet.act, collapse = "+"), ")", collapse = "")
                      proposal <- stri_paste("I(", sigmas[sample(x = length(sigmas), size = 1, prob = sigmas.prob)], "(", proposal, "))", sep = "")
                      # print(proposal)
                    } else if (deep.method == 2) {
                      cursigma <- sigmas[sample(x = length(sigmas), size = 1, prob = sigmas.prob)]
                      bet.act <- gnlr(
                        y = data.example[[fobserved]],
                        distribution = estimator.args$distribution,
                        mu = stats::as.formula(stri_paste("~", estimator.args$link, "(", cursigma, "(", "b0 +", paste0("b", 1:length(actvars), "*", fparam[actvars], collapse = "+"), "))")),
                        pmu = rep(0, length(actvars) + 1)
                      )$coefficients


                      nab <- which(is.na(bet.act))
                      if (length(nab) > 0) {
                        bet.act[nab] <- 0
                      }
                      # bet.act<-round(bet.act, digits = 2)
                      bet.act <- round(bet.act, digits = 8)
                      bet.act <- stri_paste("m(", bet.act, ",", sep = "")
                      # make a projection
                      bet.act <- stri_paste(bet.act, c("1", fparam[actvars]), ")", sep = "")
                      bet.act <- stri_paste("I(", bet.act, ")", sep = "")
                      proposal <- stri_paste("I(", stri_paste(bet.act, collapse = "+"), ")", collapse = "")
                      proposal <- stri_paste("I(", cursigma, "(", proposal, "))", sep = "")
                    } else if (deep.method == 3) {
                      cursigma <- sigmas[sample(x = length(sigmas), size = 1, prob = sigmas.prob)]
                      forstr <- stri_paste("~", estimator.args$link, "(", cursigma, "(", "m(0,1) +", paste0("m(", stats::rnorm(length(actvars), 0, 1), ",", fparam[actvars], ")", collapse = "+"), "))")
                      forstr <- stringi::stri_replace_all_fixed(str = forstr, pattern = c("m(-"), replacement = "m(")
                      forstr <- stringi::stri_replace_all_fixed(str = forstr, pattern = c("m("), replacement = "m(b_")
                      nlrr <- gnlr(
                        y = data.example[[fobserved]],
                        distribution = estimator.args$distribution,
                        mu = stats::as.formula(forstr),
                        pmu = stats::rnorm(stringi::stri_count_fixed(str = forstr, pattern = "m("), 0, 0.0001)
                      )
                      beg.rep <- stri_locate_all(str = forstr, fixed = "m(b")[[1]][, 2]
                      end.rep <- stri_locate_all(str = forstr, fixed = ",")[[1]][, 2]

                      bet.act <- nlrr$coefficients

                      nab <- which(is.na(bet.act))
                      if (length(nab) > 0) {
                        bet.act[nab] <- 0
                      }
                      bet.act <- round(bet.act, digits = 8)
                      torepl <- stringi::stri_sub(str = forstr, from = beg.rep, to = end.rep)
                      for (i in 1:length(bet.act)) {
                        forstr <- stringi::stri_replace_all_fixed(str = forstr, pattern = torepl[i], replacement = stri_paste(bet.act[i], ","))
                      }
                      forstr <- stringi::stri_replace_all_fixed(str = forstr, pattern = paste0("~", estimator.args$link), replacement = "I")

                      proposal <- forstr
                    } else {
                      bet.act <- stats::rnorm(n = (length(actvars) + 1), mean = 0, sd = 1)
                      nab <- which(is.na(bet.act))
                      if (length(nab) > 0) {
                        bet.act[nab] <- 0
                      }
                      bet.act <- round(bet.act, digits = 8)
                      bet.act <- stri_paste("m(", bet.act, ",", sep = "")
                      # make a projection
                      bet.act <- stri_paste(bet.act, c("1", fparam[actvars]), ")", sep = "")
                      bet.act <- stri_paste("I(", bet.act, ")", sep = "")
                      proposal <- stri_paste("I(", stri_paste(bet.act, collapse = "+"), ")", collapse = "")
                      proposal <- stri_paste("I(", sigmas[sample(x = length(sigmas), size = 1, prob = sigmas.prob)], "(", proposal, "))", sep = "")
                    }
                  },
                  error = function(err) {
                    # print(err)
                    proposal <- fparam.pool[sample(x = length(fparam.pool), size = 1)]
                  },
                  finally = {
                    proposal <- proposal
                  }
                )))
                # print(proposal)

                if (is.na(proposal)) {
                  print(fparam[actvars])
                  print(actvars)
                  print(bet.act)
                }
              }
            } else if (action.type == 5) {
              # reduce an operator fparam[idel]
              # print("reduction")
              # print(idel)
              if (idel > length(fparam)) {
                proposal <- fparam.pool[sample(x = length(fparam.pool), size = 1)]
              } else {
                cpm <- sum(stringi::stri_count_fixed(str = fparam[idel], pattern = c("*")))
                if (length(cpm) == 0) {
                  cpm <- 0
                }

                if (cpm > 0) {
                  t.d <- sample(size = 1, x = (cpm))
                  # print(fparam[idel])

                  loc <- c(1, stri_locate_all(str = fparam[idel], regex = "\\*")[[1]][, 1], stringi::stri_length(fparam[idel]))
                  proposal <- stri_paste(stringi::stri_sub(fparam[idel], from = 1, to = loc[t.d] - 1 + 2 * (t.d == 1)), stringi::stri_sub(fparam[idel], from = (loc[t.d + 1] + (t.d == 1)), to = stringi::stri_length(fparam[idel])))

                  if (stats::runif(n = 1, min = 0, max = 1) < del.sigma) {
                    dsigmas <- sample(size = 1, x = length(sigmas))
                    if (sigmas[dsigmas] != "") {
                      proposal <- stringi::stri_replace_all_fixed(replacement = "", str = proposal, pattern = sigmas[dsigmas])
                    }
                  }
                  so <- stringi::stri_count_fixed(str = proposal, pattern = "(")
                  sc <- stringi::stri_count_fixed(str = proposal, pattern = ")")
                  # print(proposal)
                  if (sc > so) {
                    proposal <- stri_paste(stri_paste("I", rep("(", sc - so), collapse = ""), proposal)
                  } else if (sc < so) {
                    proposal <- stri_paste(proposal, stri_paste(rep(")", so - sc), collapse = ""))
                  }
                  # print(proposal)
                } else {
                  proposal <- fparam[idel]
                }
              }
            }


            sj <- (stringi::stri_count_fixed(str = proposal, pattern = "*"))
            sj <- sj + (stringi::stri_count_fixed(str = proposal, pattern = "+"))
            sj <- sj + sum(stringi::stri_count_fixed(str = proposal, pattern = sigmas))
            sj <- sj + 1
            if (length(sj) == 0) {
              proposal <- fparam.pool[sample(x = length(fparam.pool), size = 1)]
              sj <- 0
            } else if (is.na(proposal) || is.na(sj)) {
              proposal <- fparam.pool[sample(x = length(fparam.pool), size = 1)]
              sj <- 0
            } else if (sj > max.tree.size) {
              proposal <- fparam.pool[sample(x = length(fparam.pool), size = 1)]
            }


            add <- T


            ids.lat <- integer(0)
            if (latnames[1] != "") {
              ids.lat <- which(fparam %in% latnames)
              if (sum(stringi::stri_count_fixed(str = proposal, pattern = latnames)) > 0) {
                add <- F
              }
            }
            if (length(ids.lat) == 0) {
              ids.lat <- Nvars.max + 1
            }


            # print(add)
            # print(proposal)
            if (proposal %in% latnames) {
              if ((proposal %in% fparam[ids.lat])) {
                add <- F
              }
            } else {
              if (add) {
                withRestarts(tryCatch(utils::capture.output(
                  {
                    bet.act <- do.call(.self$estimator, c(estimator.args, stats::as.formula(stri_paste(fobserved, "~ 1 +", paste0(c(fparam[-ids.lat], proposal), collapse = "+")))))$summary.fixed$mean

                    if (is.na(bet.act[length(fparam[-ids.lat]) + 2]) && (action.type != 4 && gen.prob[2] == 0 || gen.prob[2] != 0)) {
                      add <- F
                    } else {
                      idel <- idel + 1
                    }
                  },
                  error = function(err) {
                    add <- F
                  },
                  finally = {}
                )))
              }
            }

            # print(add)
            # print("next")

            if (add & Nvars < Nvars.max) # alternative restricted to correlation: if((max(cor(eval(parse(text = proposal),envir = data.example),sapply(fparam, function(x) eval(parse(text=x),envir = data.example))))<0.9999) && Nvars<Nvars.max)
              {
                fparam <<- c(fparam, proposal)
                Nvars <<- as.integer(Nvars + 1)
                p.add <<- as.array(c(p.add, p.allow.replace))
                p.post <- as.array(c(p.post, 1))



                if (printable.opt) {
                  print(paste("mutation happended ", proposal, " tree  added"))
                }
              } else if (add) # alternative restricted to correlation: if(max(abs(cor(eval(parse(text = proposal),envir = data.example),sapply(fparam, function(x) eval(parse(text=x),envir = data.example)))))<0.9999)
              {

                # if(action.type==4)
                #   print(proposal)

                if (keep.origin) {
                  to.del <- (which(p.add[(Nvars.init + 1):Nvars] < p.allow.replace) + Nvars.init)
                  lto.del <- length(x = to.del)
                } else {
                  to.del <- which(p.add < p.allow.replace)
                  lto.del <- length(x = to.del)
                }
                # to.del<-(which(p.add[(Nvars.init+1):Nvars]< p.allow.replace)+ Nvars.init)
                # lto.del<-length(x = to.del)
                if (lto.del > 0) {
                  id.replace <- to.del[round(stats::runif(n = 1, min = 1, max = lto.del))]
                  if (printable.opt) {
                    print(paste("mutation happended ", proposal, " tree  replaced ", fparam[id.replace]))
                  }
                  fparam[id.replace] <<- proposal
                  keysarr <- as.array(hash::keys(hashStat))
                  p.add[id.replace] <<- p.allow.replace
                  for (jjj in 1:length(keysarr))
                  {
                    if (length(keysarr > 0)) {
                      if (stringi::stri_sub(keysarr[jjj], from = id.replace, to = id.replace) == "1") {
                        del(x = keysarr[jjj], hash = hashStat)
                      }
                    }
                  }
                }
              }
          }
        }

        if (p.add.default < 1) {
          p.add <<- array(p.add.default, Nvars)
        }
        varcurb <- c(varcurb, array(1, dim = (Nvars - length(varcurb))))
        varcand <- c(varcand, array(1, dim = (Nvars - length(varcand))))
        varglob <- c(varglob, array(1, dim = (Nvars - length(varglob))))
        p.post <- array(1, dim = (Nvars))
        p1 <- c(p1, array(0, dim = (Nvars - length(p1))))
        p2 <- c(p1, array(1, dim = (Nvars - length(p1))))
        acc_moves <- 1
        j.a <- 1
      } else if (allow_offsprings == 4 && j %% mutation_rate == 0 && (j <= last.mutation || Nvars != Nvars.max)) {
        add.buf <- F

        # perform preliminary filtration here
        if (Nvars > Nvars.max || j == mutation_rate) {
          on.suggested <- 1
          preaccepted <- F
          # do the stuff here
          if (j == mutation_rate) {
            fparam.pool <<- unique(c(fparam.pool, filtered))
            if (!pool.cor.prob) {
              pool.probs <- array(data = 1 / length(fparam.pool), dim = length(fparam.pool))
            } else {
              fobserved.cleaned <- fobserved
              fobserved.cleaned <- stri_replace(str = fobserved.cleaned, fixed = "I(", replacement = "")
              fobserved.cleaned <- stri_replace(str = fobserved.cleaned, fixed = ")", replacement = "")
              fparam.pool.cleaned <- fparam.pool
              fparam.pool.cleaned <- stri_replace(str = fparam.pool.cleaned, fixed = "I(", replacement = "")
              fparam.pool.cleaned <- stri_replace(str = fparam.pool.cleaned, fixed = ")", replacement = "")
              pool.probs <- abs(cor(estimator.args$data[[fobserved.cleaned]], estimator.args$data[, which(fparam.pool.cleaned %in% names(estimator.args$data))])) + p.epsilon
              rm(fobserved.cleaned)
              rm(fparam.pool.cleaned)
              # print(pool.probs[1:100])
            }
          }
          to.del <- which(p.add < p.allow.tree)
          if (length(to.del) == Nvars) {
            to.del <- to.del[-sample(x = Nvars, size = sample(x = Nvars - 1, size = 1), prob = p.add + p.epsilon)]
          }
          if (length(to.del) < Nvars - Nvars.max) {
            tdl.id <- order(p.add, decreasing = T)
            to.del <- to.del[-tdl.id[1:Nvars.max]]
          }
          if (glob.model$print.freq > 0L) {
            message("Data filtered! Insignificant variables deleted!")
          }
          if (length(to.del) > 0) {
            hash::clear(hashStat)
            fparam <<- fparam[-to.del]
            if (!keep.origin) {
              pool.probs[which(fparam.pool %in% fparam)] <- 1
            }
            Nvars <<- length(fparam)
            Nvars.init <<- Nvars
            p.add <<- p.add[-to.del]
            p.post <- array(data = 1, dim = Nvars)
            varcurb <- varcurb[-to.del]
            varcand <- varcurb[-to.del]
            varglob <- varcurb[-to.del]
            p1 <- array(0, dim = (Nvars))
            p2 <- array(1, dim = (Nvars))
            acc_moves <- 1
            j.a <- 1
            super.backward <- F
          }
        } else {
          if (Nvars >= Nvars.max) {
            # delete those that are not in the active model with probability 0.5 each

            if (keep.origin) {
              idmut <- (which(p.add[(Nvars.init + 1):Nvars] <= p.allow.replace) + Nvars.init)
              lidmut <- length(idmut) # maximal number of covariates that can die out
              if (lidmut > 0) {
                p.del <- (lidmut - sum(p.add[idmut])) / lidmut
                lidmut <- stats::rbinom(n = 1, size = lidmut, prob = p.del)
              }
            } else {
              idmut <- which(p.add <= p.allow.replace)
              lidmut <- length(idmut) # maximal number of covariates that can die out
              if (lidmut > 0) {
                p.del <- (lidmut - sum(p.add[idmut])) / lidmut
                lidmut <- stats::rbinom(n = 1, size = lidmut, prob = p.del)
              }
            }
            # mod.id.old<-stats::runif(n = 1000,min = 1,max = 2^Nvars)
            if (on.suggested %% 2 == 1) {
              if (preaccepted) {
                log.mod.switch.prob.back <- sum(abs(varcurb.buf - varcurb))
                eqal.things <- which(varcurb.buf == varcurb)
                log.mod.switch.prob.back <- log.mod.switch.prob.back + length(which(fparam[eqal.things] != tmp.buf.1[eqal.things]))
                log.mod.switch.prob.back <- log.mod.switch.prob.back^prand
                if (printable.opt) print(paste0("afteracceptance ratio ", log.mod.switch.prob.back / log.mod.switch.prob))
                if (stats::runif(1, 0, 1) <= log.mod.switch.prob.back / log.mod.switch.prob) {
                  if (printable.opt) print("preaccepted mutation is accepted")
                  mlikcur <- mlikcur.buf.2
                  varcurb <- varcurb.buf.2
                  fparam <<- tmp.buf.2
                  p.add <<- p.add.buf.2
                } else {
                  if (printable.opt) print("preaccepted mutation rejected")
                  fparam <<- tmp.buf.1
                  p.add <<- p.add.buf.1
                  mlikcur <- mlikcur.buf
                  varcurb <- varcurb.buf
                }
                preaccepted <- F
              }

              tmp.buf <- fparam[which(varcurb == 1)]
              tmp.buf.1 <- fparam
              p.add.buf.1 <- p.add

              # hashStat.buf.1<-copy(hashStat)
              vect <- buildmodel(max.cpu = 1, varcur.old = varcurb, statid = -1, min.N = Nvars, max.N = Nvars, switch.type = 9)
              res.par <- lapply(X = vect, FUN = .self$fitmodel)

              mlikcur <- res.par[[1]]$mlik
              varcurb <- vect[[1]]$varcur

              mlikcur.buf <- mlikcur
              varcurb.buf <- varcurb
              if (printable.opt) print("OLDMOD")
              mso <- (mlikcur) # (sum(unlist(lapply(res.par, function(x) exp(x$mlik)))))
              if (printable.opt) print(mso)
            } else if (on.suggested %% 2 == 0) {

              # tmp.buf2<-fparam[which(varcurb==1)]
              tmp.buf.2 <- fparam
              p.add.buf.2 <- p.add
              vect <- buildmodel(max.cpu = 1, varcur.old = varcurb, statid = -1, min.N = Nvars, max.N = Nvars, switch.type = 9)
              log.mod.switch.prob <- vect[[1]]$log.mod.switch.prob
              res.par <- lapply(X = vect, FUN = .self$fitmodel)
              # print(mlikcur)

              mlikcur <- res.par[[1]]$mlik
              varcurb <- vect[[1]]$varcur

              vect <- buildmodel(max.cpu = 1, varcur.old = varcurb, statid = -1, min.N = Nvars, max.N = Nvars, switch.type = 9)
              log.mod.switch.prob <- vect[[1]]$log.mod.switch.prob
              res.par <- lapply(X = vect, FUN = .self$fitmodel)
              # print(mlikcur)

              mlikcur <- res.par[[1]]$mlik
              varcurb <- vect[[1]]$varcur

              mlikcur.buf.2 <- mlikcur
              varcurb.buf.2 <- varcurb
              msn <- (mlikcur) # (sum(unlist(lapply(res.par, function(x) exp(x$mlik)))))
              if (printable.opt) print("NEWMOD")
              if (printable.opt) print(msn)
              if (msn == Inf && mso == Inf || msn == 0 && mso == 0) {
                msn <- 1
                mso <- 1
              } #* log((sum(tmp.buf%in%fparam)==length(tmp.buf)))
              if (log(stats::runif(1, 0, 1)) >= ((msn) - (mso))) {
                if (printable.opt) print("proposal is rejected")
                fparam <<- tmp.buf.1
                p.add <<- p.add.buf.1
                mlikcur <- mlikcur.buf
                varcurb <- varcurb.buf
                preaccepted <- F
                # on.suggested<-on.suggested+1
              } else {
                if (printable.opt) print("mutation is preaccepted")
                preaccepted <- T
              }
            }
            on.suggested <- on.suggested + 1
          } else {
            idmut <- (Nvars + 1):Nvars.max
            lidmut <- Nvars.max - Nvars
          }
          # now having chosen the candidates to be deleted we can propose new variables



          idel <- 1

          while (idel <= lidmut) {

            # gen.prob<-c(1,1,1,1,1)#just uniform for now
            action.type <- sample(x = 5, size = 1, prob = gen.prob)

            if (action.type == 1) {
              # mutation (add a leave not in the search space)
              proposal <- fparam.pool[sample(x = length(fparam.pool), size = 1, prob = pool.probs)]
            } else if (action.type == 2) {
              # crossover type of a proposal
              # generate a mother
              # actvars<-which(varcurb==1)
              mother <- ifelse(stats::runif(n = 1, min = 0, max = 1) <= pool.cross, fparam[sample(x = Nvars, size = 1, prob = p.add + p.epsilon)], fparam.pool[sample(x = length(fparam.pool), size = 1, prob = pool.probs)])
              ltreem <- stringi::stri_length(mother)
              mother <- stringi::stri_sub(mother, from = 1, to = ltreem)
              # sjm<-sum(stringi::stri_count_fixed(str = mother, pattern = c("+","*")))
              # generate a father
              father <- ifelse(stats::runif(n = 1, min = 0, max = 1) <= pool.cross, fparam[sample(x = Nvars, size = 1, prob = p.add + p.epsilon)], fparam.pool[sample(x = length(fparam.pool), size = 1, prob = pool.probs)])
              ltreef <- stringi::stri_length(father)
              father <- stringi::stri_sub(father, from = 1, to = ltreef)
              # sjf<-sum(stringi::stri_count_fixed(str = father, pattern = c("+","*")))

              proposal <- stri_paste("I(", stri_paste(mother, father, sep = "*"), ")", sep = "")
            } else if (action.type == 3) {
              proposal <- ifelse(stats::runif(n = 1, min = 0, max = 1) <= pool.cross, fparam[sample(x = Nvars, size = 1, prob = p.add + p.epsilon)], fparam.pool[sample(x = length(fparam.pool), size = 1, prob = pool.probs)])
              proposal <- stri_paste("I(", sigmas[sample(x = length(sigmas), size = 1, prob = sigmas.prob)], "(", proposal, "))", sep = "")
            } else if (action.type == 4) {


              # select a sparse subset for the projection

              spad <- sum(p.add)
              actvars <- which(stats::rbinom(n = length(fparam), size = 1, prob = p.add / spad + p.epsilon) == 1)


              if (length(actvars) <= 1) {
                proposal <- fparam[1]
              } else {
                utils::capture.output(withRestarts(tryCatch(
                  {
                    # print(stats::as.formula(stri_paste(fobserved,"~ 1 +",paste0(fparam[actvars],collapse = "+"))))
                    # get the projection coefficients as the posterior mode of the fixed effects
                    if (deep.method == 1) {
                      bet.act <- do.call(.self$estimator, c(estimator.args, stats::as.formula(stri_paste(fobserved, "~ 1 +", paste0(fparam[actvars], collapse = "+")))))$summary.fixed$mean
                      nab <- which(is.na(bet.act))
                      if (length(nab) > 0) {
                        bet.act[nab] <- 0
                      }
                      bet.act <- round(bet.act, digits = 8)
                      bet.act <- stri_paste("m(", bet.act, ",", sep = "")
                      # make a projection
                      bet.act <- stri_paste(bet.act, c("1", fparam[actvars]), ")", sep = "")
                      bet.act <- stri_paste("I(", bet.act, ")", sep = "")
                      proposal <- stri_paste("I(", stri_paste(bet.act, collapse = "+"), ")", collapse = "")
                      proposal <- stri_paste("I(", sigmas[sample(x = length(sigmas), size = 1, prob = sigmas.prob)], "(", proposal, "))", sep = "")
                      # print(proposal)
                    } else if (deep.method == 2) {
                      cursigma <- sigmas[sample(x = length(sigmas), size = 1, prob = sigmas.prob)]
                      bet.act <- gnlr(
                        y = data.example[[fobserved]],
                        distribution = estimator.args$distribution,
                        mu = stats::as.formula(stri_paste("~", estimator.args$link, "(", cursigma, "(", "b0 +", paste0("b", 1:length(actvars), "*", fparam[actvars], collapse = "+"), "))")),
                        pmu = rep(0, length(actvars) + 1)
                      )$coefficients


                      nab <- which(is.na(bet.act))
                      if (length(nab) > 0) {
                        bet.act[nab] <- 0
                      }
                      # bet.act<-round(bet.act, digits = 2)
                      bet.act <- round(bet.act, digits = 8)
                      bet.act <- stri_paste("m(", bet.act, ",", sep = "")
                      # make a projection
                      bet.act <- stri_paste(bet.act, c("1", fparam[actvars]), ")", sep = "")
                      bet.act <- stri_paste("I(", bet.act, ")", sep = "")
                      proposal <- stri_paste("I(", stri_paste(bet.act, collapse = "+"), ")", collapse = "")
                      proposal <- stri_paste("I(", cursigma, "(", proposal, "))", sep = "")
                    } else if (deep.method == 3) {
                      cursigma <- sigmas[sample(x = length(sigmas), size = 1, prob = sigmas.prob)]
                      forstr <- stri_paste("~", estimator.args$link, "(", cursigma, "(", "m(0,1) +", paste0("m(", stats::rnorm(length(actvars), 0, 1), ",", fparam[actvars], ")", collapse = "+"), "))")
                      forstr <- stringi::stri_replace_all_fixed(str = forstr, pattern = c("m(-"), replacement = "m(")
                      forstr <- stringi::stri_replace_all_fixed(str = forstr, pattern = c("m("), replacement = "m(b_")
                      # print(stats::as.formula(forstr))

                      nlrr <- gnlr(
                        y = data.example[[fobserved]],
                        distribution = estimator.args$distribution,
                        mu = stats::as.formula(forstr),
                        pmu = stats::rnorm(stringi::stri_count_fixed(str = forstr, pattern = "m("), 0, 0.0001)
                      ) # $coefficients
                      beg.rep <- stri_locate_all(str = forstr, fixed = "m(b")[[1]][, 2]
                      end.rep <- stri_locate_all(str = forstr, fixed = ",")[[1]][, 2]

                      bet.act <- nlrr$coefficients

                      nab <- which(is.na(bet.act))
                      if (length(nab) > 0) {
                        bet.act[nab] <- 0
                      }
                      bet.act <- round(bet.act, digits = 8)
                      torepl <- stringi::stri_sub(str = forstr, from = beg.rep, to = end.rep)
                      for (i in 1:length(bet.act)) {
                        forstr <- stringi::stri_replace_all_fixed(str = forstr, pattern = torepl[i], replacement = stri_paste(bet.act[i], ","))
                      }
                      forstr <- stringi::stri_replace_all_fixed(str = forstr, pattern = paste0("~", estimator.args$link), replacement = "I")

                      proposal <- forstr
                    } else {
                      bet.act <- stats::rnorm(n = (length(actvars) + 1), mean = 0, sd = 1)
                      nab <- which(is.na(bet.act))
                      if (length(nab) > 0) {
                        bet.act[nab] <- 0
                      }
                      bet.act <- round(bet.act, digits = 8)
                      bet.act <- stri_paste("m(", bet.act, ",", sep = "")
                      # make a projection
                      bet.act <- stri_paste(bet.act, c("1", fparam[actvars]), ")", sep = "")
                      bet.act <- stri_paste("I(", bet.act, ")", sep = "")
                      proposal <- stri_paste("I(", stri_paste(bet.act, collapse = "+"), ")", collapse = "")
                      proposal <- stri_paste("I(", sigmas[sample(x = length(sigmas), size = 1, prob = sigmas.prob)], "(", proposal, "))", sep = "")
                    }
                  },
                  error = function(err) {
                    # print(err)
                    proposal <- fparam.pool[sample(x = length(fparam.pool), size = 1)]
                  },
                  finally = {
                    proposal <- proposal
                  }
                )))
                # print(proposal)

                if (is.na(proposal)) {
                  print(fparam[actvars])
                  print(actvars)
                  print(bet.act)
                }
              }
            } else if (action.type == 5) {
              # reduce an operator fparam[idel]
              # print("reduction")
              # print(idel)
              if (idel > length(fparam)) {
                proposal <- fparam[1]
              } else {
                cpm <- sum(stringi::stri_count_fixed(str = fparam[idel], pattern = c("*")))
                if (length(cpm) == 0) {
                  cpm <- 0
                }

                if (cpm > 0) {
                  t.d <- sample(size = 1, x = (cpm))
                  # print(fparam[idel])

                  loc <- c(1, stri_locate_all(str = fparam[idel], regex = "\\*")[[1]][, 1], stringi::stri_length(fparam[idel]))
                  proposal <- stri_paste(stringi::stri_sub(fparam[idel], from = 1, to = loc[t.d] - 1 + 2 * (t.d == 1)), stringi::stri_sub(fparam[idel], from = (loc[t.d + 1] + (t.d == 1)), to = stringi::stri_length(fparam[idel])))

                  if (stats::runif(n = 1, min = 0, max = 1) < del.sigma) {
                    dsigmas <- sample(size = 1, x = length(sigmas))
                    if (sigmas[dsigmas] != "") {
                      proposal <- stringi::stri_replace_all_fixed(replacement = "", str = proposal, pattern = sigmas[dsigmas])
                    }
                  }
                  so <- stringi::stri_count_fixed(str = proposal, pattern = "(")
                  sc <- stringi::stri_count_fixed(str = proposal, pattern = ")")
                  # print(proposal)
                  if (sc > so) {
                    proposal <- stri_paste(stri_paste("I", rep("(", sc - so), collapse = ""), proposal)
                  } else if (sc < so) {
                    proposal <- stri_paste(proposal, stri_paste(rep(")", so - sc), collapse = ""))
                  }
                  # print(proposal)
                } else {
                  proposal <- fparam[idel]
                }
              }
            }


            sj <- (stringi::stri_count_fixed(str = proposal, pattern = "*"))
            sj <- sj + (stringi::stri_count_fixed(str = proposal, pattern = "+"))
            sj <- sj + sum(stringi::stri_count_fixed(str = proposal, pattern = sigmas))
            sj <- sj + 1
            if (sj > max.tree.size || length(sj) == 0) {
              proposal <- fparam.pool[sample(x = length(fparam.pool), size = 1)]
            }
            if (is.na(proposal)) {
              proposal <- fparam.pool[sample(x = length(fparam.pool), size = 1)]
            }

            add <- T

            tryCatch(utils::capture.output(
              {
                bet.act <- do.call(.self$estimator, c(estimator.args, stats::as.formula(stri_paste(fobserved, "~ 1 +", paste0(c(fparam, proposal), collapse = "+")))))$summary.fixed$mean

                if (is.na(bet.act[length(fparam) + 2])) {
                  add <- F
                } else {
                  idel <- idel + 1
                }
              },
              error = function(err) {
                add <- F
              }
            ))

            if (add & Nvars < Nvars.max) # alternative restricted to correlation: if((max(cor(eval(parse(text = proposal),envir = data.example),sapply(fparam, function(x) eval(parse(text=x),envir = data.example))))<0.9999) && Nvars<Nvars.max)
              {
                fparam <<- c(fparam, proposal)
                Nvars <<- as.integer(Nvars + 1)
                p.add <<- as.array(c(p.add, p.allow.replace))
                p.post <- as.array(c(p.post, 1))
                # if(printable.opt)
                if (printable.opt) print(paste("mutation happended ", proposal, " tree  added"))
              } else if (add) # alternative restricted to correlation: if(max(abs(cor(eval(parse(text = proposal),envir = data.example),sapply(fparam, function(x) eval(parse(text=x),envir = data.example)))))<0.9999)
              {
                if (keep.origin) {
                  to.del <- (which(p.add[(Nvars.init + 1):Nvars] < p.allow.replace) + Nvars.init)
                  lto.del <- length(x = to.del)
                } else {
                  to.del <- which(p.add < p.allow.replace)
                  lto.del <- length(x = to.del)
                }

                if (lto.del > 0) {
                  id.replace <- to.del[round(stats::runif(n = 1, min = 1, max = lto.del))]
                  # if(printable.opt)
                  # print(paste("mutation suggested ",proposal," tree  replaced ", fparam[id.replace]))
                  fparam[id.replace] <<- proposal
                  keysarr <- as.array(hash::keys(hashStat))
                  p.add[id.replace] <<- p.allow.replace
                  for (jjj in 1:length(keysarr))
                  {
                    if (length(keysarr > 0)) {
                      if (stringi::stri_sub(keysarr[jjj], from = id.replace, to = id.replace) == "1") {
                        del(x = keysarr[jjj], hash = hashStat)
                      }
                    }
                  }
                }
              }
          }
        }
        varcurb <- c(varcurb, array(1, dim = (Nvars - length(varcurb))))
        varcand <- c(varcand, array(1, dim = (Nvars - length(varcand))))
        varglob <- c(varglob, array(1, dim = (Nvars - length(varglob))))
        p.post <- array(1, dim = (Nvars))
        p1 <- c(p1, array(0, dim = (Nvars - length(p1))))
        p2 <- c(p1, array(1, dim = (Nvars - length(p1))))
        acc_moves <- 1
        j.a <- 1
      } else if (allow_offsprings > 0 && j %% mutation_rate == 0 && j > last.mutation) {
        recalc.margin <<- 2^Nvars
      }
      # withRestarts(tryCatch({

      varcur <- varcurb



      if (LocImprove <= as.array(3)) {
        vect <- buildmodel(max.cpu = 1, varcur.old = varcurb, statid = 4 + LocImprove, min.N = min.N.glob, max.N = max.N.glob, switch.type = switch.type.glob.buf)
        max.cpu.buf <- 1
      } else {
        vect <- buildmodel(max.cpu = max.cpu.glob, varcur.old = varcurb, statid = 4 + LocImprove, min.N = min.N, max.N = max.N, switch.type = switch.type.glob.buf)
        max.cpu.buf <- max.cpu.glob
      }

      cluster <- TRUE

      flag1 <- 0

      for (mod_id in 1:max.cpu.buf)
      {
        if (is.null(vect[[mod_id]]$formula)) {
          flag1 <- flag1 + 1
        }
      }

      if (flag1 == max.cpu.glob) {
        cluster <- FALSE
        if (printable.opt) print("!!!!Models already estimated!!!!")
      } else {
        if (max.cpu.glob > 1) {
          res.par <- parallelize(X = vect, FUN = .self$fitmodel)
        } else {
          res.par <- lapply(X = vect, FUN = .self$fitmodel)
        }
      }

      if (LocImprove > 3) {
        if (printable.opt) print("!!!!Proceed with no local improvements!!!!")
        p.select.y <- array(data = 0, dim = max.cpu.glob)
        for (mod_id in 1:max.cpu.glob)
        {
          if (cluster) {
            fm <- res.par[[mod_id]]

            if (is.null(fm) && (is.na(res.par[[mod_id]]$waic))) {
              varcand <- varcurb
              if (printable.opt) print("GlobMTMCMC Model Fit Error!?")
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
          } else if (exists("hashStat") & length(hashStat) > 0) {
            iidd <- paste(varcand, collapse = "")
            waiccand <- hash::values(hashStat[iidd])[2]
            mlikcand <- hash::values(hashStat[iidd])[1]
          }

          if ((mlikcand > mlikglob)) # update the parameter of interest
            {
              if (printable.opt) print(paste("GlobMTMCMC update waic.glob = ", waiccand))
              if (printable.opt) print(paste("GlobMTMCMC update mlik.glob = ", mlikglob))
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


        ID <- sample(x = max.cpu.glob, size = 1, prob = exp(p.select.y))

        if (printable.opt) print(paste("cand ", ID, " selected"))

        varcand <- vect[[ID]]$varcur

        if (cluster) {
          waiccand <- res.par[[ID]]$waic
          mlikcand <- res.par[[ID]]$mlik
        } else if (exists("statistics1")) {
          iidd <- bittodec(varcand) + 1
          waiccand <- statistics1[iidd, 2]
          mlikcand <- statistics1[iidd, 1]
        } else if (exists("hashStat") & length(hashStat) > 0) {
          iidd <- paste(varcand, collapse = "")
          waiccand <- hash::values(hashStat[iidd])[2]
          mlikcand <- hash::values(hashStat[iidd])[1]
        }

        # p.Q.cand<- p.select.y[ID]/sum(p.select.y)

        if (printable.opt) print("do reverse step")

        p.select.z <- array(data = 0.01, dim = max.cpu.glob)


        if (max.cpu.glob != 1) {
          if (switch.type.glob.buf == 5) {
            cstm <- 6
          } else if (switch.type.glob.buf == 6) {
            cstm <- 5
          } else {
            cstm <- switch.type.glob.buf
          }
          vect1 <- buildmodel(max.cpu = max.cpu.glob - 1, varcur.old = varcand, statid = 4 + LocImprove, switch.type = cstm, min.N = min.N, max.N = max.N)

          cluster <- TRUE

          flag1 <- 0

          for (mod_id in 1:(max.cpu.glob - 1))
          {
            if (is.null(vect1[[mod_id]]$formula)) {
              flag1 <- flag1 + 1
            }
          }

          if (flag1 == (max.cpu.glob - 1)) {
            cluster <- FALSE
            if (printable.opt) print("!!!!MTMCMC reverse models already estimated!!!!")
          } else {
            res.par.back <- parallelize(X = vect1, FUN = .self$fitmodel)
          }

          for (mod_id in 1:(max.cpu.glob - 1))
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
            } else if (exists("hashStat") & length(hashStat) > 0) {
              iidd <- paste(varcand.b, collapse = "")
              waiccand.b <- hash::values(hashStat[iidd])[2]
              mlikcand.b <- hash::values(hashStat[iidd])[1]
            }

            if ((mlikcand.b > mlikglob)) {
              if (printable.opt) print(paste("GlobMTMCMC update waic.glob = ", waiccand.b))
              if (printable.opt) print(paste("GlobMTMCMC update mlik.glob = ", mlikcand.b))
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
        p.select.z[max.cpu.glob] <- (mlikcur + vect[[ID]]$log.mod.switch.prob + (lambda(c = cc, alpha = aa, g1 = -g1, g2 = -waiccand, g.domain.pos = FALSE)))

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
          ratcur <- mlikcand
          if (printable.opt) print(paste("global MTMCMC update ratcur = ", mlikcand))
          if (printable.opt) print(paste("global MTMCMC accept move with ", waiccand))
          varcurb <- varcand
          waiccur <- waiccand

          id <- bittodec(varcurb) + 1

          acc_moves <- acc_moves + 1
          if (j < glob.model$burnin) {
            distrib_of_proposals[5] <- distrib_of_proposals[5] + 1
          } else {
            if (exists("statistics1")) {
              statistics1[id, 14] <- statistics1[id, 14] + 1
              statistics1[id, 4] <- statistics1[id, 4] + 1
            }
            p.post <- (p.post + varcurb)
          }
        } else if (j < glob.model$burnin && distrib_of_proposals[5] > 0) {
          distrib_of_proposals[5] <- distrib_of_proposals[5] - 1
        }
      } else {
        if (cluster) {
          fm <- res.par[[mod_id]]

          if (is.null(fm) && (is.na(res.par[[mod_id]]$waic))) {
            varcand <- varcurb
            if (printable.opt) print("EMJMCMC Model Fit Error!?")
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

        varcur <- varcand
      }

      # try local improvements
      if (LocImprove <= as.array(3)) {
        # sa improvements
        if (LocImprove == as.array(0) || LocImprove == as.array(3)) {
          if (printable.opt) print("Try SA imptovements")
          buf.change <- array(data = 1, dim = Nvars)
          buf.change[which(varcur - varcurb != 0)] <- 0
          # buf.opt[floor(stats::runif(n = n.size,min = 1,max = Nvars+0.999999))] = 1

          if (objective == 0) {
            objold <- waiccur
            objcur <- waiccand
          } else {
            objold <- -mlikcur
            objcur <- -mlikcand
          }
          model <- list(
            statid = 4 + LocImprove, switch.type = switch.type.buf, change = buf.change, mlikcur = mlikcur,
            varcur = varcur, varold = varcurb, objcur = objcur, objold = objold, sa2 = ifelse(LocImprove == as.array(3), TRUE, FALSE), reverse = FALSE
          )

          SA.forw <- learnlocalSA(model)
          ratcand <- SA.forw$mlikcur

          if (LocImprove == as.array(0)) {
            ids <- which(buf.change == 1)
            model$varcur <- SA.forw$varcur
            model$waiccur <- SA.forw$waiccur
            model$varcur[ids] <- 1 - model$varcur[ids]

            model$mlikcur <- -Inf
            model$waiccur <- Inf
            # estimate the jump!!!
            if (objective == 0) {
              model$objold <- waiccur
            } else {
              model$objold <- -ratcand
            }

            model$reverse <- TRUE
            if (switch.type.buf == 5) {
              model$switch.type <- 6
            } else if (switch.type.buf == 6) {
              model$switch.type <- 5
            }
            SA.back <- learnlocalSA(model)
          }

          if (LocImprove == as.array(0)) {
            thact <- sum(ratcand, -ratcur, -SA.forw$log.prob.cur, SA.forw$log.prob.fix, SA.back$log.prob.cur, -SA.back$log.prob.fix, na.rm = T)
            if (log(stats::runif(n = 1, min = 0, max = 1)) <= thact) {
              ratcur <- ratcand
              mlikcur <- ratcand
              if (printable.opt) print(paste("update ratcur through SA = ", ratcur))
              varcurb <- SA.forw$varcur
              acc_moves <- acc_moves + 1
              id <- bittodec(varcurb) + 1

              if (j < glob.model$burnin) {
                distrib_of_proposals[1] <- distrib_of_proposals[1] + 1
              } else {
                if (exists("statistics1")) {
                  statistics1[id, 10] <- statistics1[id, 10] + 1
                  statistics1[id, 4] <- statistics1[id, 4] + 1
                }
                p.post <- (p.post + varcurb)
              }
            } else if (j < glob.model$burnin && distrib_of_proposals[1] > 0) {
              distrib_of_proposals[1] <- distrib_of_proposals[1] - 1
            }
            if ((SA.back$mlikglob < mlikglob)) {
              if (printable.opt) print(paste("update waic.glob sa.back = ", SA.back$waicglob))
              if (printable.opt) print(paste("update waic.glob.mlik sa.back= ", SA.back$mlikglob))
              waicglob <- SA.back$waicglob
              varglob <- SA.back$varglob
              modglob <- SA.back$modglob
            }
          } else {
            thact <- sum(ratcand, -ratcur, -SA.forw$log.prob.cur, SA.forw$log.prob.fix, vect[[mod_id]]$log.mod.switchback.prob, -vect[[mod_id]]$log.mod.switch.prob, na.rm = T)
            if (log(stats::runif(n = 1, min = 0, max = 1)) <= thact) {
              ratcur <- ratcand
              mlikcur <- ratcand
              if (printable.opt) print(paste("update ratcur through SA = ", ratcur))
              varcurb <- SA.forw$varcur

              id <- bittodec(varcurb) + 1

              acc_moves <- acc_moves + 1
              if (j < glob.model$burnin) {
                distrib_of_proposals[4] <- distrib_of_proposals[4] + 1
              } else {
                if (exists("statistics1")) {
                  statistics1[id, 13] <- statistics1[id, 13] + 1
                  statistics1[id, 4] <- statistics1[id, 4] + 1
                }
                p.post <- (p.post + varcurb)
              }
            } else if (j < glob.model$burnin && distrib_of_proposals[4] > 0) {
              distrib_of_proposals[4] <- distrib_of_proposals[4] - 1
            }
          }
          if ((SA.forw$mlikglob < mlikglob)) {
            if (printable.opt) print(paste("update waic.glob sa.forw = ", SA.forw$waicglob))
            if (printable.opt) print(paste("update waic.glob.mlik sa.forw = ", SA.forw$mlikglob))
            waicglob <- SA.forw$waicglob
            varglob <- SA.forw$varglob
            modglob <- SA.forw$modglob
          }
        } else if (LocImprove == as.array(1)) {
          if (printable.opt) print("Try MTMCMC imptovements")
          buf.change <- array(data = 1, dim = Nvars)
          buf.change[which(varcur - varcurb != 0)] <- 0
          # buf.opt[floor(stats::runif(n = n.size,min = 1,max = Nvars+0.999999))] = 1



          model <- list(
            statid = 4 + LocImprove, reverse = FALSE, change = buf.change, mlikcur = mlikcur, waiccur = waiccur,
            varcur = varcur, varold = varcurb
          )

          MTMCMC.forw <- learnlocalMCMC(model)
          ids <- which(buf.change == 1)
          model$varcur <- MTMCMC.forw$varcur


          model$varcur[ids] <- 1 - model$varcur[ids]
          model$mlikcur <- -Inf
          model$waiccur <- Inf
          model$reverse <- TRUE

          # if(printable.opt)print("learn reverse local MTMCMC")
          MTMCMC.back <- learnlocalMCMC(model)

          # if(printable.opt)print("finish reverse local MTMCMC")

          MTMCMC.p.forw <- MTMCMC.forw$log.prob.cur
          MTMCMC.p.back <- MTMCMC.back$log.prob.fix
          ratcand <- MTMCMC.forw$mlikcur


          # if(log(stats::runif(n = 1,min = 0,max = 1))<=(ratcand - ratcur - MTMCMC.forw$log.prob.cur + MTMCMC.forw$log.prob.fix + MTMCMC.back$log.prob.cur - MTMCMC.back$log.prob.fix))
          thact <- sum(ratcand, -ratcur, -MTMCMC.forw$log.prob.cur, MTMCMC.forw$log.prob.fix, MTMCMC.back$log.prob.cur, -MTMCMC.back$log.prob.fix, na.rm = T)
          if (log(stats::runif(n = 1, min = 0, max = 1)) <= thact) {
            ratcur <- ratcand
            mlikcur <- ratcand
            # if(printable.opt)print(paste("update ratcur through MTMCMC = ", ratcur))
            varcurb <- MTMCMC.forw$varcur
            acc_moves <- acc_moves + 1
            id <- bittodec(varcurb) + 1

            if (j < glob.model$burnin) {
              distrib_of_proposals[2] <- distrib_of_proposals[2] + 1
            } else {
              if (exists("statistics1")) {
                statistics1[id, 11] <- statistics1[id, 11] + 1
                statistics1[id, 4] <- statistics1[id, 4] + 1
              }
              p.post <- (p.post + varcurb)
            }
          } else if (j < glob.model$burnin && distrib_of_proposals[2] > 1) {
            distrib_of_proposals[2] <- distrib_of_proposals[2] - 1
          }

          if ((MTMCMC.forw$mlikglob < mlikglob)) {
            # if(printable.opt)print(paste("update waic.glob MTMCMC.forw = ", MTMCMC.forw$waicglob))
            # if(printable.opt)print(paste("update waic.glob.mlik MTMCMC.forw = ", MTMCMC.forw$mlikglob))
            waicglob <- MTMCMC.forw$waicglob
            varglob <- MTMCMC.forw$varglob
            modglob <- MTMCMC.forw$modglob
          }

          if ((MTMCMC.back$mlikglob < mlikglob)) {
            # if(printable.opt)print(paste("update waic.glob MTMCMC.back = ", MTMCMC.back$waicglob))
            # if(printable.opt)print(paste("update waic.glob.mlik MTMCMC.back= ", MTMCMC.back$mlikglob))
            waicglob <- MTMCMC.back$waicglob
            varglob <- MTMCMC.back$varglob
            modglob <- MTMCMC.back$modglob
          }
        } else if (LocImprove == as.array(2)) {
          if (printable.opt) print("Try greedy heuristic imptovements")
          buf.change <- array(data = 1, dim = Nvars)
          buf.change[which(varcur - varcurb != 0)] <- 0

          # buf.opt <- buf.change
          # buf.opt[floor(stats::runif(n = n.size,min = 1,max = Nvars+0.999999))] = 1
          if (objective == 0) {
            objold <- waiccur
            objcur <- waiccand
          } else {
            objold <- -mlikcur
            objcur <- -mlikcand
          }

          model <- list(
            statid = 4 + LocImprove, change = buf.change,
            varcur = varcur, varold = varcurb, switch.type = switch.type.buf, mlikcur = mlikcur, objcur = objcur, objold = objold, reverse = FALSE
          )
          GREEDY.forw <- learnlocalND(model)
          ids <- which(buf.change == 1)
          model$varcur <- GREEDY.forw$varcur

          model$varcur[ids] <- 1 - model$varcur[ids]
          model$mlikcur <- -Inf
          model$waiccur <- Inf
          model$reverse <- TRUE
          ratcand <- GREEDY.forw$mlikcur

          if (objective == 0) {
            model$objold <- waiccur
          } else {
            model$objold <- -ratcand
          }
          if (switch.type.buf == 5) {
            model$switch.type <- 6
          } else if (switch.type.buf == 6) {
            model$switch.type <- 5
          }

          GREEDY.back <- learnlocalND(model)

          GREEDY.p.forw <- GREEDY.forw$log.prob.cur
          GREEDY.p.back <- GREEDY.back$log.prob.fix



          # if(log(stats::runif(n = 1,min = 0,max = 1))<=(ratcand - ratcur - GREEDY.forw$log.prob.cur + GREEDY.forw$log.prob.fix + GREEDY.back$log.prob.cur - GREEDY.back$log.prob.fix))
          thact <- sum(ratcand, -ratcur, -GREEDY.forw$log.prob.cur, GREEDY.forw$log.prob.fix, GREEDY.back$log.prob.cur, -GREEDY.back$log.prob.fix, na.rm = T)
          if (log(stats::runif(n = 1, min = 0, max = 1)) <= thact) {
            ratcur <- ratcand
            mlikcur <- ratcand
            if (printable.opt) print(paste("update ratcur through ND = ", ratcur))
            varcurb <- GREEDY.forw$varcur
            acc_moves <- acc_moves + 1
            id <- bittodec(varcurb) + 1

            if (j < glob.model$burnin) {
              distrib_of_proposals[3] <- distrib_of_proposals[3] + 1
            } else {
              if (exists("statistics1")) {
                statistics1[id, 12] <- statistics1[id, 12] + 1
                statistics1[id, 4] <- statistics1[id, 4] + 1
              }
              p.post <- (p.post + varcurb)
            }
          } else if (j < glob.model$burnin && distrib_of_proposals[3] > 1) {
            distrib_of_proposals[3] <- distrib_of_proposals[3] - 1
          }

          if ((GREEDY.forw$mlikglob < mlikglob)) {
            if (printable.opt) print(paste("update waic.glob ND.forw = ", GREEDY.forw$waicglob))
            if (printable.opt) print(paste("update waic.glob.mlik ND.forw = ", GREEDY.forw$mlikglob))
            waicglob <- GREEDY.forw$waicglob
            varglob <- GREEDY.forw$varglob
            modglob <- GREEDY.forw$modglob
          }

          if ((GREEDY.back$mlikglob < mlikglob)) {
            if (printable.opt) print(paste("update waic.glob ND.back = ", GREEDY.back$waicglob))
            if (printable.opt) print(paste("update waic.glob.mlik ND.back= ", GREEDY.back$mlikglob))
            waicglob <- GREEDY.back$waicglob
            varglob <- GREEDY.back$varglob
            modglob <- GREEDY.back$modglob
          }
        }
      }

      # }),abort = function(){if(printable.opt)print("error");varcur<-varcurb;closeAllConnections();options(error=traceback);  onerr<-TRUE})

      if (thin_rate != -1) {
        if (acc_moves == accept_old && j > glob.model$burnin && j %% as.integer(thin_rate) == 0) # carry out smart thinning
          {
            if (!is.null(varcurb)) {
              p.post <- (p.post + varcurb)
              thin_count <- thin_count + 1
              id <- bittodec(varcurb) + 1
              if (exists("statistics1")) {
                statistics1[id, 4] <- statistics1[id, 4] + 1
              }
            }
          }
      }

      accept_old <- acc_moves

      p2 <- p.post / acc_moves
      if (j > glob.model$burnin && (recalc.margin >= 2^Nvars) && sum(p2) != 0) {
        p.add <<- as.array(p2)
      }
      eps.emp <- normprob(p1, p2)
      etm <- proc.time()
      delta.time <- (etm[3] - stm[3]) / 60.0
    }

    if (is.null(modglob)) {
      vect <- buildmodel(max.cpu = 1, varcur.old = varglob, statid = 3, switch.type = 8)
      res.par <- lapply(X = vect, FUN = fitmodel)
      modglob <- res.par[[1]]$fm
    }
    # !#print the results

    if (printable.opt) print(paste(j, " iterations completed in total"))
    # if(printable.opt)print(paste("WAIC.glob = ", waicglob))
    etm <- proc.time()
    tt <- (etm[3] - stm[3]) / 60.0

    acc_ratio <- acc_moves / j.a

    if (printable.opt) print(paste(j, " moves proposed in total, ", acc_moves, " of them accepted, acceptance ratio is ", acc_ratio))

    if (printable.opt) print(paste("posterior distribution ofproposals is", distrib_of_proposals))


    if (exists("statistics1")) {
      bayes.res <- post_proceed_results(statistics1)
      m.post <- statistics1[, 4] / j.a
    } else if (exists("hashStat")) {
      bayes.res <- post_proceed_results_hash(hashStat)
      m.post <- NULL
    } else {
      bayes.res <- NULL
      m.post <- NULL
    }

    return(list(model = modglob, vars = varglob, waic = waicglob, mlik = mlikglob, time = (tt), freqs = distrib_of_proposals, acc_ratio = acc_ratio, p.post = p.post / acc_moves, m.post = m.post, p.post.freq = p.post, eps = eps.emp, bayes.results = bayes.res))
  }
)
