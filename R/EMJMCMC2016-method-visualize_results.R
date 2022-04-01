EMJMCMC2016$methods(
  # vizualize the results graphically if the latter are available
  visualize_results = function(statistics1, template, mds_size, crit, draw_dist = FALSE) {
    ids <- NULL

    if (crit$mlik) {
      ids <- c(ids, which(statistics1[, 1] == -100000))
      ids <- c(ids, which(statistics1[, 1] == -Inf))
    }
    # ids<-c(ids,which(statistics1[,1]>0))
    if (crit$waic) {
      ids <- c(ids, which(statistics1[, 2] == Inf))
      ids <- c(ids, which(statistics1[, 3] == Inf))
    }
    if (crit$dic) {
      ids <- c(ids, which(is.na(statistics1[, 3])))
    }

    if (length(ids) != 0) {
      if (crit$mlik) {
        mlik.lim <- c(min(statistics1[-ids, 1], na.rm = TRUE), max(statistics1[-ids, 1], na.rm = TRUE))
      }
      if (crit$waic) {
        waic.lim <- c(min(statistics1[-ids, 2], na.rm = TRUE), max(statistics1[-ids, 2], na.rm = TRUE))
      }
      if (crit$dic) {
        dic.lim <- c(min(statistics1[-ids, 3], na.rm = TRUE), max(statistics1[-ids, 3], na.rm = TRUE))
      }
    } else {
      if (crit$mlik) {
        mlik.lim <- c(min(statistics1[, 1], na.rm = TRUE), max(statistics1[, 1], na.rm = TRUE))
      }
      if (crit$waic) {
        waic.lim <- c(min(statistics1[, 2], na.rm = TRUE), max(statistics1[, 2], na.rm = TRUE))
      }
      if (crit$dic) {
        dic.lim <- c(min(statistics1[, 3], na.rm = TRUE), max(statistics1[, 3], na.rm = TRUE))
      }
    }

    norm <- 0.5 * sqrt(sum(statistics1[, 4], na.rm = TRUE))

    if (crit$mlik) {
      if (printable.opt) print(paste("drawing ", workdir, template, "_legend.jpg", sep = ""))
      utils::capture.output({
        withRestarts(tryCatch(utils::capture.output({
          jpeg(file = paste(workdir, template, "_legend.jpg", sep = ""))
          plot(xlab = "posterior probability of a visit found", ylab = "visit type", c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), c(7, 7, 7, -1, -1, -1, -1), ylim = c(0, 9), xlim = c(0, 1), pch = 19, col = 7, cex = c(2, 2, 2, 0, 0, 0, 0))
          points(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), c(6, 6, 6, -1, -1, -1, -1), pch = 8, col = 5, cex = c(1, 2, 3, 0, 0, 0, 0))
          points(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), c(1, 1, 1, -1, -1, -1, -1), pch = 2, col = 2, cex = c(1, 2, 3, 0, 0, 0, 0))
          points(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), c(2, 2, 2, -1, -1, -1, -1), pch = 3, col = 3, cex = c(1, 2, 3, 0, 0, 0, 0))
          points(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), c(3, 3, 3, -1, -1, -1, -1), pch = 4, col = 4, cex = c(1, 2, 3, 0, 0, 0, 0))
          points(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), c(4, 4, 4, -1, -1, -1, -1), pch = 6, col = 6, cex = c(1, 2, 3, 0, 0, 0, 0))
          points(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), c(5, 5, 5, -1, -1, -1, -1), pch = 1, col = 1, cex = c(1, 2, 3, 0, 0, 0, 0))
          points(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), c(8, 8, 8, -1, -1, -1, -1), pch = 19, col = 2, cex = c(2, 2, 2, 0, 0, 0, 0))
          legend(x = 0.35, y = 1.5, legend = "- SA 1 local optimization accepted", bty = "n")
          legend(x = 0.35, y = 2.5, legend = "- MTMCMC local optimization accepted", bty = "n")
          legend(x = 0.35, y = 3.5, legend = "- GREEDY local optimization accepted", bty = "n")
          legend(x = 0.35, y = 4.5, legend = "- SA 2 local optimization accepted", bty = "n")
          legend(x = 0.35, y = 5.5, legend = "- NO local optimization accepted", bty = "n")
          legend(x = 0.35, y = 6.5, legend = "- Totally accepted", bty = "n")
          legend(x = 0.35, y = 7.5, legend = "- Totally explored", bty = "n")
          legend(x = 0.35, y = 8.5, legend = "- Bayes formula based posterior", bty = "n")
          dev.off()


          if (printable.opt) print(paste("drawing ", workdir, template, "_mlik.jpg", sep = ""))
          jpeg(file = paste(workdir, template, "_mlik.jpg", sep = ""))
          plot(ylim = mlik.lim, xlab = "model_id", ylab = "MLIK", statistics1[, 1], pch = 19, col = 7, cex = 1 * ((statistics1[, 9] + statistics1[, 5] + statistics1[, 6] + statistics1[, 7] + statistics1[, 8]) > 0))
          points(statistics1[, 1], pch = 8, col = ifelse(statistics1[, 4] > 0, 5, 0), cex = ifelse(statistics1[, 4] > 0, statistics1[, 4] / norm + 1, 0))
          points(statistics1[, 1], pch = 2, col = ifelse(statistics1[, 10] > 0, 2, 0), cex = ifelse(statistics1[, 10] > 0, statistics1[, 10] / norm + 1, 0))
          points(statistics1[, 1], pch = 3, col = ifelse(statistics1[, 11] > 0, 3, 0), cex = ifelse(statistics1[, 11] > 0, statistics1[, 11] / norm + 1, 0))
          points(statistics1[, 1], pch = 4, col = ifelse(statistics1[, 12] > 0, 4, 0), cex = ifelse(statistics1[, 12] > 0, statistics1[, 12] / norm + 1, 0))
          points(statistics1[, 1], pch = 6, col = ifelse(statistics1[, 13] > 0, 6, 0), cex = ifelse(statistics1[, 13] > 0, statistics1[, 13] / norm + 1, 0))
          points(statistics1[, 1], pch = 1, col = ifelse(statistics1[, 14] > 0, 1, 0), cex = ifelse(statistics1[, 14] > 0, statistics1[, 14] / norm + 1, 0))
          dev.off()
        })), abort = function() {
          onerr <- TRUE
        })
      })
    }

    if (crit$waic) {
      if (printable.opt) print(paste("drawing ", workdir, template, "_waic.jpg", sep = ""))
      utils::capture.output({
        withRestarts(tryCatch(utils::capture.output({
          jpeg(file = paste(workdir, template, "_waic.jpg", sep = ""))
          plot(ylim = waic.lim, xlab = "model_id", ylab = "WAIC", statistics1[, 2], pch = 19, col = 7, cex = 1 * ((statistics1[, 9] + statistics1[, 5] + statistics1[, 6] + statistics1[, 7] + statistics1[, 8]) > 0))
          points(statistics1[, 2], pch = 8, col = ifelse(statistics1[, 4] > 0, 5, 0), cex = ifelse(statistics1[, 4] > 0, statistics1[, 4] / norm + 1, 0))
          points(statistics1[, 2], pch = 2, col = ifelse(statistics1[, 10] > 0, 2, 0), cex = ifelse(statistics1[, 10] > 0, statistics1[, 10] / norm + 1, 0))
          points(statistics1[, 2], pch = 3, col = ifelse(statistics1[, 11] > 0, 3, 0), cex = ifelse(statistics1[, 11] > 0, statistics1[, 11] / norm + 1, 0))
          points(statistics1[, 2], pch = 4, col = ifelse(statistics1[, 12] > 0, 4, 0), cex = ifelse(statistics1[, 12] > 0, statistics1[, 12] / norm + 1, 0))
          points(statistics1[, 2], pch = 6, col = ifelse(statistics1[, 13] > 0, 6, 0), cex = ifelse(statistics1[, 13] > 0, statistics1[, 13] / norm + 1, 0))
          points(statistics1[, 2], pch = 1, col = ifelse(statistics1[, 14] > 0, 1, 0), cex = ifelse(statistics1[, 14] > 0, statistics1[, 14] / norm + 1, 0))
          dev.off()
        })), abort = function() {
          onerr <- TRUE
        })
      })
    }

    if (crit$dic) {
      if (printable.opt) print(paste("drawing ", workdir, template, "_dic.jpg", sep = ""))
      utils::capture.output({
        withRestarts(tryCatch(utils::capture.output({
          jpeg(file = paste(workdir, template, "_dic.jpg", sep = ""))
          plot(ylim = dic.lim, xlab = "model_id", ylab = "DIC", statistics1[, 3], pch = 19, col = 7, cex = 1 * ((statistics1[, 9] + statistics1[, 5] + statistics1[, 6] + statistics1[, 7] + statistics1[, 8]) > 0))
          points(statistics1[, 3], pch = 8, col = ifelse(statistics1[, 4] > 0, 5, 0), cex = ifelse(statistics1[, 4] > 0, statistics1[, 4] / norm + 1, 0))
          points(statistics1[, 3], pch = 2, col = ifelse(statistics1[, 10] > 0, 2, 0), cex = ifelse(statistics1[, 10] > 0, statistics1[, 10] / norm + 1, 0))
          points(statistics1[, 3], pch = 3, col = ifelse(statistics1[, 11] > 0, 3, 0), cex = ifelse(statistics1[, 11] > 0, statistics1[, 11] / norm + 1, 0))
          points(statistics1[, 3], pch = 4, col = ifelse(statistics1[, 12] > 0, 4, 0), cex = ifelse(statistics1[, 12] > 0, statistics1[, 12] / norm + 1, 0))
          points(statistics1[, 3], pch = 6, col = ifelse(statistics1[, 13] > 0, 6, 0), cex = ifelse(statistics1[, 13] > 0, statistics1[, 13] / norm + 1, 0))
          points(statistics1[, 3], pch = 1, col = ifelse(statistics1[, 14] > 0, 1, 0), cex = ifelse(statistics1[, 14] > 0, statistics1[, 14] / norm + 1, 0))
          dev.off()
        })), abort = function() {
          onerr <- TRUE
        })
      })
    }

    if (crit$waic && crit$mlik) {
      if (printable.opt) print(paste("drawing ", workdir, template, "_mlik-waic.jpg", sep = ""))
      utils::capture.output({
        withRestarts(tryCatch(utils::capture.output({
          jpeg(file = paste(workdir, template, "_mlik-waic.jpg", sep = ""))
          plot(ylim = mlik.lim, xlim = waic.lim, xlab = "WAIC", ylab = "MLIK", statistics1[, 2], statistics1[, 1], pch = 19, col = 7, cex = 1 * ((statistics1[, 9] + statistics1[, 5] + statistics1[, 6] + statistics1[, 7] + statistics1[, 8]) > 0))
          points(statistics1[, 2], statistics1[, 1], pch = 8, col = ifelse(statistics1[, 4] > 0, 5, 0), cex = ifelse(statistics1[, 4] > 0, statistics1[, 4] / norm + 1, 0))
          points(statistics1[, 2], statistics1[, 1], pch = 2, col = ifelse(statistics1[, 10] > 0, 2, 0), cex = ifelse(statistics1[, 10] > 0, statistics1[, 10] / norm + 1, 0))
          points(statistics1[, 2], statistics1[, 1], pch = 3, col = ifelse(statistics1[, 11] > 0, 3, 0), cex = ifelse(statistics1[, 11] > 0, statistics1[, 11] / norm + 1, 0))
          points(statistics1[, 2], statistics1[, 1], pch = 4, col = ifelse(statistics1[, 12] > 0, 4, 0), cex = ifelse(statistics1[, 12] > 0, statistics1[, 12] / norm + 1, 0))
          points(statistics1[, 2], statistics1[, 1], pch = 6, col = ifelse(statistics1[, 13] > 0, 6, 0), cex = ifelse(statistics1[, 13] > 0, statistics1[, 13] / norm + 1, 0))
          points(statistics1[, 2], statistics1[, 1], pch = 1, col = ifelse(statistics1[, 14] > 0, 1, 0), cex = ifelse(statistics1[, 14] > 0, statistics1[, 14] / norm + 1, 0))
          dev.off()
        })), abort = function() {
          onerr <- TRUE
        })
      })
    }

    if (crit$dic && crit$mlik) {
      if (printable.opt) print(paste("drawing ", workdir, template, "_mlik-dic.jpg", sep = ""))
      utils::capture.output({
        withRestarts(tryCatch(utils::capture.output({
          jpeg(file = paste(workdir, template, "_mlik-dic.jpg", sep = ""))
          plot(ylim = mlik.lim, xlim = dic.lim, xlab = "DIC", ylab = "MLIK", statistics1[, 3], statistics1[, 1], pch = 19, col = 7, cex = 1 * ((statistics1[, 9] + statistics1[, 5] + statistics1[, 6] + statistics1[, 7] + statistics1[, 8]) > 0))
          points(statistics1[, 3], statistics1[, 1], pch = 8, col = ifelse(statistics1[, 4] > 0, 5, 0), cex = ifelse(statistics1[, 4] > 0, statistics1[, 4] / norm + 1, 0))
          points(statistics1[, 3], statistics1[, 1], pch = 2, col = ifelse(statistics1[, 10] > 0, 2, 0), cex = ifelse(statistics1[, 10] > 0, statistics1[, 10] / norm + 1, 0))
          points(statistics1[, 3], statistics1[, 1], pch = 3, col = ifelse(statistics1[, 11] > 0, 3, 0), cex = ifelse(statistics1[, 11] > 0, statistics1[, 11] / norm + 1, 0))
          points(statistics1[, 3], statistics1[, 1], pch = 4, col = ifelse(statistics1[, 12] > 0, 4, 0), cex = ifelse(statistics1[, 12] > 0, statistics1[, 12] / norm + 1, 0))
          points(statistics1[, 3], statistics1[, 1], pch = 6, col = ifelse(statistics1[, 13] > 0, 6, 0), cex = ifelse(statistics1[, 13] > 0, statistics1[, 13] / norm + 1, 0))
          points(statistics1[, 3], statistics1[, 1], pch = 1, col = ifelse(statistics1[, 14] > 0, 1, 0), cex = ifelse(statistics1[, 14] > 0, statistics1[, 14] / norm + 1, 0))
          dev.off()
        })), abort = function() {
          onerr <- TRUE
        })
      })
    }

    if (crit$dic && crit$waic) {
      if (printable.opt) print(paste("drawing ", workdir, template, "_waic-dic.jpg", sep = ""))
      utils::capture.output({
        withRestarts(tryCatch(utils::capture.output({
          jpeg(file = paste(workdir, template, "_waic-dic.jpg", sep = ""))
          plot(ylim = waic.lim, xlim = dic.lim, xlab = "WAIC", ylab = "DIC", statistics1[, 3], statistics1[, 2], pch = 19, col = 7, cex = 1 * ((statistics1[, 9] + statistics1[, 5] + statistics1[, 6] + statistics1[, 7] + statistics1[, 8]) > 0))
          points(statistics1[, 3], statistics1[, 2], pch = 8, col = ifelse(statistics1[, 4] > 0, 5, 0), cex = ifelse(statistics1[, 4] > 0, statistics1[, 4] / norm + 1, 0))
          points(statistics1[, 3], statistics1[, 2], pch = 2, col = ifelse(statistics1[, 10] > 0, 2, 0), cex = ifelse(statistics1[, 10] > 0, statistics1[, 10] / norm + 1, 0))
          points(statistics1[, 3], statistics1[, 2], pch = 3, col = ifelse(statistics1[, 11] > 0, 3, 0), cex = ifelse(statistics1[, 11] > 0, statistics1[, 11] / norm + 1, 0))
          points(statistics1[, 3], statistics1[, 2], pch = 4, col = ifelse(statistics1[, 12] > 0, 4, 0), cex = ifelse(statistics1[, 12] > 0, statistics1[, 12] / norm + 1, 0))
          points(statistics1[, 3], statistics1[, 2], pch = 6, col = ifelse(statistics1[, 13] > 0, 6, 0), cex = ifelse(statistics1[, 13] > 0, statistics1[, 13] / norm + 1, 0))
          points(statistics1[, 3], statistics1[, 2], pch = 1, col = ifelse(statistics1[, 14] > 0, 1, 0), cex = ifelse(statistics1[, 14] > 0, statistics1[, 14] / norm + 1, 0))
          dev.off()
        })), abort = function() {
          onerr <- TRUE
        })
      })
    }


    norm1 <- (sum(statistics1[, 4], na.rm = TRUE))

    xyz <- which(!is.na(statistics1[, 1]))
    xyz <- intersect(xyz, which(statistics1[, 1] != -10000))
    moddee <- which(statistics1[, 1] == max(statistics1[, 1], na.rm = TRUE))[1]
    zyx <- array(data = NA, dim = 2^(Nvars) + 1)
    nconsum <- sum(exp(-statistics1[moddee, 1] + statistics1[xyz, 1]), na.rm = TRUE)
    if (nconsum > 0) {
      zyx[xyz] <- exp(statistics1[xyz, 1] - statistics1[moddee, 1]) / nconsum
      y.post.lim <- c(0, max(zyx[xyz]))
    } else {
      zyx[xyz] <- statistics1[xyz, 3] / norm1
      y.post.lim <- c(0, NaN)
    }


    if (is.nan(y.post.lim[2])) {
      y.post.lim[2] <- max(statistics1[, 3] / norm1, na.rm = TRUE)
    }
    if (is.nan(y.post.lim[2])) {
      y.post.lim[2] <- 1
    }

    if (crit$mlik) {
      if (printable.opt) print(paste("drawing ", workdir, template, "_Pr(MID).jpg", sep = ""))
      utils::capture.output({
        withRestarts(tryCatch(utils::capture.output({
          jpeg(file = paste(workdir, template, "_Pr(MID).jpg", sep = ""))
          plot(xlab = "model_id", ylab = "Pr(M(model_id)ID)", ylim = y.post.lim, statistics1[, 4] / norm1, pch = 19, col = 7, cex = 3 * ((statistics1[, 9] + statistics1[, 5] + statistics1[, 6] + statistics1[, 7] + statistics1[, 8]) > 0))
          points(statistics1[, 4] / norm1, pch = 8, col = ifelse(statistics1[, 4] > 0, 5, 0), cex = ifelse(statistics1[, 4] > 0, statistics1[, 4] / norm + 1, 0))
          points(statistics1[, 4] / norm1, pch = 2, col = ifelse(statistics1[, 10] > 0, 2, 0), cex = ifelse(statistics1[, 10] > 0, statistics1[, 10] / norm + 1, 0))
          points(statistics1[, 4] / norm1, pch = 3, col = ifelse(statistics1[, 11] > 0, 3, 0), cex = ifelse(statistics1[, 11] > 0, statistics1[, 11] / norm + 1, 0))
          points(statistics1[, 4] / norm1, pch = 4, col = ifelse(statistics1[, 12] > 0, 4, 0), cex = ifelse(statistics1[, 12] > 0, statistics1[, 12] / norm + 1, 0))
          points(statistics1[, 4] / norm1, pch = 6, col = ifelse(statistics1[, 13] > 0, 6, 0), cex = ifelse(statistics1[, 13] > 0, statistics1[, 13] / norm + 1, 0))
          points(statistics1[, 4] / norm1, pch = 1, col = ifelse(statistics1[, 14] > 0, 1, 0), cex = ifelse(statistics1[, 14] > 0, statistics1[, 14] / norm + 1, 0))
          points(zyx, pch = 20, col = 2, cex = 2)
          dev.off()
        })), abort = function() {
          onerr <- TRUE
        })
      })
    }

    if (crit$waic && crit$mlik) {
      if (printable.opt) print(paste("drawing ", workdir, template, "_waic-Pr(MID).jpg", sep = ""))
      utils::capture.output({
        withRestarts(tryCatch(utils::capture.output({
          jpeg(file = paste(workdir, template, "_waic-Pr(MID).jpg", sep = ""))
          plot(xlim = waic.lim, ylim = y.post.lim, xlab = "Wstats::AIC(M)", ylab = "Pr(MID)", statistics1[, 2], statistics1[, 4] / norm1, pch = 19, col = 7, cex = 3 * ((statistics1[, 9] + statistics1[, 5] + statistics1[, 6] + statistics1[, 7] + statistics1[, 8]) > 0))
          points(statistics1[, 2], statistics1[, 4] / norm1, pch = 8, col = ifelse(statistics1[, 4] > 0, 5, 0), cex = ifelse(statistics1[, 4] > 0, statistics1[, 4] / norm + 1, 0))
          points(statistics1[, 2], statistics1[, 4] / norm1, pch = 2, col = ifelse(statistics1[, 10] > 0, 2, 0), cex = ifelse(statistics1[, 10] > 0, statistics1[, 10] / norm + 1, 0))
          points(statistics1[, 2], statistics1[, 4] / norm1, pch = 3, col = ifelse(statistics1[, 11] > 0, 3, 0), cex = ifelse(statistics1[, 11] > 0, statistics1[, 11] / norm + 1, 0))
          points(statistics1[, 2], statistics1[, 4] / norm1, pch = 4, col = ifelse(statistics1[, 12] > 0, 4, 0), cex = ifelse(statistics1[, 12] > 0, statistics1[, 12] / norm + 1, 0))
          points(statistics1[, 2], statistics1[, 4] / norm1, pch = 6, col = ifelse(statistics1[, 13] > 0, 6, 0), cex = ifelse(statistics1[, 13] > 0, statistics1[, 13] / norm + 1, 0))
          points(statistics1[, 2], statistics1[, 4] / norm1, pch = 1, col = ifelse(statistics1[, 14] > 0, 1, 0), cex = ifelse(statistics1[, 14] > 0, statistics1[, 14] / norm + 1, 0))
          points(statistics1[, 2], zyx, pch = 20, col = 2, cex = 2)
          dev.off()
        })), abort = function() {
          onerr <- TRUE
        })
      })
    }

    if (crit$dic && crit$mlik) {
      if (printable.opt) print(paste("drawing ", workdir, template, "_dic-Pr(MID).jpg", sep = ""))
      utils::capture.output({
        withRestarts(tryCatch(utils::capture.output({
          jpeg(file = paste(workdir, template, "_dic-Pr(MID).jpg", sep = ""))
          plot(xlim = dic.lim, ylim = y.post.lim, xlab = "dic(M)", ylab = "Pr(MID)", statistics1[, 3], statistics1[, 4] / norm1, pch = 19, col = 7, cex = 3 * ((statistics1[, 9] + statistics1[, 5] + statistics1[, 6] + statistics1[, 7] + statistics1[, 8]) > 0))
          points(statistics1[, 3], statistics1[, 4] / norm1, pch = 8, col = ifelse(statistics1[, 4] > 0, 5, 0), cex = ifelse(statistics1[, 4] > 0, statistics1[, 4] / norm + 1, 0))
          points(statistics1[, 3], statistics1[, 4] / norm1, pch = 2, col = ifelse(statistics1[, 10] > 0, 2, 0), cex = ifelse(statistics1[, 10] > 0, statistics1[, 10] / norm + 1, 0))
          points(statistics1[, 3], statistics1[, 4] / norm1, pch = 3, col = ifelse(statistics1[, 11] > 0, 3, 0), cex = ifelse(statistics1[, 11] > 0, statistics1[, 11] / norm + 1, 0))
          points(statistics1[, 3], statistics1[, 4] / norm1, pch = 4, col = ifelse(statistics1[, 12] > 0, 4, 0), cex = ifelse(statistics1[, 12] > 0, statistics1[, 12] / norm + 1, 0))
          points(statistics1[, 3], statistics1[, 4] / norm1, pch = 6, col = ifelse(statistics1[, 13] > 0, 6, 0), cex = ifelse(statistics1[, 13] > 0, statistics1[, 13] / norm + 1, 0))
          points(statistics1[, 3], statistics1[, 4] / norm1, pch = 1, col = ifelse(statistics1[, 14] > 0, 1, 0), cex = ifelse(statistics1[, 14] > 0, statistics1[, 14] / norm + 1, 0))
          points(statistics1[, 3], zyx, pch = 20, col = 2, cex = 2)
          dev.off()
        })), abort = function() {
          onerr <- TRUE
        })
      })
    }


    if (crit$mlik) {
      jpeg(file = paste(workdir, template, "_mlik-Pr(MID).jpg", sep = ""))
      utils::capture.output({
        withRestarts(tryCatch(utils::capture.output({
          plot(xlim = mlik.lim, ylim = y.post.lim, xlab = "MLIK", ylab = "Pr(MID)", statistics1[, 1], statistics1[, 4] / norm1, pch = 19, col = 7, cex = 3 * ((statistics1[, 9] + statistics1[, 5] + statistics1[, 6] + statistics1[, 7] + statistics1[, 8]) > 0))
          points(statistics1[, 1], statistics1[, 4] / norm1, pch = 8, col = ifelse(statistics1[, 4] > 0, 5, 0), cex = ifelse(statistics1[, 4] > 0, statistics1[, 4] / statistics1[, 3] * 3 + 1, 0))
          points(statistics1[, 1], statistics1[, 4] / norm1, pch = 2, col = ifelse(statistics1[, 10] > 0, 2, 0), cex = ifelse(statistics1[, 10] > 0, statistics1[, 10] / statistics1[, 4] * 3 + 1, 0))
          points(statistics1[, 1], statistics1[, 4] / norm1, pch = 3, col = ifelse(statistics1[, 11] > 0, 3, 0), cex = ifelse(statistics1[, 11] > 0, statistics1[, 11] / statistics1[, 4] * 3 + 1, 0))
          points(statistics1[, 1], statistics1[, 4] / norm1, pch = 4, col = ifelse(statistics1[, 12] > 0, 4, 0), cex = ifelse(statistics1[, 12] > 0, statistics1[, 12] / statistics1[, 4] * 3 + 1, 0))
          points(statistics1[, 1], statistics1[, 4] / norm1, pch = 6, col = ifelse(statistics1[, 13] > 0, 6, 0), cex = ifelse(statistics1[, 13] > 0, statistics1[, 13] / statistics1[, 4] * 3 + 1, 0))
          points(statistics1[, 1], statistics1[, 4] / norm1, pch = 1, col = ifelse(statistics1[, 14] > 0, 1, 0), cex = ifelse(statistics1[, 14] > 0, statistics1[, 14] / statistics1[, 4] * 3 + 1, 0))
          points(statistics1[, 1], zyx, pch = 20, col = 2, cex = 2)
          dev.off()
        })), abort = function() {
          onerr <- TRUE
        })
      })
    }

    if (draw_dist) {
      if (printable.opt) print("Calculating distance matrix, may take a significant amount of time, may also produce errors if your machine does not have enough memory")
      utils::capture.output({
        withRestarts(tryCatch(utils::capture.output({
          lldd <- 2^(Nvars) + 1

          moddee <- which(zyx == max(zyx, na.rm = TRUE))
          iidr <- which(statistics1[moddee, 1] == max(statistics1[moddee, 1], na.rm = TRUE))
          iidr <- which(statistics1[moddee[iidr], 3] == max(statistics1[moddee[iidr], 3], na.rm = TRUE))
          moddee <- moddee[iidr]
          if (length(moddee) > 1) {
            moddee <- moddee[1]
          }

          vec <- dectobit.alt(moddee - 1)
          varcur <- c(array(0, dim = (Nvars - length(vec))), vec)
          df <- data.frame(varcur)


          for (i in 1:(lldd - 1))
          {
            if (i == moddee) {
              next
            } else {
              vec <- dectobit.alt(i - 1)
              varcur <- c(array(0, dim = (Nvars - length(vec))), vec)
              df <- cbind(df, varcur)
              # colnames(x = df)[i] <- paste("solution ",i)
            }
          }
          df <- t(df)


          x <- dist(x = df, method = "binary")

          dists <- c(0, x[1:lldd - 1])

          # length(dists)
          # which(dists==0)
        })), abort = function() {
          onerr <- TRUE
        })
      })

      if (crit$mlik) {
        if (printable.opt) print(paste("drawing ", workdir, template, "_distance-mlik.jpg", sep = ""))
        utils::capture.output({
          withRestarts(tryCatch(utils::capture.output({
            jpeg(file = paste(workdir, template, "_distance-MLIK.jpg", sep = ""))
            plot(xlab = "x = |M* - M| = distance from the main mode", ylab = "MLIK(M)", ylim = mlik.lim, y = c(statistics1[moddee, 1], statistics1[-moddee, 1]), x = dists, pch = 19, col = 7, cex = 1 * ((statistics1[, 9] + statistics1[, 5] + statistics1[, 6] + statistics1[, 7] + statistics1[, 8]) > 0))
            points(y = c(statistics1[moddee, 1], statistics1[-moddee, 1]), x = dists, pch = 8, col = ifelse(c(statistics1[moddee, 4], statistics1[-moddee, 4]) > 0, 5, 0), cex = ifelse(c(statistics1[moddee, 4], statistics1[-moddee, 4]) > 0, c(statistics1[moddee, 4], statistics1[-moddee, 4]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 1], statistics1[-moddee, 1]), x = dists, pch = 2, col = ifelse(c(statistics1[moddee, 10], statistics1[-moddee, 10]) > 0, 2, 0), cex = ifelse(c(statistics1[moddee, 10], statistics1[-moddee, 10]) > 0, c(statistics1[moddee, 10], statistics1[-moddee, 10]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 1], statistics1[-moddee, 1]), x = dists, pch = 3, col = ifelse(c(statistics1[moddee, 11], statistics1[-moddee, 11]) > 0, 3, 0), cex = ifelse(c(statistics1[moddee, 11], statistics1[-moddee, 11]) > 0, c(statistics1[moddee, 11], statistics1[-moddee, 11]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 1], statistics1[-moddee, 1]), x = dists, pch = 4, col = ifelse(c(statistics1[moddee, 12], statistics1[-moddee, 12]) > 0, 4, 0), cex = ifelse(c(statistics1[moddee, 12], statistics1[-moddee, 12]) > 0, c(statistics1[moddee, 12], statistics1[-moddee, 12]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 1], statistics1[-moddee, 1]), x = dists, pch = 6, col = ifelse(c(statistics1[moddee, 13], statistics1[-moddee, 13]) > 0, 6, 0), cex = ifelse(c(statistics1[moddee, 13], statistics1[-moddee, 13]) > 0, c(statistics1[moddee, 13], statistics1[-moddee, 13]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 1], statistics1[-moddee, 1]), x = dists, pch = 1, col = ifelse(c(statistics1[moddee, 14], statistics1[-moddee, 14]) > 0, 1, 0), cex = ifelse(c(statistics1[moddee, 14], statistics1[-moddee, 14]) > 0, c(statistics1[moddee, 14], statistics1[-moddee, 14]) / norm + 1, 0))
            dev.off()
          })), abort = function() {
            onerr <- TRUE
          })
        })
      }
      if (crit$waic) {
        if (printable.opt) print(paste("drawing ", workdir, template, "_distance-waic.jpg", sep = ""))
        utils::capture.output({
          withRestarts(tryCatch(utils::capture.output({
            jpeg(file = paste(workdir, template, "_distance-waic.jpg", sep = ""))
            plot(xlab = "x = |M* - M| = distance from the main mode", ylab = "Wstats::AIC(M)", ylim = waic.lim, y = c(statistics1[moddee, 2], statistics1[-moddee, 2]), x = dists, pch = 19, col = 7, cex = 1 * ((statistics1[, 9] + statistics1[, 5] + statistics1[, 6] + statistics1[, 7] + statistics1[, 8]) > 0))
            points(y = c(statistics1[moddee, 2], statistics1[-moddee, 2]), x = dists, pch = 8, col = ifelse(c(statistics1[moddee, 4], statistics1[-moddee, 4]) > 0, 5, 0), cex = ifelse(c(statistics1[moddee, 4], statistics1[-moddee, 4]) > 0, c(statistics1[moddee, 4], statistics1[-moddee, 4]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 2], statistics1[-moddee, 2]), x = dists, pch = 2, col = ifelse(c(statistics1[moddee, 10], statistics1[-moddee, 10]) > 0, 2, 0), cex = ifelse(c(statistics1[moddee, 10], statistics1[-moddee, 10]) > 0, c(statistics1[moddee, 10], statistics1[-moddee, 10]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 2], statistics1[-moddee, 2]), x = dists, pch = 3, col = ifelse(c(statistics1[moddee, 11], statistics1[-moddee, 11]) > 0, 3, 0), cex = ifelse(c(statistics1[moddee, 11], statistics1[-moddee, 11]) > 0, c(statistics1[moddee, 11], statistics1[-moddee, 11]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 2], statistics1[-moddee, 2]), x = dists, pch = 4, col = ifelse(c(statistics1[moddee, 12], statistics1[-moddee, 12]) > 0, 4, 0), cex = ifelse(c(statistics1[moddee, 12], statistics1[-moddee, 12]) > 0, c(statistics1[moddee, 12], statistics1[-moddee, 12]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 2], statistics1[-moddee, 2]), x = dists, pch = 6, col = ifelse(c(statistics1[moddee, 13], statistics1[-moddee, 13]) > 0, 6, 0), cex = ifelse(c(statistics1[moddee, 13], statistics1[-moddee, 13]) > 0, c(statistics1[moddee, 13], statistics1[-moddee, 13]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 2], statistics1[-moddee, 2]), x = dists, pch = 1, col = ifelse(c(statistics1[moddee, 14], statistics1[-moddee, 14]) > 0, 1, 0), cex = ifelse(c(statistics1[moddee, 14], statistics1[-moddee, 14]) > 0, c(statistics1[moddee, 14], statistics1[-moddee, 14]) / norm + 1, 0))
            dev.off()
          })), abort = function() {
            onerr <- TRUE
          })
        })
      }
      if (crit$dic) {
        if (printable.opt) print(paste("drawing ", workdir, template, "_distance-dic.jpg", sep = ""))
        utils::capture.output({
          withRestarts(tryCatch(utils::capture.output({
            jpeg(file = paste(workdir, template, "_distance-dic.jpg", sep = ""))
            plot(xlab = "x = |M* - M| = distance from the main mode", ylab = "DIC(M)", ylim = dic.lim, y = c(statistics1[moddee, 3], statistics1[-moddee, 3]), x = dists, pch = 19, col = 7, cex = 1 * ((statistics1[, 9] + statistics1[, 5] + statistics1[, 6] + statistics1[, 7] + statistics1[, 8]) > 0))
            points(y = c(statistics1[moddee, 3], statistics1[-moddee, 3]), x = dists, pch = 8, col = ifelse(c(statistics1[moddee, 4], statistics1[-moddee, 4]) > 0, 5, 0), cex = ifelse(c(statistics1[moddee, 4], statistics1[-moddee, 4]) > 0, c(statistics1[moddee, 4], statistics1[-moddee, 4]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 3], statistics1[-moddee, 3]), x = dists, pch = 2, col = ifelse(c(statistics1[moddee, 10], statistics1[-moddee, 10]) > 0, 2, 0), cex = ifelse(c(statistics1[moddee, 10], statistics1[-moddee, 10]) > 0, c(statistics1[moddee, 10], statistics1[-moddee, 10]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 3], statistics1[-moddee, 3]), x = dists, pch = 3, col = ifelse(c(statistics1[moddee, 11], statistics1[-moddee, 11]) > 0, 3, 0), cex = ifelse(c(statistics1[moddee, 11], statistics1[-moddee, 11]) > 0, c(statistics1[moddee, 11], statistics1[-moddee, 11]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 3], statistics1[-moddee, 3]), x = dists, pch = 4, col = ifelse(c(statistics1[moddee, 12], statistics1[-moddee, 12]) > 0, 4, 0), cex = ifelse(c(statistics1[moddee, 12], statistics1[-moddee, 12]) > 0, c(statistics1[moddee, 12], statistics1[-moddee, 12]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 3], statistics1[-moddee, 3]), x = dists, pch = 6, col = ifelse(c(statistics1[moddee, 13], statistics1[-moddee, 13]) > 0, 6, 0), cex = ifelse(c(statistics1[moddee, 13], statistics1[-moddee, 13]) > 0, c(statistics1[moddee, 13], statistics1[-moddee, 13]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 3], statistics1[-moddee, 3]), x = dists, pch = 1, col = ifelse(c(statistics1[moddee, 14], statistics1[-moddee, 14]) > 0, 1, 0), cex = ifelse(c(statistics1[moddee, 14], statistics1[-moddee, 14]) > 0, c(statistics1[moddee, 14], statistics1[-moddee, 14]) / norm + 1, 0))
            dev.off()
          })), abort = function() {
            onerr <- TRUE
          })
        })
      }

      if (crit$mlik) {
        if (printable.opt) print(paste("drawing ", workdir, template, "_distance-Pr(MID).jpg", sep = ""))
        utils::capture.output({
          withRestarts(tryCatch(utils::capture.output({
            jpeg(file = paste(workdir, template, "_distance-Pr(MID).jpg", sep = ""))
            plot(xlab = "x = |M* - M| = distance from the main mode", ylim = y.post.lim, ylab = "Pr(MID)", y = c(statistics1[moddee, 4], statistics1[-moddee, 4]) / norm1, x = dists, pch = 19, col = 7, cex = 3 * ((statistics1[, 9] + statistics1[, 5] + statistics1[, 6] + statistics1[, 7] + statistics1[, 8]) > 0))
            points(y = c(statistics1[moddee, 4], statistics1[-moddee, 4]) / norm1, x = dists, pch = 8, col = ifelse(c(statistics1[moddee, 4], statistics1[-moddee, 4]) > 0, 5, 0), cex = ifelse(c(statistics1[moddee, 4], statistics1[-moddee, 4]) > 0, c(statistics1[moddee, 4], statistics1[-moddee, 4]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 4], statistics1[-moddee, 4]) / norm1, x = dists, pch = 2, col = ifelse(c(statistics1[moddee, 10], statistics1[-moddee, 10]) > 0, 2, 0), cex = ifelse(c(statistics1[moddee, 10], statistics1[-moddee, 10]) > 0, c(statistics1[moddee, 10], statistics1[-moddee, 10]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 4], statistics1[-moddee, 4]) / norm1, x = dists, pch = 3, col = ifelse(c(statistics1[moddee, 11], statistics1[-moddee, 11]) > 0, 3, 0), cex = ifelse(c(statistics1[moddee, 11], statistics1[-moddee, 11]) > 0, c(statistics1[moddee, 11], statistics1[-moddee, 11]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 4], statistics1[-moddee, 4]) / norm1, x = dists, pch = 4, col = ifelse(c(statistics1[moddee, 12], statistics1[-moddee, 12]) > 0, 4, 0), cex = ifelse(c(statistics1[moddee, 12], statistics1[-moddee, 12]) > 0, c(statistics1[moddee, 12], statistics1[-moddee, 12]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 4], statistics1[-moddee, 4]) / norm1, x = dists, pch = 6, col = ifelse(c(statistics1[moddee, 13], statistics1[-moddee, 13]) > 0, 6, 0), cex = ifelse(c(statistics1[moddee, 13], statistics1[-moddee, 13]) > 0, c(statistics1[moddee, 13], statistics1[-moddee, 13]) / norm + 1, 0))
            points(y = c(statistics1[moddee, 4], statistics1[-moddee, 4]) / norm1, x = dists, pch = 1, col = ifelse(c(statistics1[moddee, 14], statistics1[-moddee, 14]) > 0, 1, 0), cex = ifelse(c(statistics1[moddee, 14], statistics1[-moddee, 14]) > 0, c(statistics1[moddee, 14], statistics1[-moddee, 14]) / norm + 1, 0))
            points(y = c(zyx[moddee], zyx[-moddee]), x = dists, pch = 20, col = 2, cex = 2)
            dev.off()
          })), abort = function() {
            onerr <- TRUE
          })
        })
      }


      if (crit$mlik) {
        if (printable.opt) print(paste("drawing ", workdir, template, "_mds-Pr(MID).jpg", sep = ""))
        if (printable.opt) print("Calculating distance matrix, may take a significant amount of time, may also produce errors if your machine does not have enough memory")
        utils::capture.output({
          withRestarts(tryCatch(utils::capture.output({
            # further address subset of the set of the best solution of cardinality 1024

            if (lldd > mds_size) {
              lldd <- mds_size
              quant <- (sort(statistics1[, 1], decreasing = TRUE)[lldd + 1])
              indmds <- which(statistics1[, 1] > quant)
              length(indmds)
            } else {
              quant <- -Inf
              indmds <- 1:(lldd)
            }





            vec <- dectobit.alt(moddee - 1)
            varcur <- c(array(0, dim = (Nvars - length(vec))), vec)
            df <- data.frame(varcur)


            for (i in 1:(lldd - 1))
            {
              if (i == moddee) {
                next
              } else {
                vec <- dectobit.alt(indmds[i] - 1)
                varcur <- c(array(0, dim = (Nvars - length(vec))), vec)
                df <- cbind(df, varcur)
                # colnames(x = df)[i] <- paste("solution ",i)
              }
            }
            df <- t(df)


            x <- dist(x = df, method = "binary")

            dists <- c(0, x[1:lldd - 1])

            fit.mds <- cmdscale(d = x, eig = FALSE, k = 2) # k is the number of dim


            # fit.mds # view results
            x.mds <- fit.mds[, 1]
            y.mds <- fit.mds[, 2]
            jpeg(file = paste(workdir, template, "_mds_map_posteriors.jpg", sep = ""))
            plot(x.mds[], y.mds[], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Metric MDS", type = "p", col = 5, pch = 8, cex = ifelse(c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) > 0, c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) / norm1 * 50, 0))
            points(x.mds[], y.mds[], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Metric MDS", type = "p", col = 2, pch = 10, cex = c(zyx[moddee], zyx[setdiff(indmds, moddee)]) * 50, 0)
            points(x.mds[], y.mds[], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Metric MDS", type = "p", col = 2, pch = 10, cex = c(zyx[moddee], zyx[setdiff(indmds, moddee)]) * 50, 0)
            points(x.mds[], y.mds[], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Metric MDS", type = "p", col = 5, pch = 8, cex = ifelse(c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) > 0, c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) / norm1 * 50, 0))
            points(x.mds[], y.mds[], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Metric MDS", type = "p", col = 5, pch = 8, cex = ifelse(c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) > 0, c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) / norm1 * 50, 0))
            points(x.mds[], y.mds[], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Metric MDS", type = "p", col = 5, pch = 8, cex = ifelse(c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) > 0, c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) / norm1 * 50, 0))
            points(x.mds[], y.mds[], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Metric MDS", type = "p", col = 5, pch = 8, cex = ifelse(c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) > 0, c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) / norm1 * 50, 0))
            points(x.mds[], y.mds[], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Metric MDS", type = "p", col = 5, pch = 8, cex = ifelse(c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) > 0, c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) / norm1 * 50, 0))
            points(x.mds[], y.mds[], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Metric MDS", type = "p", col = 5, pch = 8, cex = ifelse(c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) > 0, c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) / norm1 * 50, 0))
            points(x.mds[], y.mds[], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Metric MDS", type = "p", col = 5, pch = 8, cex = ifelse(c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) > 0, c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) / norm1 * 50, 0))
            points(x.mds[], y.mds[], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Metric MDS", type = "p", col = 5, pch = 8, cex = ifelse(c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) > 0, c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) / norm1 * 50, 0))
            points(x.mds[], y.mds[], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Metric MDS", type = "p", col = 5, pch = 8, cex = ifelse(c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) > 0, c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) / norm1 * 50, 0))
            points(x.mds[], y.mds[], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Metric MDS", type = "p", col = 5, pch = 8, cex = ifelse(c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) > 0, c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) / norm1 * 50, 0))
            points(x.mds[], y.mds[], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Metric MDS", type = "p", col = 5, pch = 8, cex = ifelse(c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) > 0, c(statistics1[moddee, 4], statistics1[setdiff(indmds, moddee), 4]) / norm1 * 50, 0))
            points(x.mds[], y.mds[], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Metric MDS", type = "p", col = 7, pch = 19, cex = 0.4)
            points(x.mds[], y.mds[], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Metric MDS", type = "p", col = 1, pch = 19, cex = 0.01)
            dev.off()
          })), abort = function() {
            onerr <- TRUE
          })
        })
      }
    }
  }
)
