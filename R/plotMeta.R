plotMeta1 <- function (mat, centralTend = "mean", overlay = TRUE, winsorize = c(0, 
                                                                   100), profile.names = NULL, xcoords = NULL, meta.rescale = FALSE, 
          smoothfun = NULL, line.col = NULL, dispersion = NULL, dispersion.col = NULL, 
          ylim = NULL, ylab = "average score", xlab = "bases", ...)
{
  if (!class(mat) %in% c("ScoreMatrix", "ScoreMatrixList")) 
    stop("mat is not ScoreMatrix or ScoreMatrixList\n")
  if (!centralTend %in% c("median", "mean"))
    stop("centralTend is not mean or median\n")
  disp.args <- c("se", "sd", "IQR")
  if (!is.null(dispersion) && !dispersion %in% c(disp.args)) {
    stop("dispersion is not FALSE, 'se', 'sd' or 'IQR'\n")
  }
  if (!is.null(smoothfun) & !is.function(smoothfun)) {
    stop("'smoothfun' has to be a function or NULL\n")
  }
  if (is.null(line.col) & is.null(dispersion)) 
    line.col = ifelse(is.list(mat), list(rainbow(length(mat))), 
                      "black")[[1]]
  if (is.null(line.col) & is.null(dispersion.col)) {
    dispersion.col = ifelse(is.list(mat), list(rainbow(length(mat), 
                                                       alpha = 0.4)), rainbow(1, alpha = 0.4))[[1]]
    line.col = ifelse(is.list(mat), list(rainbow(length(mat))), 
                      rainbow(1))[[1]]
  }
  if (class(mat) == "ScoreMatrix") {
    mat <- list(mat)
  }
  if (length(unique(sapply(mat, ncol))) != 1) {
    stop("ScoreMatrix number of columns do not match\n", 
         "Try using binMatrix to make matrices with high number of columns", 
         "equal\n")
  }
  if (!is.null(dispersion) && dispersion %in% disp.args) {
    bound1 <- list()
    if (dispersion == "IQR") {
      q1 <- list()
      q3 <- list()
    }
    else {
      bound2 <- list()
    }
  }
  metas <- list()
  for (i in 1:length(mat)) {
    if (winsorize[2] < 100 | winsorize[1] > 0) {
      mat[[i]] = .winsorize(mat[[i]]@.Data, winsorize)
    }
    if (centralTend == "mean") {
      if (!is.null(dispersion) && dispersion == "IQR") {
        warning("dispersion is set to show 1st and 3rd quartile and \n                confidence interval around the median, \n                but centralTend is 'mean'. Setting centralTend to 'median'..\n")
        metas[[i]] = colMedians(mat[[i]], na.rm = TRUE)
      }
      else {
        cat("here\n")
        metas[[i]] = colMeans(mat[[i]], na.rm = TRUE)
      }
    }
    else if (centralTend == "median") {
      if (!is.null(dispersion) && dispersion == "se") {
        warning("dispersion is set to standard error of the mean and 95% confidence interval for the mean, but\n                centralTend is 'median'. Setting centralTend to 'mean'..\n")
        metas[[i]] = colMeans(mat[[i]], na.rm = TRUE)
      }
      else {
        metas[[i]] = colMedians(mat[[i]], na.rm = TRUE)
      }
    }
    if (!is.null(dispersion) && dispersion %in% disp.args) {
      if (dispersion == "se") {
        bound1[[i]] <- std.error(mat[[i]], na.rm = TRUE)
        bound2[[i]] <- bound1[[i]] * 1.96
      }
      else if (dispersion == "sd") {
        bound1[[i]] <- colSds(mat[[i]], na.rm = TRUE)
        bound2[[i]] <- bound1[[i]] * 2
      }
      else if (dispersion == "IQR") {
        q <- colQuantiles(mat[[i]], probs = c(0.25, 0.75), 
                          na.rm = TRUE)
        q1[[i]] <- q[, 1]
        q3[[i]] <- q[, 2]
        n <- nrow(mat[[i]])
        bound1[[i]] <- (1.57 * (q3[[i]] - q1[[i]]))/sqrt(n)
      }
    }
    if (meta.rescale) {
      val2unit <- function(x) {
        (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - 
                                      min(x, na.rm = TRUE))
      }
      metas[[i]] = val2unit(metas[[i]])
      if (!is.null(dispersion) && dispersion %in% disp.args) {
        bound1[[i]] = val2unit(bound1[[i]])
        if (dispersion == "IQR") {
          q1[[i]] = val2unit(q1[[i]])
          q3[[i]] = val2unit(q3[[i]])
        }
        else {
          bound2[[i]] = val2unit(bound2[[i]])
        }
      }
    }
    if (!is.null(smoothfun)) {
      metas[[i]] <- smoothfun(metas[[i]])$y
      if (!is.null(dispersion) && dispersion %in% disp.args) {
        bound1[[i]] <- smoothfun(bound1[[i]])$y
        if (dispersion == "IQR") {
          q1[[i]] <- smoothfun(q1[[i]])$y
          q3[[i]] <- smoothfun(q3[[i]])$y
        }
        else {
          bound2[[i]] <- smoothfun(bound2[[i]])$y
        }
      }
    }
  }
  if (!is.null(xcoords)) {
    if (length(xcoords) == 2 & xcoords[1] < xcoords[2]) {
      xcoords = seq(xcoords[1], xcoords[2], length.out = length(metas[[1]]))
    }
    if (length(xcoords) != length(metas[[1]])) 
      stop("xcoords has wrong length: ", length(xcoords), 
           " \n", "it should be equal to the number of columns of ScoreMatrices\n", 
           "which is: ", length(metas[[1]]), "\n")
  }
  else {
    xcoords = 1:length(metas[[1]])
  }
  if (!is.null(ylim)) {
    myrange = ylim
  }
  else {
    myrange = range(unlist(metas), na.rm = TRUE)
    if (!is.null(dispersion) && dispersion %in% disp.args) {
      if (dispersion != "IQR") {
        myrange[2] <- max(unlist(metas) + unlist(bound2), 
                          na.rm = TRUE)
        myrange[1] <- min(unlist(metas) - unlist(bound2), 
                          na.rm = TRUE)
      }
      else {
        max.q3 <- max(unlist(metas) + unlist(q3), na.rm = TRUE)
        min.q1 <- min(unlist(metas) - unlist(q1), na.rm = TRUE)
        max.notch.upper <- max(unlist(metas) + unlist(bound1), 
                               na.rm = TRUE)
        min.notch.lower <- min(unlist(metas) - unlist(bound1), 
                               na.rm = TRUE)
        if (min.notch.lower < min.q1) {
          myrange[1] <- min.notch.lower
        }
        else {
          myrange[1] <- min.q1
        }
        if (max.notch.upper > max.q3) {
          myrange[2] <- max.notch.upper
        }
        else {
          myrange[2] <- max.q3
        }
      }
    }
  }
  marOrg = par()$mar
  marNew = marOrg
  marNew[4] = 8
  par(mar = marNew)
  par(xpd = TRUE)
  if (overlay & length(metas) > 1) {
    
    
    if (!is.null(dispersion) && dispersion %in% disp.args) {
      plot(xcoords, metas[[1]], type = "l", col = dispersion.col[1], 
           ylim = myrange, ylab = ylab, xlab = xlab, ...)
      for (i in 1:length(metas)) {
        if (dispersion == "IQR") {
          .dispersion2(xcoords, metas[[i]], bound1[[i]], 
                       col = dispersion.col[i], ...)
          .dispersion2(xcoords, metas[[i]], llim = q1[[i]], 
                       ulim = q3[[i]], intervals = FALSE, col = dispersion.col[i], 
                       ...)
          if (all(metas[[i]] == q1[[i]])) {
            .dispersion2(xcoords, metas[[i]], llim = bound1[[i]], 
                         ulim = metas[[i]], col = dispersion.col[i], 
                         ...)
          }
          if (all(metas[[i]] == q3[[i]])) {
            .dispersion2(xcoords, metas[[i]], ulim = bound1[[i]], 
                         llim = metas[[i]], col = dispersion.col[i], 
                         ...)
          }
        }
        else {
          .dispersion2(xcoords, metas[[i]], bound2[[i]], 
                       col = dispersion.col[i], ...)
          .dispersion2(xcoords, metas[[i]], bound1[[i]], 
                       col = dispersion.col[i], ...)
        }
      }
      for (j in 1:length(metas)) {
        lines(xcoords, metas[[j]], col = line.col[j], 
              ...)
      }
    }
    else {
      
      cat("run this step\n")
      
      plot(xcoords, metas[[1]], type = "l", col = line.col[1], 
           ylim = myrange, ylab = ylab, xlab = xlab, ...)
      for (i in 2:length(metas)) {
        lines(xcoords, metas[[i]], col = line.col[i], 
              ...)
      }
      
      metas.L <- lapply(metas, function(u){ x <- u })
      col.L <- lapply(line.col, function(u){ x <- u })
      
      res <- list(xcoords = xcoords,metas = metas.L,col=col.L,ylim = myrange, ylab = ylab, xlab = xlab)
      res
      
    }
    if (!is.null(profile.names)) 
      if (!is.null(dispersion) && dispersion %in% disp.args) {
        legend(max(xcoords) + 0.05 * max(xcoords), myrange[2], 
               legend = profile.names, fill = dispersion.col, 
               bty = "n", border = line.col)
      }
    else {
      legend(max(xcoords) + 0.05 * max(xcoords), myrange[2], 
             legend = profile.names, fill = line.col, bty = "n")
    }
  }
  else {
    for (j in 1:length(metas)) {
      plot(xcoords, metas[[j]], type = "l", col = line.col[j], 
           ylim = myrange, ylab = ylab, xlab = xlab, ...)
      re <- list(xcoords=xcoords, metas=metas[[j]],col = line.col[j], 
                 ylim = myrange, ylab = ylab, xlab = xlab)
      if (!is.null(dispersion) && dispersion %in% disp.args) {
        if (dispersion == "IQR") {
          .dispersion2(xcoords, metas[[j]], bound1[[j]], 
                       col = dispersion.col[j], ...)
          .dispersion2(xcoords, metas[[j]], llim = metas[[j]] - 
                         q1[[j]], ulim = q3[[j]] - metas[[j]], col = dispersion.col[j], 
                       ...)
          if (all(metas[[i]] == q1[[i]])) {
            .dispersion2(xcoords, metas[[i]], llim = bound1[[i]], 
                         ulim = metas[[i]], col = dispersion.col[i], 
                         ...)
          }
          if (all(metas[[i]] == q3[[i]])) {
            .dispersion2(xcoords, metas[[i]], ulim = bound1[[i]], 
                         llim = metas[[i]], col = dispersion.col[i], 
                         ...)
          }
        }
        else {
          .dispersion2(xcoords, metas[[j]], bound1[[j]], 
                       col = dispersion.col[j], ...)
          .dispersion2(xcoords, metas[[j]], bound2[[j]], 
                       col = dispersion.col[j], ...)
        }
      }
      lines(xcoords, metas[[j]], col = line.col[j], ...)
    }
    re
  }
  par(xpd = FALSE)
  par(mar = marOrg)
  invisible(do.call("rbind", metas))
  res
}