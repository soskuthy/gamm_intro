post.sim.gamm <- function (mod, newdat, n=10000, exclude=NULL) {
  lp.full <- predict(mod, newdata=newdat, type = "lpmatrix", exclude=exclude)
  coefs.full <- coef(mod)
  vc.full <- vcov(mod)
  sim.full <- mvrnorm(n, mu = coefs.full, Sigma = vc.full)
  fits.full <- lp.full %*% t(sim.full)
}

summary.lines <- function (mod, start, end) {
  s <- capture.output(summary(mod))
  cat(paste(s[start:end], collapse="\n"))
}

summary.coefs <- function (mod, digits=max(3, getOption("digits") - 3), ...) {
  s <- capture.output(print(summary(mod), digits=digits))
  start <- grep("Parametric coefficients", s)
  sigs <- grep("Signif[.] codes", s)
  end <- sigs[length(sigs)]
  cat(paste(s[start:end][-grep("(Signif[.] codes|---)", s[start:end])], collapse="\n"))
}


get_difference <- function (model, comp, cond = NULL, rm.ranef = NULL, se = TRUE, 
          f = 1.96, print.summary = getOption("itsadug_print")) 
{
  newd <- NULL
  su <- model$var.summary
  dat <- model$model
  if (is.null(names(comp))) {
    stop("Predictor specified in 'comp' unknown. Please provide a named list for 'comp', in the form of 'comp=list(Predictor=c('level1', 'level2'))'.")
  }
  if (all(names(comp) %in% colnames(dat))) {
    for (i in 1:length(comp)) {
      if (length(comp[[i]]) < 2) {
        stop(sprintf("Provide two levels for %s to calculate difference.", 
                     names(comp)[i]))
      } else if (length(comp[[i]]) > 2) {
        warning(sprintf("More than two levels provided for predictor %s. Only first two levels are being used.", 
                        names(comp)[i]))
      }
    }
  } else {
    errname <- paste(which(!names(comp) %in% colnames(dat)), 
                     collapse = ", ")
    stop(sprintf("Grouping predictor(s) not found in model: %s.", 
                 errname))
  }
  if (any(names(cond) %in% names(comp))) {
    for (i in names(cond)[names(cond) %in% names(comp)]) {
      cond[[i]] <- NULL
      warning(sprintf("Predictor %s specified in comp and cond. (The value in cond will be ignored.)", 
                      i))
    }
  }
  new.cond1 <- list()
  new.cond2 <- list()
  for (i in names(su)) {
    if (i %in% names(comp)) {
      new.cond1[[i]] <- comp[[i]][1]
      new.cond2[[i]] <- comp[[i]][2]
    } else if (i %in% names(cond)) {
      new.cond1[[i]] <- new.cond2[[i]] <- cond[[i]]
    } else {
      if (class(su[[i]]) == "factor") {
        new.cond1[[i]] <- as.character(su[[i]][1])
        new.cond2[[i]] <- as.character(su[[i]][1])
      } else if (class(su[[i]]) == "numeric") {
        new.cond1[[i]] <- su[[i]][2]
        new.cond2[[i]] <- su[[i]][2]
      }
    }
  }
  newd1 <- expand.grid(new.cond1)
  newd2 <- expand.grid(new.cond2)
  p1 <- predict(model, newd1, type = "lpmatrix")
  p2 <- predict(model, newd2, type = "lpmatrix")
  newd <- as.data.frame(newd1)
  newd.names <- colnames(newd)
  for (nn in newd.names) {
    if (nn %in% names(comp)) {
      newd[, nn] <- NULL
    }
  }
  mysummary <- summary_data(newd, print = FALSE)
  if (class(rm.ranef) == "logical") {
    if (rm.ranef[1] == FALSE) {
      rm.ranef <- NULL
    }
  }
  if (!is.null(rm.ranef)) {
    smoothlabels.table <- as.data.frame(do.call("rbind", 
                                                lapply(model$smooth, function(x) {
                                                  data.frame(Label = x[["label"]], Dim = x[["null.space.dim"]], 
                                                             Class = attr(x, "class")[1], stringsAsFactors = FALSE)
                                                })))
    smoothlabels <- as.vector(smoothlabels.table[smoothlabels.table$Class %in% 
                                                   c("random.effect", "fs.interaction"), "Label"])
    if (class(rm.ranef) == "logical") {
      if (rm.ranef[1] == TRUE) {
        rm.ranef <- smoothlabels
      } else {
        rm.ranef <- ""
      }
    } else if (inherits(rm.ranef, c("numeric", "integer"))) {
      smoothlabels.table <- smoothlabels.table[rm.ranef, 
                                               ]
      smoothlabels <- as.vector(smoothlabels.table[smoothlabels.table$Class %in% 
                                                     c("random.effect", "fs.interaction"), "Label"])
    }
    rm.col <- unlist(lapply(rm.ranef, function(x) {
      colnames(p1)[grepl(x, colnames(p1), fixed = TRUE)]
    }))
    rm.col <- unlist(lapply(smoothlabels, function(x) {
      rm.col[grepl(x, rm.col, fixed = TRUE)]
    }))
    p1[, rm.col] <- 0
    p2[, rm.col] <- 0
    predictors <- do.call("rbind", lapply(model$smooth, 
                                          function(x) {
                                            data.frame(Label = x[["label"]], Terms = x[["term"]])
                                          }))
    test <- table(predictors$Terms) - table(predictors[predictors$Label %in% 
                                                         rm.ranef, ]$Terms)
    for (pred in names(test[test == 0])) {
      if (pred %in% names(mysummary)) {
        mysummary[[pred]] <- paste(mysummary[[pred]], 
                                   "(Might be canceled as random effect, check below.)")
      }
    }
    if (length(rm.col) > 0) {
      mysummary[["NOTE"]] = sprintf("The following random effects columns are canceled: %s\n", 
                                    paste(smoothlabels, collapse = ","))
    } else {
      warning("No random effects to cancel.\n")
    }
  }
  p <- p1 - p2
  newd$difference <- as.vector(p %*% coef(model))
  if (se) {
    newd$CI <- f * sqrt(rowSums((p %*% vcov(model)) * 
                                  p))
  }
  if (print.summary == TRUE) {
    print_summary(mysummary)
  }
  return(newd)
}





plot_diff <- function (model, view, comp, cond = NULL, plotCI = TRUE, f = 1.96, 
          eegAxis = FALSE, col = "black", shade = TRUE, n.grid = 100, 
          add = FALSE, print.summary = getOption("itsadug_print"), 
          plot = TRUE, rm.ranef = NULL, main = NULL, ylab = NULL, xlab = NULL, 
          xlim = NULL, ylim = NULL, transform.view = NULL, mark.diff = TRUE, 
          hide.label = FALSE, ...) 
{
  dat = model$model
  xvar <- NULL
  by_predictor <- NULL
  if (length(view) > 1) {
    warning("Only first element of 'view' is being used. Use plot_diff2 for plotting difference surfaces.")
  }
  else {
    xvar <- view[1]
    if (xvar %in% names(cond)) {
      warning(sprintf("Predictor %s specified in view and cond. Values in cond being used, rather than the whole range of %s.", 
                      xvar))
    } else {
      cond[[xvar]] <- seq(min(na.exclude(dat[, xvar])), 
                          max(na.exclude(dat[, xvar])), length = n.grid)
    }
  }
  if (!is.null(xlim)) {
    if (length(xlim) != 2) {
      warning("Invalid xlim values specified. Argument xlim is being ignored.")
    }
    else {
      cond[[xvar]] <- seq(xlim[1], xlim[2], length = n.grid)
    }
  }
  newd <- c()
  newd <- get_difference(model, comp = comp, cond = cond, print.summary = print.summary, 
                         rm.ranef = rm.ranef, f = f)
  errormessage <- function() {
    return("Error: the function specified in transformation.view cannot be applied to x-values, because infinite or missing values are not allowed.")
  }
  if (!is.null(transform.view)) {
    tryCatch(newd[, xvar] <- sapply(newd[, xvar], transform.view), 
             error = function(x) {
             }, warning = function(x) {
             })
    if (any(is.infinite(newd[, xvar])) | any(is.nan(newd[, 
                                                         xvar])) | any(is.na(newd[, xvar]))) {
      stop(errormessage())
    }
    if (print.summary) {
      cat("\t* Note: x-values are transformed.\n")
    }
  }
  if (is.null(main)) {
    levels1 <- paste(sapply(comp, function(x) x[1]), collapse = ".")
    levels2 <- paste(sapply(comp, function(x) x[2]), collapse = ".")
    main = sprintf("Difference between %s and %s", levels1, 
                   levels2)
  }
  if (is.null(ylab)) {
    ylab = sprintf("Est. difference in %s", as.character(model$formula[[2]]))
  }
  if (is.null(xlab)) {
    xlab = xvar
  }
  if (is.null(ylim)) {
    ylim <- range(newd$difference)
    if (plotCI) {
      ylim <- with(newd, range(c(difference + CI, difference - 
                                   CI)))
    }
  }
  out <- data.frame(est = newd$difference, x = newd[, xvar])
  names(out)[2] <- xvar
  if (plotCI) {
    out$CI <- newd$CI
    out$f <- f
  }
  out$comp = list2str(names(comp), comp)
  if (plot == TRUE) {
    if (add == FALSE) {
      emptyPlot(range(newd[, xvar]), ylim, main = main, 
                xlab = xlab, ylab = ylab, h0 = 0, eegAxis = eegAxis, 
                ...)
      if (hide.label == FALSE) {
        addlabel = "difference"
        if (!is.null(rm.ranef)) {
          if (rm.ranef != FALSE) {
            addlabel = paste(addlabel, "excl. random", 
                             sep = ", ")
          }
        }
        mtext(addlabel, side = 4, line = 0, adj = 0, 
              cex = 0.75, col = "gray35", xpd = TRUE)
      }
    }
    if (plotCI == TRUE) {
      plot_error(newd[, xvar], newd$difference, newd$CI, 
                 shade = shade, col = col, ...)
    }
    else {
      lines(newd[, xvar], newd$difference, col = col, ...)
    }
    if (mark.diff == TRUE) {
      diff <- find_difference(newd$difference, newd$CI, 
                              newd[, xvar])
      if (length(diff$start) > 0) {
        addInterval(pos = getFigCoords("p")[3], diff$start, 
                    diff$end, col = "red", lwd = 2 * par()$lwd, 
                    length = 0, xpd = TRUE)
        abline(v = c(diff$start, diff$end), lty = 3, 
               col = "red")
      }
    }
    if (print.summary) {
      if (length(diff$start) > 0) {
        tmp <- c(sprintf("%s window(s) of significant difference(s):", 
                         xvar), sprintf("\t%f - %f", diff$start, diff$end))
      }
      else {
        tmp <- "Difference is not significant."
      }
      cat("\n")
      cat(paste(tmp, collapse = "\n"))
      cat("\n")
    }
    invisible(out)
  }
  else {
    return(out)
  }
}

plot_smooth.cont <- function (x, plot_all.c = NULL, cond = list(), ylim=NULL, add=F, ...) {
  dat <- x$model
  levels <- unique(x$model[,plot_all.c])
  cols <- rainbow(length(levels))
  if (is.null(ylim) & add==F) {
    fvs <- list()
    for (l in 1:length(levels)) {
      cnd <- list(levels[l])
      names(cnd) <- plot_all.c
      cond2 <- cond
      cond2[[plot_all.c]] <- levels[l]
      pdf(file=NULL)
      fvs[[l]] <- plot_smooth(x, cond=cond2, print.summary=F, ...)$fv
      dev.off()
    }
    ymin <- min(do.call(rbind, fvs)$ll)
    ymax <- max(do.call(rbind, fvs)$ul)
    ylim <- c(ymin - (ymax - ymin)*0.05, ymax + (ymax - ymin)*0.05)
  }
  cond2 <- cond
  cond2[[plot_all.c]] <- levels[1]
  plot_smooth(x, cond=cond2, ylim=ylim, col=cols[1], print.summary=F, add=add, ...)
  mtext(as.character(levels[1]), 3, line=1, adj=1, col=cols[1])
  for (l in 2:length(levels)) {
    cond2 <- cond
    cond2[[plot_all.c]] <- levels[l]
    plot_smooth(x, cond=cond2, add=T, col=cols[l], , print.summary=F, ...)
    mtext(as.character(levels[l]), 3, line=2-l, adj=1, col=cols[l])
  }
}

