# Fit a CTMC to dive data

#' Compute matrices required to fit temporal smooth
#'
#' @param dat data frame (see fitCTMCdive)
#' @param nknots number of basis functions / knots
#'
#' @return INLA finite element matrices required to fit smooth
#' @export
MakeSmooth <- function(dat, nknots = NULL) {
  # make mesh over time
  knot.locations <- seq(min(dat$time), max(dat$time), length = nknots + 1)
  mesh <- inla.mesh.1d(knot.locations, interval = c(min(dat$time), max(dat$time)), degree = 2)
  # make penalty matrices
  ind <- unique(dat$ID)
  nind <- length(ind)
  Cind <- vector(mode = "list", length = nind)
  G1ind <- vector(mode = "list", length = nind)
  G2ind <- vector(mode = "list", length = nind)
  Aind <- vector(mode = "list", length = nind)
  for (i in 1:nind) {
    S <- inla.mesh.1d.fem(mesh)
    Cind[[i]] <- S$c1
    G1ind[[i]] <- S$g1
    G2ind[[i]] <- S$g2
    Aind[[i]] <- inla.mesh.1d.A(mesh, dat$time)
  }
  res <- list(C = bdiag(Cind),
              G1 = bdiag(G1ind),
              G2 = bdiag(G2ind),
              A = bdiag(Aind))
  return(res)
}

#' Compute log-likelihood of continuous-time Markov chain
#'
#' @param par parameter vector on link scale
#' @param Xs design matrices
#' @param dat data frame (see fitCTMCdive)
#' @param len indices that divide par into dive, surfacing parameters
#' @param model type of model "iid" is no temporal correlation, "cor" is temporal correlation
#' @param model.args contains nknots (number of knots), smooth (output of MakeSmooth function), and smoothing
#' parameters (lambda)
#' @param print TRUE if useful output to be printed
#'
#' @return log-likelihood value
#' @export
CalcLlk <- function(par, Xs, dat, len, model = "iid", model.args = NULL, print = FALSE) {
  ## Compute parameters
  # linear predictors
  nu_surf <- Xs[[1]] %*% par[(1):(len[1])]
  nu_dive <- Xs[[2]] %*% par[(len[1] + 1):(len[1] + len[2])]
  # if random walk smooth added
  npar <- sum(len)
  if (model == "cor") {
    smooth <- model.args$S
    lambda <- model.args$lambda
    nknots <- model.args$nknots
    if (is.null(nknots)) nknots <- 10
    A <- smooth$A
    if (!is.na(lambda[1])) {
      x_surf_par <- par[(npar + 1):(npar + nknots)]
      x_surf <- A %*% x_surf_par
      nu_surf <- nu_surf + x_surf
      npar <- npar + nknots
    }
    if (!is.na(lambda[2])) {
      x_dive_par <- par[(npar + 1):(npar + nknots)]
      x_dive <- A %*% x_dive_par
      nu_dive <- nu_dive + x_dive
      npar <- npar + nknots
    }
  }
  # event intensity for diving and surfacing
  l_dive <- exp(nu_dive)
  l_surf <- exp(nu_surf)
  ## Compute likelihood
  llk <- 0
  # for observed dive times
  llk <- llk + sum(nu_surf - dat$dive * l_surf)
  # for observed surface times
  llk <- llk + sum(nu_dive - dat$surface * l_dive)
  # if RW smooth then add penalties
  if (model == "cor") {
    C <- smooth$C
    G1 <- smooth$G1
    G2 <- smooth$G2
    if (!is.na(lambda[1])) llk <- llk - as.numeric(lambda[1] * t(x_surf_par) %*% (lambda[1]^4 * C + 2 * lambda[1]^2 * G1 + G2) %*% x_surf_par)
    if (!is.na(lambda[2])) llk <- llk - as.numeric(lambda[2] * t(x_dive_par) %*% (lambda[2]^4 * C + 2 * lambda[2]^2 * G1 + G2) %*% x_dive_par)
  }
  ## Print likelihood if desired
  if (print) cat(" llk: ", llk, "\n")
  return(llk)
}

#' Compute negative log-likelihood
#'
#' @inheritParams CalcLlk
#'
#' @return negative log-likelihood
#' @export
CalcNegllk <- function(par, Xs, dat, len, model = "iid", model.args = NULL, print = FALSE) {
  llk <- CalcLlk(par, Xs, dat, len, model = model, model.args, print)
  return(-llk)
}

#' Fits continuous-time Markov chain to dive and surface duration data
#'
#' @param forms a list with formulae for "dive" and "surface" variables
#' @param dat a data frame with at least three columns named "dive" (dive durations), "surface" (surface durations),
#' and "time" (start time of dive); all of these must be numeric.
#' @param model "iid" fits a CTMC where durations are independent over time, "cor" fits temporal smooths so that
#' durations can be temporally correlated.
#' @param model.args if model = "cor" then model.args is a list with "nknots" set to the number of knots in the
#' temporal smooth; it also must contain a vector "lambda" with two elements that correspond to "dive" and "surface"
#' processes, if set to NA then no temporal smooth is fit, otherwise must be a positive number. This lambda is
#' the smoothing penalty: the larger lambda is, the smoother the estimated temporal curve.
#' @param print if TRUE, useful output is printed
#' @param iterlim iteration limit for nlm optimiser (default is 500)
#'
#' @return a CTMCdive model object: a list of the estimated results (res), variance matrix (var), fitted model returned from nlm (mod),
#' design matrices (Xs), formulae (forms), indices that divide par between dive and surface parameters (len),
#' data frame (dat), model arguments (model.args) including smoothing matrices, and the
#' model type (model).
#' @export
FitCTMCdive <- function(forms, dat, model = "iid", model.args = NULL, print = TRUE, iterlim = 500) {
  ## Make design matrices and parameter vector
  if (print) cat("Computing design matrices.......")
  Xs <- vector(mode = "list", length = 2)
  par <- NULL
  len <- rep(0, 2)
  b <- 0
  for (i in 1:2) {
    ind <- switch (as.character(forms[[i]][[2]]), "surface" = 2, "dive" = 1)
    ter <- terms(forms[[i]])
    ter <- delete.response(ter)
    Xs[[ind]] <- model.matrix(ter, data = dat)
    len[ind] <- ncol(Xs[[ind]])
    par <- c(par, rep(0, len[ind]))
    b <- b + len[ind]
  }
  par[1] <- 1.0 / mean(dat$dive)
  par[len[1] + 1] <- 1.0 / mean(dat$surface)
  if(print) cat("done\n")

  # starting par given?
  if (!is.null(model.args$ini_par)) par <- model.args$ini_par

  ## If model = "cor" then compute smooth information
  if (model == "cor") {
    if(print) cat("Computing smoothing matrices.......")
    nknots <- model.args$nknots
    if (is.null(nknots)) nknots <- 10
    model.args$S <- MakeSmooth(dat, nknots)
    # add parameters
    lambda <- model.args$lambda
    if (!is.na(lambda[1])) par <- c(par, rep(0, nknots))
    if (!is.na(lambda[2])) par <- c(par, rep(0, nknots))
    if(print) cat("done\n")
  }

  ## Fit model
  if (print) cat("Fitting model.......\n")
  t0 <- Sys.time()
  mod <- suppressWarnings(nlm(CalcNegllk,
                              par,
                              Xs = Xs,
                              dat = dat,
                              len = len,
                              model = model,
                              model.args = model.args,
                              print = print,
                              hessian = TRUE,
                              iterlim = iterlim))
  if (mod$code > 2) stop("Convergence failed with nlm code ", mod$code, ".\n")
  t1 <- Sys.time()
  diff <- difftime(t1, t0)
  if(print) cat("Model fit in ", signif(diff[[1]], 2), attr(diff, "units"), "\n")

  ## Compute estimates
  if(print) cat("Estimating variance.......")
  npar <- sum(len)
  est <- mod$estimate[1:npar]
  var <- try(solve(mod$hessian[1:npar, 1:npar]))
  if ("try-error" %in% class(var)) warning("Could not invert Hessian to compute variance matrix.")
  sds <- sqrt(diag(var))
  if(print) cat("done\n")

  # confidence intervals
  if(print) cat("Computing confidence intervals.......")
  LCL <- est - qnorm(0.975) * sds
  UCL <- est + qnorm(0.975) * sds
  if(print) cat("done\n")

  # pvalues
  pval <- 2 * pnorm(-abs(est), 0, sds)

   # surface result
  if(print) cat("Formatting results.......")
  s <- 1
  e <- len[1]
  res_surf <- data.frame(Est. = est[s:e],
                          SD. = sds[s:e],
                          LCL. = LCL[s:e],
                          UCL. = UCL[s:e],
                          p_value = pval[s:e])
  res_surf <- signif(res_surf, 4)
  rownames(res_surf) <- colnames(Xs[[1]])

  # dive result
  s <- len[1] + 1
  e <- len[1] + len[2]
  res_dive <- data.frame(Est. = est[s:e],
                         SD. = sds[s:e],
                         LCL. = LCL[s:e],
                         UCL. = UCL[s:e],
                         p_value = pval[s:e])
  res_dive <- signif(res_dive, 4)
  rownames(res_dive) <- colnames(Xs[[2]])

  # full result
  res <- list(surface = res_surf, dive = res_dive)

  ans <- list(res = res,
              var = var,
              mod = mod,
              Xs = Xs,
              forms = forms,
              len = len,
              dat = dat,
              model.args = model.args,
              model = model)
  if(print) cat("done\n")
  class(ans) <- "CTMCdive"
  return(ans)
}

#' Print summary of CTMCdive object
#'
#' @param mod the CTMCdive fitted model object
#'
#' @return prints a summary
#' @export
summary.CTMCdive <- function(mod) {
  cat("Continuous-time Markov Chain Diving Model\n")
  cat("Number of Dives: ", nrow(mod$dat), "\n")
  cat("\n")
  cat("Model fit:\n")
  print(mod$forms[[1]])
  print(mod$forms[[2]])
  if (mod$model == "cor") {
    nknots <- ifelse(is.null(mod$model.args$nknots), 10, mod$model.args$nknots)
    cat("with 1D SPDE smooth with", nknots, "knots for ")
    if (!is.na(mod$model.args$lambda[1])) cat("dive duration")
    if (!is.na(mod$model.args$lambda[1]) & !is.na(mod$model.args$lambda[2])) cat(" and ")
    if (!is.na(mod$model.args$lambda[2])) cat("surfacing duration")
    cat(".")
  }
  cat("\n")
  cat(rep("-", 30), "\n")
  cat("DIVE DURATION\n")
  print(mod$res$surface)
  cat("\n")
  cat(rep("-", 30), "\n")
  cat("SURFACE DURATION\n")
  print(mod$res$dive)
  cat("\n")
  invisible(mod)
 }

print.CTMCdive <- summary.CTMCdive

#' Predict mean duration from fitted CTMC model for dives and surfacing
#'
#' @param mod CTMCdive fitted model object
#' @param newdata  new data frame to make predictions for
#'
#' @return a list of surface and dive mean predictions
#' @export
predict.CTMCdive <- function(mod, newdata = NULL) {
  par <- mod$mod$estimate
  # get design matrices
  if (!is.null(newdata)) {
    forms <- mod$forms
    Xs <- vector(mode = "list", length = 2)
    for (i in 1:2) {
      ind <- switch (as.character(forms[[i]][[2]]), "surface" = 2, "dive" = 1)
      ter <- terms(forms[[i]])
      ter <- delete.response(ter)
      Xs[[ind]] <- model.matrix(ter, data = newdata)
    }
  } else {
    Xs <- mod$Xs
  }
  # linear predictors
  len <- mod$len
  nu_surf <- Xs[[1]] %*% par[(1):(len[1])]
  nu_dive <- Xs[[2]] %*% par[(len[1] + 1):(len[1] + len[2])]
  # add RW smooths if necessary
  if (mod$model == "cor") {
    npar <- sum(len)
    smooth <- mod$model.args$S
    lambda <- mod$model.args$lambda
    nknots <- mod$model.args$nknots
    if (is.null(nknots)) nknots <- 10
    if (!is.na(lambda[1])) {
      x_surf <- par[(npar + 1):(npar + nknots)]
      nu_surf <- nu_surf + smooth$A %*% x_surf
      npar <- npar + nknots
    }
    if (!is.na(lambda[2])) {
      x_dive <- + par[(npar + 1):(npar + nknots)]
      nu_dive <- nu_dive + smooth$A %*% x_dive
      npar <- npar + nknots
    }
  }
  # predicted durations for diving and surfacing (1/intensity)
  mu_dive <- exp(-nu_surf)
  mu_surf <- exp(-nu_dive)

  # return predictions
  res <- list(surface = mu_surf, dive = mu_dive)
}

#' Plot fitted CTMCdive model
#'
#' @param mod fitted CTMCdive model object
#' @param quant the quantile to plot the observed data up to, for example if set to 0.95 then
#' only the durations upto the 95% quantile are plotted (prevents outliers dominating the plot)
#' @param pick if not NULL then either "dive" or "surface" if only want a plot of one of the fitted
#' processes
#' @param pred if predicted means already computed can supply them here rather than have them recomputed
#' @param xlim limits for x axis if you want to restrict range
#'
#' @return plots of dive and surface durations or only one of if "pick" is used
#' @export
plot.CTMCdive <- function(mod, quant = 1, pick = NULL, pred = NULL, xlim = NULL) {
  if (is.null(pick)) {
    par(mfrow=c(2, 1))
    on.exit(par(mfrow=c(1,1)))
    pick <- "all"
  }
  # get predicted values
  if (is.null(pred)) pred <- predict(mod)
  # get maximum time if scaled
  max.time <- 1
  if (!is.null(attributes(mod$dat)$max.time)) max.time <- attributes(mod$dat)$max.time
  time <- dat$time * max.time
  # plot fitted values over observed
  if (pick == "all" | pick == "surface") {
    q <- quantile(dat$surface, prob = quant)
    plot(time, dat$surface, col = "grey60", xlab = "Time", ylab = "Surface duration",  ylim = c(min(dat$surface), q), xlim = xlim)
    lines(time, pred$surface, col = "red")
  }
  if (pick == "all" | pick == "dive") {
    q <- quantile(dat$dive, prob = quant)
    plot(time, dat$dive, col = "grey60", xlab = "Time", ylab = "Dive duration",  ylim = c(min(dat$dive), q), xlim = xlim)
    lines(time, pred$dive, col = "blue")
  }
  invisible(list(mod = mod, pred = pred))
}

#' Returns log-likelihood with degrees of freedom
#'
#' @param object fitted CTMCdive object
#' @param ...
#'
#' @return logLik with df attribute
#' @note This does not account for degrees of freedom reduction with smooths (i.e
#' if lambda > 0)
#' @export
logLik.CTMCdive <- function(object, ...) {
  val <- -object$mod$minimum
  npar <- length(object$mod$estimate)
  llk <- val
  attributes(llk)$df <- npar
  return(llk)
}

#' Simulate dive and surface durations for a single individual
#'
#' @param forms formulae for dive and surface processes
#' @param dat data frame with covariates
#' @param par true parameter values for linear predictors with "dive" then "surface"
#' parameters
#' @param n number of records to simulate
#' @param model see fitCTMCdive
#' @param model.args seefitCTMCdive
#' @param print see fitCTMCdive
#'
#' @return data frame with "dive", "surface", and time columns added
#' @export
simulateCTMCdive <- function(forms, dat, par, n = 1, model = "iid", model.args = NULL, print = TRUE) {
  ## Make design matrices and parameter vector
  Xs <- vector(mode = "list", length = 2)
  len <- rep(0, 2)
  for (i in 1:2) {
    ind <- switch (as.character(forms[[i]][[2]]), "surface" = 2, "dive" = 1)
    ter <- terms(forms[[i]])
    ter <- delete.response(ter)
    Xs[[ind]] <- model.matrix(ter, data = dat)
    len[ind] <- ncol(Xs[[ind]])
  }
  ## If model = "cor" then compute smooth information
  if (model == "cor") {
    nknots <- model.args$nknots
    if (is.null(nknots)) nknots <- 10
    model.args$S <- MakeSmooth(dat, nknots)
    # add parameters
    lambda <- model.args$lambda
    if (!is.na(lambda[1])) par <- c(par, rep(0, nknots))
    if (!is.na(lambda[2])) par <- c(par, rep(0, nknots))
  }
  # linear predictors
  nu_surf <- Xs[[1]] %*% par[(1):(len[1])]
  nu_dive <- Xs[[2]] %*% par[(len[1] + 1):(len[1] + len[2])]
  # event intensity for diving and surfacing
  l_dive <- exp(nu_dive)
  l_surf <- exp(nu_surf)
  # simulate dives and surfacings
  dive <- rexp(n, l_dive)
  surface <- rexp(n, l_surf)
  time <- cumsum(dive + surface) - dive[1] - surface[1]
  res <- data.frame(time = time, dive = dive, surface = surface)
  return(cbind(res, dat))
}


