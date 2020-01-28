# Fit a CTMC to dive data

#' Compute matrices required to fit temporal smooth
#'
#' @param forms formulae
#' @param dat data frame (see fitCTMCdive)
#' @param nint number of integration points
#'
#' @return list of mesh, SPDE objects, point A matrices, and integral A matrices
#' @export
#' @importFrom mgcv gam predict.gam
#' @importFrom methods as
MakeSmooth <- function(forms, dat, nint = 10000) {

  # results list
  res <- list()

  ## dive model
  # GAM setup
  gam_dive <- gam(forms[["dive"]], data = dat, method = "REML")
  res$dive_lambda <- 0
  # extract smoothing matrix
  if(length(gam_dive$smooth) > 0){
    res$S_dive <- as(gam_dive$smooth[[1]]$S[[1]], "sparseMatrix")
    res$dive_lambda <- gam_dive$sp
  }
  # build design matrix
  res$A_dive <- predict(gam_dive, newdata = data.frame(time = dat$time),
                    type = "lpmatrix")[, -1]
  # fixed effects design matrix
  res$Xs_dive <- model.matrix(gam_dive$pterms, data = dat)

  ## surface model
  # GAM setup
  gam_surface <- gam(forms[["surface"]], data = dat, method = "REML")
  res$surface_lambda <- 0
  # extract smoothing matrix
  if(length(gam_surface$smooth) > 0){
    res$S_surface <- as(gam_surface$smooth[[1]]$S[[1]], "sparseMatrix")
    res$surface_lambda <- gam_surface$sp
  }
  # build design matrix
  res$A_surf <- predict(gam_surface,
                        newdata = data.frame(time = dat$time + dat$dive),
                        type = "lpmatrix")[, -1]
  # fixed effects design matrix
  res$Xs_surface <- model.matrix(gam_surface$pterms, data = dat)

  # construct prediction grid
  n <- nrow(dat)
  ints <- seq(dat$time[1],
              dat$time[n] + dat$dive[n] + dat$surface[n],
              length = nint)
  # spacing on prediction grid
  dt <- mean(diff(ints))

  # predictor matrix
  res$A_grid_dive <- predict(gam_dive, newdata = data.frame(time = ints),
                        type = "lpmatrix")[, -1]
  res$A_grid_surface <- predict(gam_surface, newdata = data.frame(time = ints),
                        type = "lpmatrix")[, -1]

  # get integration matrices
  indD <- indS <- matrix(0, nrow = n, ncol = nint)
  dive_end <- dat$time + dat$dive
  for (i in 1:nint) {
    wh <- which(ints[i] <= dat$time[-1] & ints[i] >= dive_end[-n])
    if (length(wh) != 0) {
      indD[wh, i] <- 1
    }
    wh <- which(ints[i] >= dat$time & ints[i] <= dive_end)
    if (length(wh) != 0) {
      indS[wh, i] <- 1
    }
  }
  res$indD <- as(indD, "sparseMatrix")
  res$indS <- as(indS, "sparseMatrix")

  res$ints <- ints
  res$dt <- dt

  return(res)
}

#' Fits continuous-time Markov chain to dive and surface duration data
#'
#' @param forms a \code{list} with formulae for \code{dive} and \code{surface} variables
#' @param dat a \code{data.frame} with at least three columns named \code{dive} (dive durations), \code{surface} (surface durations), and \code{time} (start time of dive); all of these must be numeric.
#' @param print if \code{TRUE}, useful output is printed
#'
#' @return a CTMCdive model object: a list of the estimated results (\code{res}), variance matrix (\code{var}), fitted model returned from \code{optim} (\code{mod}), output from \code{sdreport} (\code{res}), design matrices (\code{Xs}), smoothing data (\code{sm}), formulae (\code{forms{), indices that divide par between dive and surface parameters (\code{len}), data (\code{dat}), indicators for smooths used (\code{lambda}), and model type (\code{model}).
#' @export
#' @importFrom stats delete.response model.matrix optim
#'             pnorm predict qnorm quantile terms
#' @importFrom TMB MakeADFun sdreport
#' @useDynLib ctmc_dive
FitCTMCdive <- function(forms, dat, print = TRUE) {
  ## Make design matrices and parameter vector
  if (print) cat("Computing design matrices.......")
  par <- NULL

  # name that list
  names(forms) <- c(as.character(forms[[1]][[2]]),
                    as.character(forms[[2]][[2]]))

  # fixed effects starting values
  par <- c(mean(dat$dive), mean(dat$surface))
  par <- -log(par)
  names(par) <- c("dive", "surface")

  # smoothing data
  random <- NULL
  map <- list()

  sm <- MakeSmooth(forms, dat)

  len <- c(ncol(sm$Xs_dive), ncol(sm$Xs_surface))
  names(len) <- c("dive", "surface")
  # all setup done now
  if(print) cat("done\n")

  ## Setup TMB parameter list
  tmb_parameters <- list(par_dive = par[["dive"]],
                         par_surf = par[["surface"]],
                         log_lambda_dive = 0,
                         log_lambda_surf = 0)

  # if there are smooths of dive or surface, set up the TMB
  # model objects correctly
  if (is.null(sm$S_dive)) {
    map <- c(map, list(log_lambda_dive = as.factor(NA),
                       s_dive = factor(NA)))
    tmb_parameters$s_dive <- 0
    sm$S_dive <- as(matrix(0, 1, 1), "sparseMatrix")
    sm$A_grid_dive <- matrix(1, length(sm$indD), 1)
  } else {
    random <- c(random, "s_dive")
    tmb_parameters$s_dive <- rep(0, ncol(sm$S_dive))
  }
  if (is.null(sm$S_surface)) {
    map <- c(map, list(log_lambda_surf = as.factor(NA),
                       s_surf = factor(NA)))
    tmb_parameters$s_surf <- 0
    sm$S_surface <- as(matrix(0, 1, 1), "sparseMatrix")
    sm$A_grid_surface <- matrix(1, length(sm$indS), 2)
  } else {
    random <- c(random, "s_surf")
    tmb_parameters$s_surf <- rep(0, ncol(sm$S_surface))
  }

  ## Setup TMB data list
  tmb_dat <- list(Xdive = sm$Xs_dive,
                  Xsurf = sm$Xs_surface,
                  S_dive = sm$S_dive,
                  S_surface = sm$S_surface,
                  A_dive = sm$A_dive,
                  A_surf = sm$A_surf,
                  A_grid_surface = sm$A_grid_surface,
                  A_grid_dive = sm$A_grid_dive,
                  indD = sm$indD,
                  indS = sm$indS,
                  dt = sm$dt)

  ## Create object
  obj <- MakeADFun(tmb_dat, tmb_parameters, random = random, map = map,
                   hessian = TRUE, DLL = "ctmc_dive")

  ## Fit Model
  if (print) cat("Fitting model.......\n")
  t0 <- Sys.time()
  mod <- do.call(optim, obj)
  t1 <- Sys.time()
  diff <- difftime(t1, t0)
  if(print) cat("Model fit in ", signif(diff[[1]], 2), attr(diff, "units"), "\n")

  ## Compute estimates
  if(print) cat("Estimating variance.......")
  npar <- sum(len)
  rep <- sdreport(obj)
  est <- rep$par.fixed
  var <- rep$cov.fixed
  sds <- sqrt(diag(var))
  if(print) cat("done\n")

  # confidence intervals
  if(print) cat("Computing confidence intervals.......")
  LCL <- est - qnorm(0.975) * sds
  UCL <- est + qnorm(0.975) * sds
  if(print) cat("done\n")

  # pvalues
  pval <- 2 * pnorm(-abs(est), 0, sds)

  #  dive result
  if(print) cat("Formatting results.......")
  s <- 1
  e <- len[1]
  res_dive <- data.frame(Est. = est[s:e],
                         SD. = sds[s:e],
                         LCL. = LCL[s:e],
                         UCL. = UCL[s:e],
                         p_value = pval[s:e])
  rownames(res_dive) <- colnames(sm$Xs_dive)

  # surface result
  s <- len[1] + 1
  e <- len[1] + len[2]
  res_surf <- data.frame(Est. = est[s:e],
                         SD. = sds[s:e],
                         LCL. = LCL[s:e],
                         UCL. = UCL[s:e],
                         p_value = pval[s:e])
  rownames(res_surf) <- colnames(sm$Xs_surface)


  # full result
  res <- list(surface = res_surf, dive = res_dive)

  ans <- list(res = res,
              var = var,
              mod = mod,
              rep = rep,
              Xs_dive = sm$Xs_dive,
              Xs_surface = sm$Xs_surface,
              sm = sm,
              forms = forms,
              len = len,
              dat = dat
             )
  if(print) cat("done\n")
  class(ans) <- "CTMCdive"
  return(ans)
}

#' Print summary of CTMCdive object
#'
#' @param object the CTMCdive fitted model object
#' @param \dots unused (for S3 compatability)
#'
#' @return prints a summary
#' @export
summary.CTMCdive <- function(object, ...) {
  cat("Continuous-time Markov Chain Diving Model\n")
  cat("Number of Dives: ", nrow(object$dat), "\n")
  cat("\n")
  cat("Model fit:\n")
  print(object$forms[[1]])
  print(object$forms[[2]])
  # TODO: rewrite this!!!
  #if (object$model != "iid") {
  #  cat("with cubic shrinkage smooth for ")
  #  if (object$lambda[1] > 1e-10) cat("dive intensity")
  #  if (object$lambda[1] > 1e-10 & object$lambda[2] > 1e-10) cat(" and ")
  #  if (object$lambda[2] > 1e-10) cat("surfacing intensity")
  #  cat(".")
  #}
  cat("\n")
  cat(rep("-", 30), "\n")
  cat("DIVE INTENSITY\n")
  print(signif(object$res$dive, 4))
  cat("\n")
  cat(rep("-", 30), "\n")
  cat("SURFACE INTENSITY\n")
  print(signif(object$res$surface, 4))
  cat("\n")
  invisible(object)
 }

print.CTMCdive <- summary.CTMCdive

#' Predict mean duration from fitted CTMC model for dives and surfacing
#'
#' @param object CTMCdive fitted model object
#' @param newdata  new data frame to make predictions for
#' @param \dots unused (for S3 compatability)
#'
#' @return a list of surface and dive mean predictions
#' @export
predict.CTMCdive <- function(object, newdata = NULL, ...) {
  par <- object$rep$par.fixed
  len <- object$len
  par_dive <- par[1:len[1]]
  par_surf <- par[(len[1] + 1):(len[1] + len[2])]

  # linear predictors
  nu_dive <- object$Xs_dive %*% par_dive
  nu_surf <- object$Xs_surface %*% par_surf

  # create time grid
  ints <- object$sm$ints
  nints <- length(ints)
  lambda_dive <- rep(0, nints)
  lambda_surf <- rep(0, nints)
  for (i in 1:(nrow(object$dat) - 1)) {
    lambda_dive[ints >= object$dat$time[i] & ints < object$dat$time[i + 1]] <- nu_dive[i]
    lambda_surf[ints >= object$dat$time[i] & ints < object$dat$time[i + 1]] <- nu_surf[i]
  }
  lambda_dive[ints > object$dat$time[nrow(object$dat)]] <- nu_dive[length(nu_dive)]
  lambda_surf[ints > object$dat$time[nrow(object$dat)]] <- nu_surf[length(nu_surf)]

  # TODO: this needs to be fixed!!
  if (any(names(object$rep$par.random) == "s_surf")) {
    x_surf <- object$rep$par.random[names(object$rep$par.random) == "s_surf"]
    lambda_surf <- lambda_surf + (object$sm$A_grid_surface %*% x_surf)
  }
  if (any(names(object$rep$par.random) == "s_dive")) {
    x_dive <- object$rep$par.random[names(object$rep$par.random) == "s_dive"]
    lambda_dive <- lambda_dive + (object$sm$A_grid_dive %*% x_dive)
  }
  lambda_dive <- exp(lambda_dive)
  lambda_surf <- exp(lambda_surf)

  # get expectations
  dt <- mean(diff(ints))
  exp_dive <- exp_surf <- rep(0, nrow(object$dat))
  r_dive <- r_surf <- rep(0, nrow(object$dat))
  for (i in 1:nrow(object$dat)) {
    Ldive <- cumsum(lambda_dive[ints > object$dat$time[i] + object$dat$dive[i]])
    Lsurf <- cumsum(lambda_surf[ints > (object$dat$time[i])])
    Sdive <- exp(-Ldive * dt)
    Ssurf <- exp(-Lsurf * dt)
    # E(dive) comes from surface intensity and vice-versa
    exp_dive[i] <- sum(Ssurf * dt)
    exp_surf[i] <- sum(Sdive * dt)
    # get dive p-value
    r_dive[i] <- Ssurf[max(floor(object$dat$dive[i] / dt), 1)]
    r_surf[i] <- Sdive[max(floor(object$dat$surface[i] / dt), 1)]
  }
  # return predictions
  res <- list(surface = exp_surf, dive = exp_dive, rdive = r_dive, rsurf = r_surf)
  return(res)
}

#' Plot fitted CTMCdive model
#'
#' Plot models fitted by \code{CTMCdive}.
#'
#' @param x fitted CTMCdive model object
#' @param quant the quantile to plot the observed data up to, for example if set to 0.95 then
#' only the durations up to the 95\% quantile are plotted (prevents outliers dominating the plot)
#' @param pick if not \code{NULL} then either "dive" or "surface" if only want a plot of one of the fitted
#' processes
#' @param pred if predicted means already computed can supply them here rather than have them recomputed
#' @param xlim limits for x axis if you want to restrict range
#' @param \dots unused (for S3 compatability)
#'
#' @return plots of dive and surface durations or only one of if "\code{pick}" is specified
#' @export
#' @importFrom graphics lines par plot
plot.CTMCdive <- function(x, quant = 1, pick = NULL, pred = NULL, xlim = NULL, ...) {
  if (is.null(pick)) {
    par(mfrow=c(2, 1))
    on.exit(par(mfrow=c(1,1)))
    pick <- "all"
  }
  # get predicted values
  if (is.null(pred)) pred <- predict(x)
  # get maximum time if scaled
  time <- x$dat$time
  # plot fitted values over observed
  if (pick == "all" | pick == "surface") {
    q <- quantile(x$dat$surface, prob = quant)
    plot(time, x$dat$surface, col = "grey60", xlab = "Time", ylab = "Surface duration",  ylim = c(min(x$dat$surface), q), xlim = xlim)
    lines(time, pred$surface, col = "red")
  }
  if (pick == "all" | pick == "dive") {
    q <- quantile(x$dat$dive, prob = quant)
    plot(time, x$dat$dive, col = "grey60", xlab = "Time", ylab = "Dive duration",  ylim = c(min(x$dat$dive), q), xlim = xlim)
    lines(time, pred$dive, col = "blue")
  }
  invisible(list(mod = x, pred = pred))
}

#' Returns log-likelihood with degrees of freedom
#'
#' @param object fitted CTMCdive object
#' @param \dots unused (for S3 compatability)
#'
#' @return logLik with df attribute
#' @note This does not account for degrees of freedom reduction with smooths (i.e
#' if lambda > 0)
#' @export
logLik.CTMCdive <- function(object, ...) {
  val <- -object$mod$value
  npar <- length(object$mod$par)
  llk <- val
  attributes(llk)$df <- npar
  return(llk)
}
