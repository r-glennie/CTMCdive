# Fit a CTMC to dive data

#' Compute matrices required to fit temporal smooth
#'
#' @param forms formulae
#' @param dat data frame (see fitCTMCdive)
#' @param nk number of knots (UNUSED)
#' @param nint number of integration points
#'
#' @return list of mesh, SPDE objects, point A matrices, and integral A matrices
#' @export
#' @importFrom mgcv gam predict.gam
#' @importFrom methods as
MakeSmooth <- function(forms, dat, nk = 100, nint = 10000) {
  # get gam data
  gamdat <- gam(dive ~ s(time, bs = "cs"), data = dat, method = "REML")
  # extract smoothing matrix
  S <- as(gamdat$smooth[[1]]$S[[1]], "sparseMatrix")

  # construct prediction grid
  n <- nrow(dat)
  ints <- seq(dat$time[1], dat$time[n] + dat$dive[n] + dat$surface[n],
              length = nint)
  # spacing on prediction grid
  dt <- mean(diff(ints))
  # predictor matrices
  A_dive <- predict(gamdat, newdata = data.frame(time = dat$time),
                    type = "lpmatrix")[,-1]
  A_surf <- predict(gamdat, newdata = data.frame(time = dat$time + dat$dive),
                    type = "lpmatrix")[,-1]
  A_grid <- predict(gamdat, newdata = data.frame(time = ints),
                    type = "lpmatrix")[,-1]

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
  indD <- as(indD, "sparseMatrix")
  indS <- as(indS, "sparseMatrix")
  res <- list(S = S,
              A_dive = A_dive,
              A_surf = A_surf,
              A_grid = A_grid,
              indD = indD,
              indS = indS,
              ints = ints,
              dt = dt)
  return(res)
}

#' Fits continuous-time Markov chain to dive and surface duration data
#'
#' @param forms a \code{list} with formulae for \code{dive} and \code{surface} variables
#' @param dat a \code{data.frame} with at least three columns named \code{dive} (dive durations), \code{surface} (surface durations), and \code{time} (start time of dive); all of these must be numeric.
#' @param model \code{"iid"} fits a CTMC where durations are independent over time, \code{"d"} fits a model where dive durations are correlated, \code{"s"} where surface durations are correlated, and \code{"ds"} fits where both dive and surface are correlated
#' @param print if \code{TRUE}, useful output is printed
#'
#' @return a CTMCdive model object: a list of the estimated results (\code{res}), variance matrix (\code{var}), fitted model returned from \code{optim} (\code{mod}), output from \code{sdreport} (\code{res}), design matrices (\code{Xs}), smoothing data (\code{sm}), formulae (\code{forms{), indices that divide par between dive and surface parameters (\code{len}), data (\code{dat}), indicators for smooths used (\code{lambda}), and model type (\code{model}).
#' @export
#' @importFrom stats delete.response model.matrix optim
#'             pnorm predict qnorm quantile terms
#' @importFrom TMB MakeADFun sdreport
#' @useDynLib ctmc_dive
FitCTMCdive <- function(forms, dat, model = "iid", print = TRUE) {
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
  par[1] <- -log(mean(dat$dive))
  par[len[1] + 1] <- -log(mean(dat$surface))
  if(print) cat("done\n")

  # smoothing data
  random <- NULL
  map <- list()

  sm <- MakeSmooth(forms, dat)
  n_mesh <- ncol(sm$S)
  lambda <- rep(0, 2)
  if (model != "iid") {
    if (model == "d") lambda[1] <- 1
    if (model == "s") lambda[2] <- 1
    if (model == "ds" | model == "sd") lambda[1:2] <- 1
  }

  ## Setup TMB parameter list
  tmb_parameters <- list(par_dive = par[1:len[1]],
                         par_surf = par[(len[1] + 1):length(par)],
                         log_lambda_dive = 0,
                         log_lambda_surf = 0,
                         s_dive = rep(0, n_mesh),
                         s_surf = rep(0, n_mesh))

  ## Create TMB model object
  if (!(lambda[1] < 1e-10)) {
    random <- c(random, "s_dive")
    tmb_parameters$s_dive <- rep(0, n_mesh)
  } else {
    map <- c(map, list(log_lambda_dive = as.factor(NA),
                       s_dive = factor(rep(NA, n_mesh))))
  }
  if (!(lambda[2] < 1e-10)) {
    random <- c(random, "s_surf")
    tmb_parameters$s_surf <- rep(0, n_mesh)
  } else {
    map <- c(map, list(log_lambda_surf = as.factor(NA),
                       s_surf = factor(rep(NA, n_mesh))))
  }

  ## Setup TMB data list
  tmb_dat <- list(Xdive = Xs[[1]],
                  Xsurf = Xs[[2]],
                  S = sm$S,
                  A_dive = sm$A_dive,
                  A_surf = sm$A_surf,
                  A_grid = sm$A_grid,
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
  rownames(res_dive) <- colnames(Xs[[1]])

  # surface result
  s <- len[1] + 1
  e <- len[1] + len[2]
  res_surf <- data.frame(Est. = est[s:e],
                         SD. = sds[s:e],
                         LCL. = LCL[s:e],
                         UCL. = UCL[s:e],
                         p_value = pval[s:e])
  rownames(res_surf) <- colnames(Xs[[2]])


  # full result
  res <- list(surface = res_surf, dive = res_dive)

  ans <- list(res = res,
              var = var,
              mod = mod,
              rep = rep,
              Xs = Xs,
              sm = sm,
              forms = forms,
              len = len,
              dat = dat,
              lambda = lambda,
              model = model)
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
  if (object$model != "iid") {
    cat("with cubic shrinkage smooth for ")
    if (object$lambda[1] > 1e-10) cat("dive intensity")
    if (object$lambda[1] > 1e-10 & object$lambda[2] > 1e-10) cat(" and ")
    if (object$lambda[2] > 1e-10) cat("surfacing intensity")
    cat(".")
  }
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
  Xs <- object$Xs
  # linear predictors
  nu_dive <- Xs[[1]] %*% par_dive
  nu_surf <- Xs[[2]] %*% par_surf
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
  # add RW smooths if necessary
  if (object$model != "iid") {
    npar <- sum(len)
    lambda <- object$lambda
    if (lambda[2]>1e-10) {
      x_surf <- object$rep$par.random[names(object$rep$par.random) == "s_surf"]
      lambda_surf <- lambda_surf + (object$sm$A_grid %*% x_surf)
    }
    if (lambda[1]>1e-10) {
      x_dive <- object$rep$par.random[names(object$rep$par.random) == "s_dive"]
      lambda_dive <- lambda_dive + (object$sm$A_grid %*% x_dive)
    }
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
    r_dive[i] <- Ssurf[floor(object$dat$dive[i] / dt)]
    r_surf[i] <- Sdive[floor(object$dat$surface[i] / dt)]
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
