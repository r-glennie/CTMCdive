# Fit a CTMC to dive data

#' Compute matrices required to fit the model
#'
#' @param forms formulae
#' @param dat data frame (see fitCTMCdive)
#' @param min_dwell start point for the integrals
#' @param series TRUE if data is a series rather than dive-by-dive
#' @param nint number of integration points
#'
#' @return list of matrices required to fit the model
#' @export
#' @importFrom mgcv gam predict.gam
#' @importFrom methods as
#' @importFrom Matrix bdiag
MakeMatrices <- function(forms, dat, min_dwell, series = FALSE, nint = 10000) {

  # results list
  res <- list()

  ## dive model
  # GAM setup
  gam_dive <- gam(forms[["dive"]], data = dat, method = "REML")
  res$gam_dive <- gam_dive
  # accumulate smoothing matrix, block diagonal
  if(length(gam_dive$smooth) > 0){
    S_dive_list <- list()
    for(i in seq_along(gam_dive$smooth)){
      smoo <- gam_dive$smooth[[i]]
      for(j in seq_along(smoo$S)){
        S_dive_list <- c(S_dive_list, as(smoo$S[[j]], "sparseMatrix"))
        res$S_dive_n <- c(res$S_dive_n, nrow(smoo$S[[j]]))
        res$s_dive_k <- c(res$s_dive_k, smoo$bs.dim)
        res$s_dive_names <- c(res$s_dive_names, attr(smoo$sp, "names"))
      }
    }
    # build a block diagonal matrix
    res$S_dive <- bdiag(S_dive_list)
  }

  # reduce data if series to dive-by-dive
  if (series) {
    divedat <- dat[dat$start == 1,] 
  } else {
    divedat <- dat
  }
  
  ## build design matrix
  mm <- predict(gam_dive, newdata = divedat, type = "lpmatrix")
  # fixed effects design matrix
  res$Xs_dive <- mm[, 1:gam_dive$nsdf, drop=FALSE]
  # smooth design matrix
  res$A_dive <- mm[, -c(1:gam_dive$nsdf), drop=FALSE]
  # weibull design matrix 
  res$W_dive <- matrix(c(0, log(divedat$surface[-1] + 1e-10)), nr = nrow(divedat), nc = 1)
  
  ## surface model
  # GAM setup
  gam_surface <- gam(forms[["surface"]], data = dat, method = "REML")
  res$gam_surface <- gam_surface 
  # accumulate smoothing matrix, block diagonal
  if(length(gam_surface$smooth) > 0){
    S_surface_list <- list()
    for(i in seq_along(gam_surface$smooth)){
      smoo <- gam_surface$smooth[[i]]
      for(j in seq_along(smoo$S)){
        S_surface_list <- c(S_surface_list, as(smoo$S[[j]], "sparseMatrix"))
        res$S_surface_n <- c(res$S_surface_n, nrow(smoo$S[[j]]))
        res$s_surface_k <- c(res$s_surface_k, smoo$bs.dim)
        res$s_surface_names <- c(res$s_surface_names, attr(smoo$sp, "names"))
      }
    }
    # build a block diagonal matrix
    res$S_surface <- bdiag(S_surface_list)
  }
  if (series) {
    surfdat <- dat[dat$start == 1,]
  } else {
    surfdat <- dat
  }
  surfdat$time <- surfdat$time + surfdat$dive

  ## build design matrix
  mm <- predict(gam_surface, newdata = surfdat, type = "lpmatrix")
  # fixed effects design matrix
  res$Xs_surface <- mm[, 1:gam_surface$nsdf, drop=FALSE]
  # smooth design matrix
  res$A_surf <- mm[, -c(1:gam_surface$nsdf), drop=FALSE]
  # weibull design matrix
  res$W_surf <- matrix(log(surfdat$dive + 1e-10), nr = nrow(surfdat), nc = 1)
  
  # construct prediction grid
  n <- nrow(divedat)
  ints <- seq(divedat$time[1],
              divedat$time[n] + divedat$dive[n] + divedat$surface[n],
              length = nint)
  # spacing on prediction grid
  dt <- mean(diff(ints))

  # construct predictor matrix
  pred_mat_maker <- function(model, dat, ints){
    # storage
    newdat <- data.frame(time=ints)

    # what covars do we need?
    covs <- as.character(attr(delete.response(terms(model)),
                              "variables"))[-1]
    covs <- covs[covs!="time"]
    # as.factor shennanigans
    covs <- gsub("as.factor\\(", "", covs)
    covs <- gsub("\\)", "", covs)
    # for these covars, do an interpolation
    for(cc in covs){
      if(is.factor(dat[[cc]])){
        newdat[[cc]] <- approx(dat$time, dat[[cc]], ints, method="constant", rule=2)$y
        newdat[[cc]] <- factor(newdat[[cc]],
                                 1:length(unique(dat[[cc]])),
                                 levels(dat[[cc]]))
      }else{
        newdat[[cc]] <- approx(dat$time, dat[[cc]], ints, method="linear", rule=2)$y
      }
    }
    # build the Lp matrix!
    predict(model, newdata = newdat, type = "lpmatrix")
  }

  # do this for each part of the model
  pred_mat_dive <- pred_mat_maker(gam_dive, dat, ints)
  res$Xs_grid_dive <- pred_mat_dive[, 1:gam_dive$nsdf, drop=FALSE]
  res$A_grid_dive <- pred_mat_dive[, -c(1:gam_dive$nsdf), drop=FALSE]

  # weibull for integration grid 
  t_dive <- c(-1e-10, divedat$time + divedat$dive, max(divedat$time + divedat$surface + divedat$dive))
  br <- sort(unique(t_dive))
  cut_dive <- as.numeric(cut(ints, breaks = br))
  dt_dive <- ints - br[cut_dive]
  res$W_grid_dive <- log(dt_dive + 1e-10)
  
  pred_mat_surface <- pred_mat_maker(gam_surface, dat, ints)
  res$Xs_grid_surface <- pred_mat_surface[, 1:gam_surface$nsdf, drop=FALSE]
  res$A_grid_surface <- pred_mat_surface[, -c(1:gam_surface$nsdf), drop=FALSE]

  # weibull for integration grid 
  t_surf <- c(-1e-10, divedat$time, max(divedat$time + divedat$surface + divedat$dive))
  br <- sort(unique(t_surf))
  cut_surf <- as.numeric(cut(ints, breaks = br))
  dt_surf <- ints - br[cut_surf] 
  res$W_grid_surf <- log(dt_surf + 1e-10)
  
  # get integration matrices
  indD <- indS <- matrix(0, nrow = n, ncol = nint)
  dive_end <- divedat$time + divedat$dive
  for (i in 1:nint) {
    wh <- which((ints[i] < divedat$time[-1] + 1e-10) &
                (ints[i] > (dive_end[-n]+min_dwell$surface) - 1e-10))
    if (length(wh) != 0) {
      indD[wh, i] <- 1
    }
    wh <- which(ints[i] > (divedat$time+min_dwell$dive-1e-10) & ints[i] < (dive_end + 1e-10))
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
#' @param min_dwell Minimum dwell time in a state. Useful if, for example, dives are only categorised as such if they are longer than a certain interval. Named list, needs to be in the same units as \code{time}, \code{surface} and \code{dive}.
#' @param series if TRUE then series data, otherwise dive-by-dive data 
#' @param rf if TRUE then fit a discrete random effect on dive 
#' @param dt set time-step in integration, otherwise set to total_time / 10000 
#' 
#' @return a CTMCdive model x: a list of the estimated results (\code{res}), variance matrix (\code{var}), fitted model returned from \code{optim} (\code{mod}), output from \code{sdreport} (\code{res}), design matrices (\code{Xs}), smoothing data (\code{sm}), formulae (\code{forms{), indices that divide par between dive and surface parameters (\code{len}), data (\code{dat}), indicators for smooths used (\code{lambda}), and model type (\code{model}).
#' @export
#' @importFrom stats delete.response model.matrix optim
#'             pnorm predict qnorm quantile terms
#' @importFrom TMB MakeADFun sdreport
#' @useDynLib ctmc_dive
FitCTMCdive <- function(forms, dat, print = TRUE,
                        min_dwell=list(dive=0, surface=0), 
                        series = FALSE, 
                        re = "none",
                        dt = NULL) {

  ## Check series
  if ("start" %in% colnames(dat) & !series) warning("Did you mean to fit a series model? Series = FALSE but dat has a start column.")
  
  ## Make design matrices and parameter vector
  if (print) cat("Computing design matrices.......")
  
  ## Set time step
  if (is.null(dt)) {
    nint <- 10000
  } else {
    nint <- round(max(dat$time) / dt)
  }
  

  # name that list
  names(forms) <- c(as.character(forms[[1]][[2]]),
                    as.character(forms[[2]][[2]]))

  # smoothing data
  random <- NULL
  map <- list()
  sm <- MakeMatrices(forms, dat, min_dwell = min_dwell, series = series, nint = nint)

  len <- c(ncol(sm$Xs_dive), ncol(sm$Xs_surface))
  names(len) <- c("dive", "surface")

  # fixed effects starting values
  par_dive <- c(-log(mean(dat$dive)), rep(0, len[1]-1))
  par_surf <- c(-log(mean(dat$surface)), rep(0, len[2]-1))

  # all setup done now
  if(print) cat("done\n")

  ## Setup TMB parameter list
  tmb_parameters <- list(par_dive = par_dive,
                         par_surf = par_surf,
                         log_kappa_dive = 0, 
                         log_kappa_surf = 0, 
                         log_lambda_dive = rep(0, length(sm$S_dive_n)),
                         log_lambda_surf = rep(0, length(sm$S_surface_n)), 
                         log_rf_sd = rep(0, 2))

  # if there are smooths of dive or surface, set up the TMB
  # model xs correctly
  if (is.null(sm$S_dive)) {
    map <- c(map, list(log_lambda_dive = as.factor(NA),
                       s_dive = factor(NA)))
    tmb_parameters$s_dive <- 0
    tmb_parameters$log_lambda_dive <- 0 
    sm$S_dive <- as(matrix(0, 1, 1), "sparseMatrix")
    sm$S_dive_n <- 0 
    sm$A_grid_dive <- matrix(1, ncol(sm$indD), 1)
  } else {
    random <- c(random, "s_dive")
    tmb_parameters$s_dive <- rep(0, ncol(sm$S_dive))
  }
  if (is.null(sm$S_surface)) {
    map <- c(map, list(log_lambda_surf = as.factor(NA),
                       s_surf = factor(NA)))
    tmb_parameters$s_surf <- 0
    tmb_parameters$log_lambda_surf <- 0 
    sm$S_surface <- as(matrix(0, 1, 1), "sparseMatrix")
    sm$S_surface_n <- 0
    sm$A_grid_surface <- matrix(1, ncol(sm$indS), 2)
  } else {
    random <- c(random, "s_surf")
    tmb_parameters$s_surf <- rep(0, ncol(sm$S_surface))
  }
  
  # Discrete random effect
  tmb_parameters$rf_dive <- rep(0, nrow(dat))
  tmb_parameters$rf_surf <- rep(0, nrow(dat))
  if (re ==  "none") {
    map <- c(map, list(log_rf_sd = factor(c(NA, NA)), 
                       rf_dive = factor(rep(NA, nrow(dat))), 
                       rf_surf = factor(rep(NA, nrow(dat)))))
    random2 <- NULL
  } else {
    random <- c(random, "rf_dive", "rf_surf")
    random2 <- c("rf_dive", "rf_surf")
    if (re == "corr") {
      tmb_parameters$log_rf_sd <- rep(0, 3)
    }
  }

  ## Setup TMB data list
  tmb_dat <- list(Xdive = sm$Xs_dive,
                  Xsurf = sm$Xs_surface,
                  S_dive = sm$S_dive,
                  S_dive_n = as.integer(sm$S_dive_n),
                  S_surface = sm$S_surface,
                  S_surface_n = as.integer(sm$S_surface_n),
                  include_smooths = 1, 
                  re = switch(re, none = -1, ind = 1, corr = 2), 
                  A_dive = sm$A_dive,
                  A_surf = sm$A_surf,
                  A_grid_surface = sm$A_grid_surface,
                  A_grid_dive = sm$A_grid_dive,
                  W_dive = sm$W_dive, 
                  W_surf = sm$W_surf, 
                  W_grid_surf = sm$W_grid_surf, 
                  W_grid_dive = sm$W_grid_dive, 
                  Xs_grid_surface = sm$Xs_grid_surface,
                  Xs_grid_dive = sm$Xs_grid_dive,
                  indD = sm$indD,
                  indS = sm$indS,
                  tindD = t(sm$indD),
                  tindS = t(sm$indS), 
                  flag = 1L,
                  dt = sm$dt)

  ## Create x
  if (print) cat("Making AD fun.......")
  obj <- MakeADFun(tmb_dat, tmb_parameters, random = random, map = map,
                   hessian = TRUE, DLL = "ctmc_dive", silent=!print)
  tmb_dat$include_smooths <- -1 
  obj_full <- MakeADFun(tmb_dat, tmb_parameters, map = map,
                        random = random2, 
                        hessian = TRUE, DLL = "ctmc_dive", silent=!print)
  if(print) cat("done\n")

  # find the normalization for the GMRF (smooths)
  if (print) cat("Finding normalization.......")
  obj <- normalize(obj, flag="flag")
  if(print) cat("done\n")

  ## Fit Model
  if (print) cat("Fitting model.......\n")
  t0 <- Sys.time()
  mod <- nlminb(obj$par, obj$fn, gradient = obj$gr)
  t1 <- Sys.time()
  diff <- difftime(t1, t0)
  if(print) cat("Model fit in ", signif(diff[[1]], 2), attr(diff, "units"), "\n")

  ## Compute estimates
  if(print) cat("Estimating variance.......\n")
  npar <- sum(len)
  # want to return the full vcov (well, precision)
  rep <- sdreport(obj, getJointPrecision=TRUE)
  # let's invert to get variance
  if(!is.null(rep$jointPrecision)) {
    vcov <- solve(rep$jointPrecision)
  } else {
    vcov <- NULL
  }

  est <- rep$par.fixed
  var <- rep$cov.fixed
  sds <- sqrt(diag(var))

  ## Compute log-likelihood under full model
  nmsrand <- names(rep$par.random)
  par_all <- c(est, rep$par.random[!(nmsrand %in% c("rf_dive", "rf_surf"))])
  llk_full <- -obj_full$fn(par_all)
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
              vcov = vcov, 
              llk_full = llk_full,
              Xs_dive = sm$Xs_dive,
              Xs_surface = sm$Xs_surface,
              sm = sm,
              forms = forms,
              len = len,
              dat = dat, 
              min_dwell = min_dwell, 
              series = series,
              dt = dt)
  if(print) cat("done\n")
  class(ans) <- "CTMCdive"
  return(ans)
}

#' Print summary of CTMCdive x
#'
#' @param x the CTMCdive fitted model x
#' @param \dots unused (for S3 compatability)
#'
#' @return prints a summary
#' @aliases print.CTMCdive
#' @export
summary.CTMCdive <- function(x, ...) {
  cat("Continuous-time Markov Chain Diving Model\n")
  cat("Number of Dives: ", nrow(x$dat), "\n")
  cat("\n")
  cat("Model formulae:\n")
  cat("  ")
  print(x$forms[[1]])
  cat("  ")
  print(x$forms[[2]])
  cat("\n")

  cat("AIC:", AIC(x)$AIC,"\n")
  cat("\n")

  cat(rep("-", 30), "\n")
  cat("DIVE INTENSITY\n")
  cat("\nFixed effects:\n")
  print(signif(x$res$dive, 4))
  cat("\n")

  if(!is.null(x$sm$s_dive_names)) {
    sdive <- data.frame(Term = x$sm$s_dive_names)
    sdive[["k'"]] <- NA
    sdive$EDF <- NA

    starti <- 1
    lambda <- x$rep$par.fixed
    lambda <- exp(lambda[grepl("^log_lambda_dive", names(lambda))])
    for(i in seq_along(x$sm$S_dive_n)){
      Sn <- x$sm$S_dive_n[i]
      S <- x$sm$S_dive[starti:(starti+Sn-1),
                            starti:(starti+Sn-1), drop=FALSE]
      X <- x$sm$A_dive[,starti:(starti+Sn-1), drop=FALSE]

      sdive[["k'"]][i] <- nrow(S)
      sdive$EDF[i] <- EDF_f(X, S, lambda[i], Sn)
      starti <- starti + Sn
    }

    sdive$EDF <- signif(sdive$EDF, 4)
    cat("\nSmooth terms:\n")
    print(sdive)
    cat("\n")
  }

  cat("\nkappa: ",
      signif(exp(x$rep$par.fixed[names(x$rep$par.fixed) == "log_kappa_dive"]), 4))
  cat("\n")
  cat("\n")

  cat(rep("-", 30), "\n")

  cat("SURFACE INTENSITY\n")
  cat("\nFixed effects:\n")
  print(signif(x$res$surface, 4))

  if(!is.null(x$sm$s_surface_names)) {
    ssurface <- data.frame(Term = x$sm$s_surface_names)
    ssurface[["k'"]] <- NA
    ssurface$EDF <- NA

    starti <- 1
    lambda <- x$rep$par.fixed
    lambda <- exp(lambda[grepl("^log_lambda_surf", names(lambda))])
    for(i in seq_along(x$sm$S_surface_n)){
      Sn <- x$sm$S_surface_n[i]
      S <- x$sm$S_surface[starti:(starti+Sn-1),
                               starti:(starti+Sn-1), drop=FALSE]
      X <- x$sm$A_surf[,starti:(starti+Sn-1), drop=FALSE]

      ssurface[["k'"]][i] <- nrow(S)
      ssurface$EDF[i] <- EDF_f(X, S, lambda[i], Sn)
      starti <- starti + Sn
    }

    ssurface$EDF <- signif(ssurface$EDF, 4)
    cat("\nSmooth terms:\n")
    print(ssurface)
    cat("\n")
  }
  cat("\nkappa: ",
      signif(exp(x$rep$par.fixed[names(x$rep$par.fixed) == "log_kappa_surf"]), 4))
  cat("\n")
  cat("\n")



  cat("\n")
  invisible(x)
}

#' @export
print.CTMCdive <- function(x)
  summary.CTMCdive(x)

#' Predict mean duration from fitted CTMC model for dives and surfacing
#'
#' @param x CTMCdive fitted model x
#' @param newdata  new data frame to make predictions for
#' @param \dots unused (for S3 compatability)
#'
#' @return a list of surface and dive mean predictions
#' @export
predict.CTMCdive <- function(x, newdata = NULL, ...) {
  par <- x$rep$par.fixed
  len <- x$len
  par_dive <- par[1:len[1]]
  par_surf <- par[(len[1] + 1):(len[1] + len[2])]

  if (is.null(newdata)) {
    # linear predictors
    lambda_dive <- x$sm$Xs_grid_dive %*% par_dive
    lambda_surf <- x$sm$Xs_grid_surface %*% par_surf
    if (any(names(x$rep$par.random) == "s_surf")) {
      x_surf <- x$rep$par.random[names(x$rep$par.random) == "s_surf"]
      lambda_surf <- lambda_surf + (x$sm$A_grid_surface %*% x_surf)
    }
    if (any(names(x$rep$par.random) == "s_dive")) {
      x_dive <- x$rep$par.random[names(x$rep$par.random) == "s_dive"]
      lambda_dive <- lambda_dive + (x$sm$A_grid_dive %*% x_dive)
    }
    # weibull
    log_kappa_dive <- x$rep$par.fixed["log_kappa_dive"]
    log_kappa_surf <- x$rep$par.fixed["log_kappa_surf"]
    kappa_dive <- exp(log_kappa_dive)
    kappa_surf <- exp(log_kappa_surf)
    lambda_dive <- kappa_dive * lambda_dive + log_kappa_dive
    lambda_surf <- kappa_surf * lambda_surf + log_kappa_surf
  } else {
    lambda_dive <- newdata$lambda_dive
    lambda_surf <- newdata$lambda_surf 
    kappa_dive <- newdata$kappa_dive
    kappa_surf <- newdata$kappa_surf 
  }

  # correct data for tag 
  if (x$series) {
    dat <- x$dat[x$dat$start == 1,] 
  } else {
    dat <- x$dat
  }

  # get expectations
  ints <- x$sm$ints
  dt <- mean(diff(ints))
  exp_dive <- exp_surf <- rep(0, nrow(dat))
  r_dive <- r_surf <- rep(0, nrow(dat))
  for (i in 1:nrow(dat)) {
    nearest_int_time <- which.min((ints - dat$time[i])^2)
    t <- ints[nearest_int_time]
    Ldive <- lambda_dive[ints > t + dat$dive[i]] + (kappa_dive - 1) * 
      log(ints[ints > t + dat$dive[i]] - t - dat$dive[i])
    Ldive <- cumsum(exp(Ldive))
    Lsurf <- lambda_surf[ints > t] + (kappa_surf - 1) * log(ints[ints > t] - t)
    Lsurf <- cumsum(exp(Lsurf))
    Sdive <- exp(-Ldive * dt)
    Ssurf <- exp(-Lsurf * dt)
    # E(dive) comes from surface intensity and vice-versa
    exp_dive[i] <- sum(Ssurf * dt)
    exp_surf[i] <- sum(Sdive * dt)
    # get dive p-value
    if (sum(ints > t) > 1) {
      Ssurf_approx <- approxfun(c(0, ints[ints > t] - t), c(1, Ssurf))
      r_dive[i] <- qnorm(Ssurf_approx(dat$dive[i]))
    } 
    if (sum(ints > t + dat$dive[i])) {
      Sdive_approx <- approxfun(c(0, ints[ints > t + dat$dive[i]]) - t - dat$dive[i], c(1, Sdive))
      r_surf[i] <- qnorm(Sdive_approx(dat$surface[i]))
    }
  }
  # return predictions
  res <- list(surface = exp_surf, 
              dive = exp_dive, 
              diveI = exp(lambda_dive), 
              surfI = exp(lambda_surf), 
              rdive = r_dive, 
              rsurf = r_surf)
  return(res)
}

#' Plot fitted CTMCdive model
#'
#' Plot models fitted by \code{CTMCdive}.
#'
#' @param x fitted CTMCdive model x
#' @param quant the quantile to plot the observed data up to, for example if set to 0.95 then only the durations up to the 95\% quantile are plotted (prevents outliers dominating the plot)
#' @param pick if not \code{NULL} then either "dive" or "surface" if only want a plot of one of the fitted
#' processes
#' @param pred if predicted means already computed can supply them here rather than have them recomputed
#' @param xlim limits for x axis if you want to restrict range
#' @param se display uncertainty bands around the lines (takes time to compute)
#' @param n_samp number of posterior samples to use if \code{se=TRUE}, otherwise ignored
#' @param \dots unused (for S3 compatability)
#'
#' @return plots of dive and surface durations or only one of if "\code{pick}" is specified
#' @export
#' @importFrom graphics lines par plot
plot.CTMCdive <- function(x, quant = 1, pick = NULL, pred = NULL, xlim = NULL, se=FALSE, n_samp=200, pch=19, cex=0.6, ...) {
  if (is.null(pick)) pick <- "all"
  if (pick == "all") {
    par(mfrow=c(2, 1))
    on.exit(par(mfrow=c(1,1)))
    pick <- "all"
  }
  # get predicted values
  if (is.null(pred)) pred <- predict(x)
  # correct data for tag 
  if (x$series) {
    dat <- x$dat[x$dat$start == 1,] 
  } else {
    dat <- x$dat
  }
  time <- dat$time

  if(se){
    # get uncertainty
    ss <- get_samples(x, n_samp)
    # function to apply
    afn <- function(pars, m, nms.fixed, nms.random){
      lamnm <- grepl("^log_lambda", nms.fixed)
      len <- length(m$rep$par.fixed) - sum(lamnm)
      m$rep$par.fixed[!lamnm] <- pars[1:len]
      m$rep$par.random <- pars[-(1:len)]
      names(m$rep$par.fixed) <- nms.fixed
      names(m$rep$par.random) <- nms.random
      pp <- predict(m)
      return(c(pp$surface, pp$dive))
    }
    nms.fixed <- names(x$rep$par.fixed)
    nms.random <- names(x$rep$par.random)
    samples <- apply(ss, 1, afn, m=x, nms.fixed = nms.fixed, nms.random = nms.random)
    surface_samples <- samples[1:length(time), ]
    dive_samples <- samples[-(1:length(time)), ]
  }

  # plot fitted values over observed
  if (pick == "all" | pick == "surface") {
    q <- quantile(dat$surface, prob = quant)
    plot(time, dat$surface, xlab = "Time", ylab = "Surface duration",  ylim = c(min(dat$surface), q), xlim = xlim, type="n", ...)
    if(se){
      surface_upper <- apply(surface_samples, 1, quantile, probs=0.975)
      surface_lower <- apply(surface_samples, 1, quantile, probs=0.025)
      polygon(c(time, rev(time)), c(surface_lower, rev(surface_upper)), col="grey80", border=NA)
    }
    points(time, dat$surface, col = "grey60", pch=pch, cex=cex)
    lines(time, pred$surface)
  }
  if (pick == "all" | pick == "dive") {
    q <- quantile(dat$dive, prob = quant)
    plot(time, dat$dive, xlab = "Time", ylab = "Dive duration",  ylim = c(min(dat$dive), q), xlim = xlim, type="n", ...)
    if(se){
      dive_upper <- apply(dive_samples, 1, quantile, probs=0.975)
      dive_lower <- apply(dive_samples, 1, quantile, probs=0.025)
      polygon(c(time, rev(time)), c(dive_lower, rev(dive_upper)), col="grey80", border=NA)
    }
    points(time, dat$dive, col = "grey60", pch=pch, cex=cex)
    lines(time, pred$dive)
  }
  invisible(list(mod = x, pred = pred))
}

#' internal function to get posterior samples
#'
#' @param mod a fitted model x
#' @param n number of samples to take
#' @importFrom mgcv rmvn
get_samples <- function(mod, n=200){

  # extract the joint precision
  prec <- mod$rep$jointPrecision
  # remove smoopars
  prec <- prec[!grepl("log_lambda_", colnames(prec)),
               !grepl("log_lambda_", colnames(prec)), drop=FALSE]
  # solve to get variance-covariance matrix
  vc <- solve(prec)

  # get pars
  pars <- c(mod$rep$par.fixed, mod$rep$par.random)
  pars <- pars[!grepl("log_lambda_", names(pars))]

  samps <- rmvn(n, pars, vc)
  colnames(samps) <- names(pars)
  return(samps)
}

#' Returns log-likelihood with degrees of freedom
#'
#' @param x fitted CTMCdive x
#' @param \dots unused (for S3 compatability)
#'
#' @return logLik with df attribute
#' @note This does not account for degrees of freedom reduction with smooths (i.e if lambda > 0)
#' @export
logLik.CTMCdive <- function(x, full = TRUE, ...) {
  if (full) {
    val <- x$llk_full 
  } else {
    val <- -x$mod$objective
  }
  npar <- length(x$mod$par)
  llk <- val
  attributes(llk)$df <- npar
  return(llk)
}

#' Akaike's An Information Criterion for CTMCdive models
#'
#' Calculate the AIC from a fitted model.
#'
#' @param x a fitted detection function x
#' @param k penalty per parameter to be used; the default \code{k = 2} is the "classical" AIC
#' @param \dots optionally more fitted model xs.
#' @author David L Miller, Richard Glennie 
#' @export
#' @importFrom stats logLik
AIC.CTMCdive <- function(x, ..., k=2){

  # get the models
  models <- list(x, ...)

  # build the table
  aics <- matrix(NA, nrow=length(models), ncol=2)
  for(i in seq_along(models)){
    this_mod <- models[[i]]
    ll <- this_mod$llk_full
    
    # get joint covariance 
    vcov <- this_mod$vcov
    vcov_nms <- colnames(this_mod$rep$jointPrecision)
    
    # dive component
    vdive <- vcov[vcov_nms == "s_dive", vcov_nms == "s_dive"]
    if (!is.null(vdive)) {
      Xdive <- this_mod$sm$A_dive
      Idive <- t(Xdive) %*% Xdive 
      dive_edf <- sum(diag(vdive %*% Idive))
    } else {
      dive_edf <- 0 
    }
    
    # surface component
    vsurf <- vcov[vcov_nms == "s_surf", vcov_nms == "s_surf"]
    if (!is.null(vsurf)) {
      Xsurf <- this_mod$sm$A_surf
      Isurf <- t(Xsurf) %*% Xsurf 
      surface_edf <- sum(diag(vsurf %*% Isurf))
    } else {
      surface_edf <- 0
    }
    
    # DF for fixed effects (not smoopars)
    parnms <- names(this_mod$rep$par.fixed)
    fixed_edf <- sum(!(parnms %in% c("log_lambda_dive","log_lambda_surf")))
    # total
    total_edf <- dive_edf + surface_edf + fixed_edf

    aics[i, 1] <- total_edf
    aics[i, 2] <- -2*ll + k*aics[i, 1]
  }
  # make it a data.frame
  aics <- as.data.frame(aics)
  names(aics) <- c("df", "AIC")
  # add row names
  call <- match.call(expand.dots=TRUE)
  rownames(aics) <- as.character(call)[-1]
  # sort by AIC
  aics <- aics[order(aics[,2]),]
  return(aics)

}

#' EDF for smooth terms = trace(F)
#'
#' @param X design matrix 
#' @param S smoothing matrix 
#' @param lambda smoothing parameter 
#' @param Sn dimension of smoothing matrix 
#'
#' @return trace of F = (Xt X + sp*S)^-1 Xt X
#' @importFrom Matrix t solve diag 
EDF_f <- function(X, S, lambda, Sn){

  if(length(lambda)==0){
    return(0)
  }
  # duplicate lambda enough times
  lambda <- rep(lambda, Sn)
  # calculate lambda*S
  Sbig <- S * lambda

  # calculate the hat matrix
  XtX <- t(X) %*% X
  Fi <- solve(XtX + Sbig)
  F <- Fi %*% XtX

  # return the trace
  sum(diag(F))
}


#' Estimate exposure effect for dive and surface 
#'
#' @param mod fitted CTMC model 
#' @param predgrid data frame with same columns as data but with a user-defined time grid to estimate the effect over 
#'                 must have a column which is 1 when exposure is one and zero otherwise
#' @param exp name of exposure indicator 
#' @parama val column name of variable used as magnitude for exposure (e.g. seconds since) if time not used as smooth
#' @param nsims number of posterior simulations, default is 1000
#'
#' @return list with dive and surface elements which each contain the mean exposure effect and the credible interval bands
#' @importFrom Matrix solve t rowMeans
#' @export
GetExposureEff <- function(mod, predgrid, basedat, nsims = 1000) {
  rand <- mod$rep$par.random
  fixed <- mod$rep$par.fixed
  Xdive <- predict(mod$sm$gam_dive, predgrid, type = "lpmatrix")
  Xdivebase <- predict(mod$sm$gam_dive, basedat, type = "lpmatrix")
  Xsurf <- predict(mod$sm$gam_surf, predgrid, type = "lpmatrix")
  Xbasesurf <- predict(mod$sm$gam_surf, basedat, type = "lpmatrix")
  betas <- get_samples(mod, n = nsims)
  betadive <- t(betas[, c(colnames(betas) %in% c("par_dive", "s_dive"))])
  betasurf <- t(betas[, c(colnames(betas) %in% c("par_surf", "s_surf"))])
  simdive <- exp(Xdive %*% betadive)
  simbasedive <- exp(Xdivebase %*% betadive)
  simsurf <- exp(Xsurf %*% betasurf)
  simbasesurf <- exp(Xbasesurf %*% betasurf)
  simdiveexp <- (simdive - simbasedive) 
  simsurfexp <- (simsurf - simbasesurf)  
  dive <- surf <- NULL
  dive$mean <- rowMeans(simdiveexp)
  dive$ci <- apply(simdiveexp, 1, quantile, prob = c(0.025, 0.975))
  surf$mean <- rowMeans(simsurfexp)
  surf$ci <- apply(simsurfexp, 1, quantile, prob = c(0.025, 0.975))
  dive$sims <- simdiveexp
  surf$sims <- simsurfexp
  res <- list(dive = dive, surf = surf, predgrid = predgrid)
  class(res) <- "ExposureEffect"
  return(res) 
}

#' Expand covariates across a time grid 
#'
#' @param dat data frame with covariates in columns and time column
#' @param tgrid a vector with new time grid
#'
#' @return a expanded data frame where covariates are copied from nearest time point to grid point 
#' @export
ExpandCovs <- function(dat, tgrid) {
  subdat <- dat[dat$time > min(tgrid) & dat$time <= max(tgrid),]
  freq <- as.numeric(table(sapply(tgrid, FUN = function(x) which.min(abs(x - subdat$time)))))
  predgrid <- subdat[rep(seq(nrow(subdat)), freq), ] 
  predgrid$time <- tgrid
  return(predgrid)
}

#' Plot exposure effect 
#'
#' @param expeff exposure effect estimates from GetExposureEff function
#' @param pick if "both" (default) then for dive and surface, otherwise for "dive" or "surface" 
#'
#' @return plots of exposure effect with vertical line at point where no evidence exposure has effect from baseline
#' @export
plot.ExposureEffect <- function(expeff, pick = "all") {
  if (is.null(pick)) pick <- "all"
  if (pick == "all") {
    par(mfrow=c(2, 1))
    on.exit(par(mfrow=c(1,1)))
    pick <- "all"
  }
  if (pick == "dive" | pick == "all") {
    sig <- ifelse(expeff$dive$mean > 0, expeff$dive$ci[1,] > 0, expeff$dive$ci[2,] < 0)
    cols <- ifelse(sig, "blue", "grey80")
    plot(expeff$predgrid$time, expeff$dive$mean, type = "n", ylim = range(expeff$dive$ci), xlab = "Time", ylab = "Diving intensity exposure effect")
    lines(expeff$predgrid$time, expeff$dive$ci[1,], lty = "dashed", lwd = 1.5)
    lines(expeff$predgrid$time, expeff$dive$ci[2,], lty = "dashed", lwd = 1.5)
    n <- length(expeff$predgrid$time)
    segments(expeff$predgrid$time[-n], expeff$dive$mean[-n], 
             expeff$predgrid$time[-1], expeff$dive$mean[-1], 
             col = cols)
    abline(h = 0, col = "red", lty = "dotted")
  } 
  if (pick == "surf" | pick == "all") {
    sig <- ifelse(expeff$surf$mean > 0, expeff$surf$ci[1,] > 0, expeff$surf$ci[2,] < 0)
    cols <- ifelse(sig, "blue", "grey80")
    plot(expeff$predgrid$time, expeff$surf$mean, type = "n", ylim = range(expeff$surf$ci), xlab = "Time", ylab = "Surfacing intensity exposure effect")
    lines(expeff$predgrid$time, expeff$surf$ci[1,], lty = "dashed", lwd = 1.5)
    lines(expeff$predgrid$time, expeff$surf$ci[2,], lty = "dashed", lwd = 1.5)
    n <- length(expeff$predgrid$time)
    segments(expeff$predgrid$time[-n], expeff$surf$mean[-n], 
             expeff$predgrid$time[-1], expeff$surf$mean[-1], 
             col = cols, lwd = 1.5)
    abline(h = 0, col = "red", lty = "dotted", lwd = 1.5)
  }
}

#' Add or remove terms from a CTMC model and fit the new model 
#'
#' @param mod model to be amended 
#' @param change change to formula to make, see ?update.formula 
#' @param which which formula to update, 1 = dive, 2 = surface; if which = 0
#'   then three models are fit where change is made to each formula individually
#'   and then to both, an AIC table is printed. 
#'
#' @return fitted model if which = 1,2 or list of three fitted models if which = 0
#' @export
update.CTMCdive <- function(mod, change, which = 0) {
  if (which == 0) {
    ms <- vector(mode = "list", length = 3)
    f <- mod$forms  
    f1 <- f 
    f1[["dive"]] <- update(f[["dive"]], change)
    dive <- ms[[1]] <- try(FitCTMCdive(f1, mod$dat, min_dwell = mod$min_dwell, series = mod$series, dt = mod$dt))
    f1 <- f
    f1[["surface"]] <- update(f[["surface"]], change)
    surf <- ms[[2]] <- try(FitCTMCdive(f1, mod$dat, min_dwell = mod$min_dwell, series = mod$series, dt = mod$dt))
    f1 <- lapply(f, FUN = function(fi) {update(fi, change)})
    both <- ms[[3]] <- try(FitCTMCdive(f1, mod$dat, min_dwell = mod$min_dwell, series = mod$series, dt = mod$dt))
    aics <- try(AIC(mod, dive, surf, both))
    names(ms) <- c("dive", "surf", "both")
    print(aics)
    return(ms)
  } else {
    f <- mod$forms
    f[[which]] <- update(f[[which]], change)
    m <- FitCTMCdive(f, mod$dat, min_dwell = mod$min_dwell, series = mod$series, dt = mod$dt)
    return(m)
  }
}

#' Formala translated to a character 
#' Borrow from formula.tools session 
#'  Christopher Brown (2018). formula.tools: Programmatic Utilities for Manipulating Formulas,
#'  Expressions, Calls, Assignments and Other R Objects. R package version 1.7.1.
#'  https://CRAN.R-project.org/package=formula.tools
#'
#' @param x 
#' @param ... 
#'
#' @return character 
#' @export
for2char <- function (x, ...) {
  form <- paste(deparse(x), collapse = " ")
  form <- gsub("\\s+", " ", form, perl = FALSE)
  return(form)
}



