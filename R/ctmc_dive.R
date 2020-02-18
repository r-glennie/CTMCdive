# Fit a CTMC to dive data

#' Compute matrices required to fit the model
#'
#' @param forms formulae
#' @param dat data frame (see fitCTMCdive)
#' @param min_dwell start point for the integrals
#' @param nint number of integration points
#'
#' @return list of matrices required to fit the model
#' @export
#' @importFrom mgcv gam predict.gam
#' @importFrom methods as
#' @importFrom Matrix bdiag
MakeMatrices <- function(forms, dat, min_dwell, nint = 10000) {

  # results list
  res <- list()

  ## dive model
  # GAM setup
  gam_dive <- gam(forms[["dive"]], data = dat, method = "REML")
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

  ## build design matrix
  mm <- predict(gam_dive, newdata = dat, type = "lpmatrix")
  # fixed effects design matrix
  res$Xs_dive <- mm[, 1:gam_dive$nsdf, drop=FALSE]
  # smooth design matrix
  res$A_dive <- mm[, -c(1:gam_dive$nsdf), drop=FALSE]
  # weibull design matrix 
  res$W_dive <- matrix(c(0, log(dat$surface[-1] + 1e-10)), nr = nrow(dat), nc = 1)
  
  ## surface model
  # GAM setup
  gam_surface <- gam(forms[["surface"]], data = dat, method = "REML")
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
  surfdat <- dat
  surfdat$time <- dat$time + dat$dive

  ## build design matrix
  mm <- predict(gam_surface, newdata = surfdat, type = "lpmatrix")
  # fixed effects design matrix
  res$Xs_surface <- mm[, 1:gam_surface$nsdf, drop=FALSE]
  # smooth design matrix
  res$A_surf <- mm[, -c(1:gam_surface$nsdf), drop=FALSE]
  # weibull design matrix
  res$W_surf <- matrix(log(dat$dive + 1e-10), nr = nrow(dat), nc = 1)
  
  # construct prediction grid
  n <- nrow(dat)
  ints <- seq(dat$time[1],
              dat$time[n] + dat$dive[n] + dat$surface[n],
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
  t_dive <- c(-1e-10, dat$time + dat$dive, max(dat$time + dat$surface + dat$dive))
  cut_dive <- as.numeric(cut(ints, breaks = t_dive))
  dt_dive <- ints - t_dive[cut_dive]
  res$W_grid_dive <- log(dt_dive + 1e-10)
  
  pred_mat_surface <- pred_mat_maker(gam_surface, dat, ints)
  res$Xs_grid_surface <- pred_mat_surface[, 1:gam_surface$nsdf, drop=FALSE]
  res$A_grid_surface <- pred_mat_surface[, -c(1:gam_surface$nsdf), drop=FALSE]

  # weibull for integration grid 
  t_surf <- c(-1e-10, dat$time, max(dat$time + dat$surface + dat$dive))
  cut_surf <- as.numeric(cut(ints, breaks = t_surf))
  dt_surf <- ints - t_surf[cut_surf] 
  res$W_grid_surf <- log(dt_surf + 1e-10)
  
  # get integration matrices
  indD <- indS <- matrix(0, nrow = n, ncol = nint)
  dive_end <- dat$time + dat$dive
  for (i in 1:nint) {
    wh <- which(ints[i] <= dat$time[-1] &
                ints[i] >= (dive_end[-n]+min_dwell$surface))
    if (length(wh) != 0) {
      indD[wh, i] <- 1
    }
    wh <- which(ints[i] >= (dat$time+min_dwell$dive) & ints[i] <= dive_end)
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
#'
#' @return a CTMCdive model object: a list of the estimated results (\code{res}), variance matrix (\code{var}), fitted model returned from \code{optim} (\code{mod}), output from \code{sdreport} (\code{res}), design matrices (\code{Xs}), smoothing data (\code{sm}), formulae (\code{forms{), indices that divide par between dive and surface parameters (\code{len}), data (\code{dat}), indicators for smooths used (\code{lambda}), and model type (\code{model}).
#' @export
#' @importFrom stats delete.response model.matrix optim
#'             pnorm predict qnorm quantile terms
#' @importFrom TMB MakeADFun sdreport
#' @useDynLib ctmc_dive
FitCTMCdive <- function(forms, dat, print = TRUE,
                        min_dwell=list(dive=0, surface=0)) {

  ## Make design matrices and parameter vector
  if (print) cat("Computing design matrices.......")

  # name that list
  names(forms) <- c(as.character(forms[[1]][[2]]),
                    as.character(forms[[2]][[2]]))

  # smoothing data
  random <- NULL
  map <- list()
  sm <- MakeMatrices(forms, dat, min_dwell)

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
                         log_lambda_surf = rep(0, length(sm$S_surface_n)))

  # if there are smooths of dive or surface, set up the TMB
  # model objects correctly
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

  ## Setup TMB data list
  tmb_dat <- list(Xdive = sm$Xs_dive,
                  Xsurf = sm$Xs_surface,
                  S_dive = sm$S_dive,
                  S_dive_n = as.integer(sm$S_dive_n),
                  S_surface = sm$S_surface,
                  S_surface_n = as.integer(sm$S_surface_n),
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
                  flag = 1L,
                  dt = sm$dt)

  ## Create object
  if (print) cat("Making AD fun.......")
  obj <- MakeADFun(tmb_dat, tmb_parameters, random = random, map = map,
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
#' @aliases print.CTMCdive
#' @export
summary.CTMCdive <- function(object, ...) {
  cat("Continuous-time Markov Chain Diving Model\n")
  cat("Number of Dives: ", nrow(object$dat), "\n")
  cat("\n")
  cat("Model formulae:\n")
  cat("  ")
  print(object$forms[[1]])
  cat("  ")
  print(object$forms[[2]])
  cat("\n")

  cat("AIC:", AIC(object)$AIC,"\n")
  cat("\n")

  cat(rep("-", 30), "\n")
  cat("DIVE INTENSITY\n")
  cat("\nFixed effects:\n")
  print(signif(object$res$dive, 4))
  cat("\n")

  if(!is.null(object$sm$S_dive)) {
    sdive <- data.frame(Term = object$sm$s_dive_names)
    sdive[["k'"]] <- NA
    sdive$EDF <- NA

    starti <- 1
    for(i in seq_along(object$sm$S_dive_n)){
      lambda <- object$rep$par.fixed
      lambda <- exp(lambda[grepl("^log_lambda_dive", names(lambda))])
      Sn <- object$sm$S_dive_n[i]
      S <- object$sm$S_dive[starti:(starti+Sn-1),
                            starti:(starti+Sn-1), drop=FALSE]
      X <- object$sm$A_dive[starti:(starti+Sn-1),
                            starti:(starti+Sn-1), drop=FALSE]

      sdive[["k'"]][i] <- nrow(S)
      sdive$EDF[i] <- EDF_f(X, S, lambda, Sn)
      starti <- starti + Sn
    }

    sdive$EDF <- signif(sdive$EDF, 4)
    cat("\nSmooth terms:\n")
    print(sdive)
    cat("\n")
  }

  cat(rep("-", 30), "\n")

  cat("SURFACE INTENSITY\n")
  cat("\nFixed effects:\n")
  print(signif(object$res$surface, 4))

  if(!is.null(object$sm$S_surface)) {
    ssurface <- data.frame(Term = object$sm$s_surface_names)
    ssurface[["k'"]] <- NA
    ssurface$EDF <- NA

    starti <- 1
    for(i in seq_along(object$sm$S_surface_n)){
      lambda <- object$rep$par.fixed
      lambda <- exp(lambda[grepl("^log_lambda_surf", names(lambda))])
      Sn <- object$sm$S_surface_n[i]
      S <- object$sm$S_surface[starti:(starti+Sn-1),
                               starti:(starti+Sn-1), drop=FALSE]
      X <- object$sm$A_surf[starti:(starti+Sn-1),
                            starti:(starti+Sn-1), drop=FALSE]

      ssurface[["k'"]][i] <- nrow(S)
      ssurface$EDF[i] <- EDF_f(X, S, lambda, Sn)
      starti <- starti + Sn
    }

    ssurface$EDF <- signif(ssurface$EDF, 4)
    cat("\nSmooth terms:\n")
    print(ssurface)
    cat("\n")
  }



  cat("\n")
  invisible(object)
}

#' @export
print.CTMCdive <- function(x)
  summary.CTMCdive(x)

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

  if(!is.null(newdata)){
    stop("New data not supported yet!")
  }

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

  if (any(names(object$rep$par.random) == "s_surf")) {
    x_surf <- object$rep$par.random[names(object$rep$par.random) == "s_surf"]
    lambda_surf <- lambda_surf + (object$sm$A_grid_surface %*% x_surf)
  }
  if (any(names(object$rep$par.random) == "s_dive")) {
    x_dive <- object$rep$par.random[names(object$rep$par.random) == "s_dive"]
    lambda_dive <- lambda_dive + (object$sm$A_grid_dive %*% x_dive)
  }
  # weibull
  log_kappa_dive <- mod$rep$par.fixed["log_kappa_dive"]
  log_kappa_surf <- mod$rep$par.fixed["log_kappa_surf"]
  kappa_dive <- exp(log_kappa_dive)
  kappa_surf <- exp(log_kappa_surf)
  lambda_dive <- kappa_dive * lambda_dive + log_kappa_dive
  lambda_surf <- kappa_surf * lambda_surf + log_kappa_surf

  # get expectations
  dt <- mean(diff(ints))
  exp_dive <- exp_surf <- rep(0, nrow(object$dat))
  r_dive <- r_surf <- rep(0, nrow(object$dat))
  for (i in 1:nrow(object$dat)) {
    Ldive <- lambda_dive[ints > object$dat$time[i] + object$dat$dive[i]] + (kappa_dive - 1) * 
      log(ints[ints > object$dat$time[i] + object$dat$dive[i]] - object$dat$time[i] - object$dat$dive[i])
    Ldive <- cumsum(exp(Ldive))
    Lsurf <- lambda_surf[ints > object$dat$time[i]] + (kappa_surf - 1) * 
      log(ints[ints > object$dat$time[i]] - object$dat$time[i])
    Lsurf <- cumsum(exp(Lsurf))
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
  if (is.null(pick)) {
    par(mfrow=c(2, 1))
    on.exit(par(mfrow=c(1,1)))
    pick <- "all"
  }
  # get predicted values
  if (is.null(pred)) pred <- predict(x)
  time <- x$dat$time

  if(se){
    # get uncertainty
    ss <- get_samples(mod, n_samp)
    # function to apply
    afn <- function(pars, m){
      m$rep$par.fixed <- pars[1:length(m$rep$par.fixed)]
      m$rep$par.random <- pars[-(1:length(m$rep$par.fixed))]
      pp <- predict(m)
      return(c(pp$surface, pp$dive))
    }
    samples <- apply(ss, 1, afn, m=mod)
    surface_samples <- samples[1:length(time), ]
    dive_samples <- samples[-(1:length(time)), ]
  }

  # plot fitted values over observed
  if (pick == "all" | pick == "surface") {
    q <- quantile(x$dat$surface, prob = quant)
    plot(time, x$dat$surface, xlab = "Time", ylab = "Surface duration",  ylim = c(min(x$dat$surface), q), xlim = xlim, type="n", ...)
    if(se){
      surface_upper <- apply(surface_samples, 1, quantile, probs=0.975)
      surface_lower <- apply(surface_samples, 1, quantile, probs=0.025)
      polygon(c(time, rev(time)), c(surface_upper, surface_lower), col="grey80", border=NA)
    }
    points(time, x$dat$surface, col = "grey60", pch=pch, cex=cex)
    lines(time, pred$surface)
  }
  if (pick == "all" | pick == "dive") {
    q <- quantile(x$dat$dive, prob = quant)
    plot(time, x$dat$dive, xlab = "Time", ylab = "Dive duration",  ylim = c(min(x$dat$dive), q), xlim = xlim, type="n", ...)
    if(se){
      dive_upper <- apply(dive_samples, 1, quantile, probs=0.975)
      dive_lower <- apply(dive_samples, 1, quantile, probs=0.025)
      polygon(c(time, rev(time)), c(dive_upper, dive_lower), col="grey80", border=NA)
    }
    points(time, x$dat$dive, col = "grey60", pch=pch, cex=cex)
    lines(time, pred$dive)
  }
  invisible(list(mod = x, pred = pred))
}

#' internal function to get posterior samples
#'
#' @param mod a fitted model object
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

  rmvn(n, pars, vc)
}

#' Returns log-likelihood with degrees of freedom
#'
#' @param object fitted CTMCdive object
#' @param \dots unused (for S3 compatability)
#'
#' @return logLik with df attribute
#' @note This does not account for degrees of freedom reduction with smooths (i.e if lambda > 0)
#' @export
logLik.CTMCdive <- function(object, ...) {
  val <- -object$mod$objective
  npar <- length(object$mod$par)
  llk <- val
  attributes(llk)$df <- npar
  return(llk)
}

#' Akaike's An Information Criterion for CTMCdive models
#'
#' Calculate the AIC from a fitted model.
#'
#' @param object a fitted detection function object
#' @param k penalty per parameter to be used; the default \code{k = 2} is the "classical" AIC
#' @param \dots optionally more fitted model objects.
#' @author David L Miller
#' @export
#' @importFrom stats logLik
AIC.CTMCdive <- function(object, ..., k=2){

  # get the models
  models <- list(object, ...)

  # build the table
  aics <- matrix(NA, nrow=length(models), ncol=2)
  for(i in seq_along(models)){
    this_mod <- models[[i]]
    ll <- logLik(this_mod)

    # dive component
    dive_lambda <- exp(this_mod$rep$par.fixed[grepl("^log_lambda_dive", names(this_mod$rep$par.fixed))])
    dive_edf <- EDF_f(this_mod$sm$A_dive, this_mod$sm$S_dive, dive_lambda, this_mod$sm$S_dive_n)
    # surface component
    surface_lambda <- exp(this_mod$rep$par.fixed[grepl("^log_lambda_surf", names(this_mod$rep$par.fixed))])
    surface_edf <- EDF_f(this_mod$sm$A_surf, this_mod$sm$S_surface, surface_lambda, this_mod$sm$S_surface_n)
    # DF for fixed effects
    fixed_edf <- length(this_mod$rep$par.fixed) - (length(dive_lambda) + length(surface_lambda))
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

  return(aics)

}


# EDF for smooth terms = trace(F)
# F = (Xt X + sp*S)^-1 Xt X
# Wood 2017 p 212
EDF_f <- function(X, S, lambda, Sn){
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
