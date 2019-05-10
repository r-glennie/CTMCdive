# Fit a CTMC to dive data

#' Compute matrices required to fit temporal smooth
#'
#' @param dat data frame (see fitCTMCdive)
#' @param nint number of integration points / knots 
#'
#' @return list of mesh, SPDE objects, point A matrices, and integral A matrices 
#' @export
MakeSmooth <- function(dat) {
  # make mesh over time
  nk <- nrow(dat)
  knots <- seq(min(dat$dive_start), max(dat$dive_end), length = nk)
  mesh <- inla.mesh.1d(knots, degree = 1)
  # get projector matrices
  A_dive <- inla.spde.make.A(mesh, loc = dat$dive_start)
  A_surf <- inla.spde.make.A(mesh, loc = dat$dive_end)
  # make SPDEs 
  spde <- inla.spde2.matern(mesh, alpha = 2)
  res <- list(mesh = mesh, 
              A_dive = A_dive, 
              A_surf = A_surf, 
              C  = spde$param.inla["M0"][[1]], 
              G1 = spde$param.inla["M1"][[1]], 
              G2 = spde$param.inla["M2"][[1]])
  return(res)
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
FitCTMCdive <- function(forms, dat, model = "iid", correlated = FALSE, model.args = NULL, print = TRUE) {
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
  par[1] <- log(mean(dat$dive))
  par[len[1] + 1] <- log(mean(dat$surface))
  if(print) cat("done\n")

  # starting par given?
  if (!is.null(model.args$ini_par)) par <- model.args$ini_par

  ## If model = "cor" then compute smooth information
  random <- NULL
  map <- list()
  sm <- MakeSmooth(dat)
  n_mesh <- sm$mesh$n
  lambda <- rep(0, 4)
  if (model != "iid") {
    lambda <- rep(0, 4)
    if (model == "d") lambda[1:2] <- 1
    if (model == "s") lambda[3:4] <- 1 
    if (model == "ds" | model == "sd") lambda[1:4] <- 1 
  } 
  
  ## Setup TMB parameter list 
  tmb_parameters <- list(par_dive = par[1:len[1]], 
                         par_surf = par[(len[1] + 1): length(par)], 
                         log_sigma_surf = log(sd(dat$dive)), 
                         log_sigma_dive = log(sd(dat$surface)), 
                         cor_surfdive = 0, 
                         log_kappa_dive = 0, 
                         log_kappa_surf = 0,
                         log_tau_dive = 0, 
                         log_tau_surf = 0, 
                         s_dive = rep(0, n_mesh), 
                         s_surf = rep(0, n_mesh)) 
  
  ## Create TMB model object
  if (!(lambda[1] < 1e-10)) {
    random <- c(random, "s_dive")
    # fix intercept 
    par_int <- seq(length(tmb_parameters$par_dive))
    par_int[1] <- NA
    par_int <- as.factor(par_int)
    map <- c(map, list(par_dive = par_int))
    Xs[[1]][,1] <- 0
    tmb_parameters$s_dive <- rep(tmb_parameters$par_dive[1], n_mesh)
  } else {
    map <- c(map,  list(log_tau_dive = as.factor(NA), 
                        log_kappa_dive = as.factor(NA), 
                        s_dive = factor(rep(NA, n_mesh))))
  }
  if (!(lambda[3] < 1e-10)) {
    random <- c(random, "s_surf")
    # remove intercept 
    par_int <- seq(length(tmb_parameters$par_surf))
    par_int[1] <- NA
    par_int <- as.factor(par_int)
    map <- c(map, list(par_surf = par_int))
    Xs[[2]][,1] <- 0 
    tmb_parameters$s_surf <- rep(tmb_parameters$par_surf[1], n_mesh)
  } else {
    map <- c(map,  list(log_tau_surf = as.factor(NA), 
                        log_kappa_surf = as.factor(NA), 
                        s_surf = factor(rep(NA, n_mesh))))
  }
  if (!correlated) map <- c(map, list(cor_surfdive = as.factor(NA)))
  
  ## Setup TMB data list 
  tmb_dat <- list(ldive = log(dat$dive), 
                  lsurf = log(dat$surface), 
                  A_dive = sm$A_dive, 
                  A_surf = sm$A_surf, 
                  Xdive = Xs[[1]],
                  Xsurf = Xs[[2]],
                  C = sm$C, 
                  G1 = sm$G1, 
                  G2 = sm$G2)
  
  ## Create object
  obj <- MakeADFun(tmb_dat, tmb_parameters, random = random, map = map, hessian = TRUE, DLL = "ctmc_dive")
                         
  ## Fit Model 
  if (print) cat("Fitting model.......\n")
  t0 <- Sys.time()
  mod <- nlminb(obj$par, obj$fn, gradient = obj$gr)
  t1 <- Sys.time()
  diff <- difftime(t1, t0)
  if(print) cat("Model fit in ", signif(diff[[1]], 2), attr(diff, "units"), "\n")
  # 
  ## Compute estimates
  if(print) cat("Estimating variance.......")
  npar <- sum(len)
  rep <- sdreport(obj)
  est <- rep$par.fixed
  var <- rep$cov.fixed
  sds <- sqrt(diag(var))
  if (model != "iid") {
    if (lambda[1] > 1e-10) {
      est <- c(rep$value[1], est) 
      sds <- c(rep$sd[1], sds)
    }
    if (lambda[3] > 1e-10) {
      est <- c(est[1:len[1]], rep$value[2], est[(len[1]+1):length(est)])
      sds <- c(sds[1:len[1]], rep$sd[2], sds[(len[1]+1):length(sds)])
    }
  }
  if(print) cat("done\n")
  # 
  # confidence intervals
  if(print) cat("Computing confidence intervals.......")
  LCL <- est - qnorm(0.975) * sds
  UCL <- est + qnorm(0.975) * sds
  if(print) cat("done\n")
  # 
  # pvalues
  pval <- 2 * pnorm(-abs(est), 0, sds)
  # 
  #  dive result
  if(print) cat("Formatting results.......")
  s <- 1
  e <- len[1]
  res_dive <- data.frame(Est. = est[s:e],
                           SD. = sds[s:e],
                           LCL. = LCL[s:e],
                           UCL. = UCL[s:e],
                           p_value = pval[s:e])
  res_dive <- signif(res_dive, 4)
  rownames(res_dive) <- colnames(Xs[[1]])

  # surface result
  s <- len[1] + 1
  e <- len[1] + len[2]
  res_surf <- data.frame(Est. = est[s:e],
                         SD. = sds[s:e],
                         LCL. = LCL[s:e],
                         UCL. = UCL[s:e],
                         p_value = pval[s:e])
  res_surf <- signif(res_surf, 4)
  rownames(res_surf) <- colnames(Xs[[2]])
  
  # variance 
  s <- len[1] + len[2] + 1
  e <- length(est)
  res_cor <- data.frame(Est. = est[s:e],
                         SD. = sds[s:e],
                         LCL. = LCL[s:e],
                         UCL. = UCL[s:e],
                         p_value = pval[s:e])
  res_cor <- signif(res_cor, 4)

  # full result
  res <- list(surface = res_surf, dive = res_dive, cor = res_cor)

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
  if (mod$model != "iid") {
    cat("with 1D Matern smooth for ")
    if (mod$lambda[1] > 1e-10) cat("dive duration")
    if (mod$lambda[1] > 1e-10 & mod$lambda[3] > 1e-10) cat(" and ")
    if (mod$lambda[3] > 1e-10) cat("surfacing duration")
    cat(".")
  }
  cat("\n")
  cat(rep("-", 30), "\n")
  cat("DIVE DURATION\n")
  print(mod$res$dive)
  cat("\n")
  cat(rep("-", 30), "\n")
  cat("SURFACE DURATION\n")
  print(mod$res$surface)
  cat("\n")
  cat(rep("-", 30), "\n")
  cat("VARIANCE MODEL\n")
  print(mod$res$cor)
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
  par <- mod$rep$par.fixed
  len <- mod$len
  if (mod$lambda[1] > 1e-10) par <- c(0, par)
  if (mod$lambda[3] > 1e-10) par <- c(par[1:len[1]], 0, par[(len[1]+1):length(par)])
  par_dive <- par[1:len[1]]
  par_surf <- par[(len[1] + 1):(len[1] + len[2])]
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
  nu_dive <- Xs[[1]] %*% par_dive
  nu_surf <- Xs[[2]] %*% par_surf
  # add RW smooths if necessary
  if (mod$model != "iid") {
    npar <- sum(len)
    lambda <- mod$lambda
    if (lambda[3]>1e-10) {
      x_surf <- mod$rep$par.random[names(mod$rep$par.random) == "s_surf"]
      nu_surf <- nu_surf + (mod$sm$A_surf %*% x_surf)
    }
    if (lambda[1]>1e-10) {
      x_dive <- mod$rep$par.random[names(mod$rep$par.random) == "s_dive"]
      nu_dive <- nu_dive + (mod$sm$A_dive %*% x_dive)
    }
  }
  # predicted durations for diving and surfacing (1/intensity)
  mu_dive <- exp(nu_dive)
  mu_surf <- exp(nu_surf)

  # return predictions
  res <- list(surface = mu_surf, dive = mu_dive)
  return(res)
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
  val <- -object$mod$objective
  npar <- length(object$mod$par)
  llk <- val
  attributes(llk)$df <- npar
  return(llk)
}
