
#' Simulate CTMCdive data
#'
#' @param dive_fn function that given a vector of times computes dive intensity at those times
#' @param surf_fn function that given a vector of times computes surface intensity at those times
#' @param dt time step for simulation 
#'
#' @return CTMC dataset 
#' @export
simulateCTMC <- function(dive_fn, surf_fn, T, dt, sd = NULL, print = TRUE) {
  tgr <- seq(0, T, by = dt)
  # predict intensity on fine grid 
  diveI <- dive_fn(tgr)
  surfI <- surf_fn(tgr)
  nt <- length(tgr)
  s <- rep(0, nt)
  s[1] <- 1 # always start in dive
  diving <- TRUE
  # create empty dataset
  dat <- data.frame(ID = 1, dive = 0, surface = 0, time = 1)
  # record number 
  cur <- 1 
  if (is.null(sd)) {
    rf <- c(0, 0)
  } else {
    rf <- rnorm(2, 0, sd)
  }
  # fill out data
  for (i in 2:nt) {
    if(print) cat(round(100 * (i / nt), 0), "% done \r", sep = "")
    isurf <- exp(log(surfI[i]) + rf[2])
    idive <- exp(log(diveI[i]) + rf[1])
    trm <- matrix(c(-isurf, idive, isurf, -idive), nc = 2)
    tpm <- expm(trm * dt)
    s[i] <- sample(1:2, size = 1, prob = tpm[s[i - 1],])
    if (s[i] == 1) {
      dat$dive[cur] <- dat$dive[cur] + dt
    }
    if (s[i] == 2) {
      dat$surface[cur] <- dat$surface[cur] + dt
      diving <- FALSE
    }
    if (s[i] == 1 & !diving) {
      dat <- rbind(dat, data.frame(ID = 1, dive = 0, surface = 0, time = i*dt))
      cur <- cur + 1
      diving <- TRUE
      if(!is.null(sd)) rf <- rnorm(2, 0, sd)
    }
  }
  dat <- dat[-1,]
  return(dat)
}


#' Simulate CTMCdive data
#'  
#' @param n number of dive-surfacings to simulate 
#' @param dive_fn function that given a vector of times computes dive intensity at those times
#' @param surf_fn function that given a vector of times computes surface intensity at those times
#' @param dt time step for simulation 
#'
#' @return CTMC dataset 
#' @export
simulateCTMC2 <- function(dive_fn, surf_fn, T, dt, sd = NULL, print = TRUE) {
  tgr <- seq(0, T, by = dt)
  # predict intensity on fine grid 
  odiveI <- dive_fn(tgr)
  osurfI <- surf_fn(tgr)
  # current time 
  t <- 0 
  # create empty dataset
  dat <- data.frame(ID = 1, dive = 0, surface = 0, time = 1)
  cur <- 1 
  if (!is.null(sd)) {
    if (length(sd) < 3) {
      rf <- rnorm(2, 0, sd)
    } else {
      V <- matrix(c(sd[1]^2, sd[3], sd[3], sd[2]^2), nr = 2, nc = 2)
      rf <- as.numeric(rmvn(n = 1, mu = c(0, 0), V = V))
    }
    
  } else {
    rf <- c(0, 0)
  }
  diveI <- odiveI * exp(rf[1])
  surfI <- osurfI * exp(rf[2])
  while (t < T) {
    #cat("t = ", t, " T = ", T, "\r")
    Lsurf <- cumsum(surfI[tgr >= t]) * dt
    cdf_surf <- 1 - exp(-Lsurf)
    u_dive <- runif(1)
    t_dive <- tgr[tgr >= t][which.min(abs(cdf_surf - u_dive))]
    if (t_dive < t + 1e-10) t_dive <- t + dt / 2
    if (t_dive > T - dt) break 
    Ldive <- cumsum(diveI[tgr >= t_dive]) * dt 
    cdf_dive <- 1 - exp(-Ldive)
    u_surf <- runif(1)
    t_surf <- tgr[tgr >= t_dive][which.min(abs(cdf_dive - u_surf))]
    if (t_surf < t_dive + 1e-10) t_surf <- t_dive + dt / 2
    if (t_surf > T - dt) break
    dat <- rbind(dat, data.frame(ID = 1, dive = t_dive - t, surface = t_surf - t_dive, time = t))
    t <- t_surf
    if (!is.null(sd)) {
      if (length(sd) < 3) {
        rf <- rnorm(2, 0, sd)
      } else {
        rf <- as.numeric(rmvn(n = 1, mu = c(0, 0), V = V))
      }
    }
    diveI <- odiveI * exp(rf[1])
    surfI <- osurfI * exp(rf[2])
  }
  dat <- dat[-1,]
  return(dat)
}


