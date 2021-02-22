
#' Simulate CTMCdive data
#'
#' @param dive_fn function that given a vector of times computes dive intensity at those times
#' @param surf_fn function that given a vector of times computes surface intensity at those times
#' @param dt time step for simulation 
#'
#' @return CTMC dataset 
#' @export
simulateCTMC <- function(dive_fn, surf_fn, T, dt, print = TRUE) {
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
  # fill out data
  for (i in 2:nt) {
    if(print) cat(round(100 * (i / nt), 0), "% done \r", sep = "")
    trm <- matrix(c(-surfI[i], diveI[i], surfI[i], -diveI[i]), nc = 2)
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
simulateCTMC2 <- function(dive_fn, surf_fn, T, dt, print = TRUE) {
  tgr <- seq(0, T, by = dt)
  # predict intensity on fine grid 
  diveI <- dive_fn(tgr)
  surfI <- surf_fn(tgr)
  # current time 
  t <- 0 
  # create empty dataset
  dat <- data.frame(ID = 1, dive = 0, surface = 0, time = 1)
  cur <- 1 
  while (t < T) {
    #cat("t = ", t, " T = ", T, "\r")
    Lsurf <- cumsum(surfI[tgr >= t]) * dt
    cdf_surf <- 1 - exp(-Lsurf)
    u_dive <- runif(1)
    t_dive <- tgr[tgr >= t][which.min(abs(cdf_surf - u_dive))]
    if (t_dive > T - dt) break 
    Ldive <- cumsum(diveI[tgr >= t_dive]) * dt 
    cdf_dive <- 1 - exp(-Ldive)
    u_surf <- runif(1)
    t_surf <- tgr[tgr >= t_dive][which.min(abs(cdf_dive - u_surf))]
    if (t_surf > T - dt) break
    dat <- rbind(dat, data.frame(ID = 1, dive = t_dive - t, surface = t_surf - t_dive, time = t))
    t <- t_surf
  }
  dat <- dat[-1,]
  return(dat)
}


