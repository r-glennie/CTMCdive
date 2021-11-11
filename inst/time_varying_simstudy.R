library(CTMCdive)

# Setup -------------------------------------------------------------------
set.seed(46644)
nsims <- 100

# total observation time 
T <- 24 * 60 * 7

# time step 
dt <- 0.01

#time-varying intensities  
tgr <- seq(0, T, by = dt)
dive_I <- function(t) {
  #return(rep(100, length(t)))
  return(0.01 + 0.2 * (t/T - 1/2)^2)
}
surf_I <- function(t) {
  x <- t / T
  f <- 0.2 * x^11 * (10 * (1 - x))^6 + 4 * 
    (10 * x)^3 * (1 - x)^10
  return(0.1 - f / 50)
}

# kappa
kappa <- list(dive = 3, surf = 3)

# mean durations given start time 
divei <- dive_I(tgr)
surfi <- surf_I(tgr)

surf_dur <- function(t, divei, tgr, dt) {
  surv <- exp(-cumsum(divei[tgr >= t - 1e-10]) * dt)
  est_duration <- sum(surv) * dt
  return(est_duration)
}
dive_dur <- function(t, surfi, tgr, dt) {
  surv <- exp(-cumsum(surfi[tgr >= t - 1e-10]) * dt)
  est_duration <- sum(surv) * dt
  return(est_duration)
}
#surfdurs <- sapply(tgr, surf_dur, divei = divei, tgr = tgr, dt = dt)
#divedurs <- sapply(tgr, dive_dur, surfi = surfi, tgr = tgr, dt = dt)

# plot truth 
plot(tgr, dive_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Dive Intensity")
plot(tgr, surf_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Surface Intensity")

# setup model
forms <- list(surface ~ s(time, bs = "cs"),
              dive ~ s(time, bs = "cs"))

# space to save fitted models
mods <- vector(mode = "list", length = nsims)

# Simstudy ----------------------------------------------------------------

pb <- txtProgressBar(min = 0, max = nsims, style = 3)
for (sim in 1:nsims) {
  # set progress
  setTxtProgressBar(pb, sim)
  
  # simulate data
  dat <- simulateCTMC2(dive_I, surf_I, T, dt = 1, kappa = kappa, print = FALSE)
  
  # fit model
  mods[[sim]] <- FitCTMCdive(forms, dat, print = FALSE)
  
}

# Predictions -------------------------------------------------------------
# compute predicted dive times 
preds <- lapply(mods, FUN = predict)
divepred <- sapply(preds, FUN = function(x) {x$diveI})
surfpred <- sapply(preds, FUN = function(x) {x$surfI})
ints <- mods[[1]]$sm$ints

# estimated intensities 
plot(tgr, dive_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Dive Intensity")
matlines(ints, divepred, col = "grey80")
lines(ints, rowMeans(divepred), col = "firebrick", lwd = 1.5)
plot(tgr, surf_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Surface Intensity")
matlines(ints, surfpred, col = "grey80")
lines(ints, rowMeans(surfpred), col = "firebrick", lwd = 1.5)

# predicted durations
estD <- sapply(preds, FUN = function(x) {x$dive})
estS <- sapply(preds, FUN = function(x) {x$surface})
trueD <- lapply(mods, FUN = function(m) {
    sapply(m$dat$time, dive_dur, surfi = surfi, tgr = tgr, dt = dt)
}) 
trueS <- lapply(mods, FUN = function(m) {
  sapply(m$dat$time + m$dat$dive, surf_dur, divei = divei, tgr = tgr, dt = dt)
})

# dive predicted durations
reldistD <- sapply(1:length(trueD), FUN = function(i) {
  reldist <- abs(estD[[i]] - trueD[[i]]) #/ trueD[[i]]
  reldist[1]
})
summary(reldistD)
hist(reldistD)

j <- 72
plot(mods[[j]]$dat$time, estD[[j]], col = "red")
points(mods[[j]]$dat$time, trueD[[j]])

# surface predicted durations
reldistS <- sapply(1:length(trueS), FUN = function(i) {
  reldist <- (estS[[i]] - trueS[[i]]) / trueS[[i]]
  100 * mean(reldist)
})
summary(reldistS)
hist(reldistS)

j <- 10
plot(mods[[j]]$dat$time, estS[[j]], col = "red", pch = 20)
points(mods[[j]]$dat$time, trueS[[j]], pch = 20)





