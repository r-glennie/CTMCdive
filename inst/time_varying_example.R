library(CTMCdive)

# Simulate data -----------------------------------------------------------

# total observation time 
T <- 24 * 60 * 7

# time step 
dt <- 1

# time-varying intensities  
tgr <- seq(0, T, by = dt)
dive_I <- function(t) {
  return(0.01 + 0.2 * (t/T - 1/2)^2)
}
surf_I <- function(t) {
  x <- t / T
  f <- 0.2 * x^11 * (10 * (1 - x))^6 + 10 * 
    (10 * x)^3 * (1 - x)^10
  return(0.02 + f / 150)
}

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
surfdurs <- sapply(tgr, surf_dur, divei = divei, tgr = tgr, dt = dt)
divedurs <- sapply(tgr, dive_dur, surfi = surfi, tgr = tgr, dt = dt)

# plot truth 
plot(tgr, dive_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Dive Intensity")
plot(tgr, surf_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Surface Intensity")

# simulate data
set.seed(sample(1:65555, size = 1))
dat <- simulateCTMC2(dive_I, surf_I, T, dt)

# plot data
plot(dat$dive, pch = 19, xlab = "Time of Dive Start", ylab = "Dive Duration")
plot(dat$surf, pch = 19, xlab = "Time of Dive Start", ylab = "Surface Duration")

# Fit Model ---------------------------------------------------------------

# setup model
forms <- list(surface ~ s(time, bs = "cs"),
              dive ~ s(time, bs = "cs"))
# fit model
mod <- FitCTMCdive(forms, dat, print = TRUE)

# see results
mod
exp(mod$res$surface[,1])
exp(mod$res$dive[,1])

# plot fitted model
plot(mod)

# get predicted values
pred <- predict(mod)

# plot residuals
rdive <- pred$rdive
plot(dat$time, rdive, pch = 19)
hist(rdive)
qqnorm(rdive); qqline(rdive)
rsurf <- pred$rsurf
plot(dat$time, rsurf, pch = 19)
hist(rsurf)
qqnorm(rsurf); qqline(rsurf)

# test residuals for normality 
ks.test(rsurf, "pnorm")
ks.test(rdive, "pnorm")

# plot estimated intensities against truth
plot(tgr, dive_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Dive Intensity")
lines(mod$sm$ints, pred$diveI, lwd = 1.5, col = "firebrick")
plot(tgr, surf_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Surface Intensity")
lines(mod$sm$ints, pred$surfI, lwd = 1.5, col = "firebrick")

# compare true and estimated mean durations 
estD <- pred$dive
trueD <- sapply(dat$time, dive_dur, surfi = surfi, tgr = tgr, dt = dt)
plot(estD, trueD)
abline(a = 0, b = 1)
plot(dat$time, estD)
lines(dat$time, trueD)

estS <- pred$surface
trueS <- sapply(dat$time + dat$dive, surf_dur, divei = divei, tgr = tgr, dt = dt)
plot(estS, trueS)
abline(a = 0, b = 1)
plot(dat$time + dat$dive, estS)
lines(dat$time + dat$dive, trueS)
