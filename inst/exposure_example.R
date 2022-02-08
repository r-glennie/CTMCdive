library(CTMCdive)

# Simulate data -----------------------------------------------------------

# total observation time 
T <- 24 * 60 * 7 

# exposure time 
exp_T <- T / 2

# time step 
dt <- 1

# time-varying intensities  
tgr <- seq(0, T, by = dt)
dive_I <- function(t) {
  eff <- 0.06 * (sin(2 * pi * t / T - pi / 2)) + 6 * 0.05
  return(eff)
}
surf_I <- function(t, exp = TRUE) {
  eff <- 0.02 * (sin(2 * pi * t / T) + 1.2) + 0.06
  a <- exp_T - 24 * 60 
  b <- exp_T + 24 * 60
  if(exp) eff <- ifelse(t > exp_T, eff - 0.03 * (t - a)/(exp_T - a) * (t - b)/(exp_T - b) * (t > a) * (t < b), eff)
  return(eff)
}

# kappa 
kappa <- list(dive = 3, surf = 3)

# plot truth 
plot(tgr, dive_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Dive Intensity")
lines(tgr, dive_I(tgr), col = "steelblue")
plot(tgr, surf_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Surface Intensity")
abline(v = c(exp_T, exp_T + 24 * 60), col = "firebrick", lty = "dashed")

# simulate data
seed <- sample(1:65555, size = 1)
set.seed(seed)
dat <- simulateCTMC2(dive_I, surf_I, T, dt, kappa = kappa)

# add exposure data
#dat$expt <- ifelse(dat$time >= exp_T & dat$time < exp_T + 24 * 60, dat$time - exp_T, 0)
#dat$exp <- ifelse(dat$time >= exp_T & dat$time < exp_T + 24 * 60, 1, 0)
dat$expt <- ifelse(dat$time >= exp_T, dat$time - exp_T, 0)
dat$exp <- ifelse(dat$time >= exp_T, 1, 0)

# plot data
plot(dat$time, dat$dive, pch = 19, xlab = "Time of Dive Start", ylab = "Dive Duration")
abline(v = c(exp_T, exp_T + 24 * 60), col = "firebrick", lty = "dashed")
plot(dat$time, dat$surf, pch = 19, xlab = "Time of Dive Start", ylab = "Surface Duration")
abline(v = c(exp_T, exp_T + 24 * 60), col = "firebrick", lty = "dashed")

# Fit Model ---------------------------------------------------------------

# setup model
forms <- list(surface ~ s(time, bs = "ts"),
              dive ~ s(time, bs = "ts"))
# fit model
m0 <- FitCTMCdive(forms, dat, dt = 1, print = TRUE, exp_time = "expt")

# fit exposure effect
mexp <- update(m0, ~.+ s(expt, by = exp, bs = "bs", k = 10, m = c(2, 1)))

# select model
mod <- mexp$surf

# see results
summary(mod)
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

# estimated exposure effect
expeff <- GetExposureEff(mod, exp_var = "exp")
plotExposureEffect(expeff, pick = "dive")
abline(v = c(exp_T, exp_T + 24 * 60))

# plot exposure against baseline prediction 
limits <- c(exp_T - 24*60, exp_T + 3*24*60)
plot(dat$time, dat$dive, pch = 19, xlab = "Time of Dive Start", ylab = "Surface Duration", xlim = limits)
abline(v = c(exp_T, exp_T + 24 * 60), col = "firebrick", lty = "dashed")
lines(dat$time[dat$time < exp_T], pred$dive[dat$time < exp_T], col = "steelblue", lwd = 1.5)
lines(dat$time[dat$time > exp_T + 24 * 60], pred$dive[dat$time > exp_T + 24 * 60], col = "steelblue", lwd = 1.5)
matlines(expeff$time, expeff$dive$pred[,1:20], col = "firebrick", lwd = 1.5, alpha = 0.7)
matlines(expeff$time, expeff$dive$base[,1:20], col = "steelblue")

# plot estimated intensities against truth
plot(tgr, dive_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Dive Intensity")
lines(mod$sm$ints, pred$diveI, lwd = 1.5, col = "firebrick")
plot(tgr, surf_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Surface Intensity")
lines(mod$sm$ints, pred$surfI, lwd = 1.5, col = "firebrick")
lines(tgr, surf_I(tgr, exp = FALSE), col = "steelblue")



