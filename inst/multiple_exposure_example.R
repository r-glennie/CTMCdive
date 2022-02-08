library(CTMCdive)

# Simulate data -----------------------------------------------------------

# total observation time 
T <- 24 * 60 * 7 

# exposure time 
exp_T <- c(T / 4, 3 * T / 4)

# time step 
dt <- 0.01

# time-varying intensities  
tgr <- seq(0, T, by = dt)
dive_I <- function(t) {
  eff <- 0.06 * (sin(2 * pi * t / T - pi / 2)) + 6 * 0.05
  return(eff)
}
surf_I <- function(t, exp = TRUE) {
  eff <- 0.02 * (sin(2 * pi * t / T) + 1.2) + 0.06
  a <- exp_T[1] - 24 * 60 
  b <- exp_T[1] + 24 * 60
  if(exp) eff <- ifelse(t > exp_T[1], eff - 0.05 * (t - a)/(exp_T[1] - a) * (t - b)/(exp_T[1] - b) * (t > a) * (t < b), eff)
  a <- exp_T[2] - 24 * 60 
  b <- exp_T[2] + 12 * 60
  if(exp) eff <- ifelse(t > exp_T[2], eff - 0.05 * (t - a)/(exp_T[2] - a) * (t - b)/(exp_T[2] - b) * (t > a) * (t < b), eff)
  return(eff)
}

# kappa 
kappa <- list(dive = 3, surf = 3)

# plot truth 
plot(tgr, dive_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Dive Intensity")
lines(tgr, dive_I(tgr), col = "steelblue")
plot(tgr, surf_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Surface Intensity")
abline(v = c(exp_T[1], exp_T[1] + 24 * 60), col = "firebrick", lty = "dashed")
abline(v = c(exp_T[2], exp_T[2] + 24 * 60), col = "firebrick", lty = "dashed")

# simulate data
seed <- sample(1:65555, size = 1)
set.seed(seed)
dat <- simulateCTMC2(dive_I, surf_I, T, dt, kappa = kappa)

# add exposure data
dat$exp1 <- ifelse(dat$time >= exp_T[1] & dat$time < exp_T[2], 1, 0)
dat$exp2 <- ifelse(dat$time >= exp_T[2], 1, 0)
dat$expt1 <- (dat$time - exp_T[1]) * dat$exp1
dat$expt2 <- (dat$time - exp_T[2]) * dat$exp2


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
m0 <- FitCTMCdive(forms, dat, dt = 1, print = TRUE, exp_time = c("expt1", "expt2"))

# fit exposure effect
mexp1 <- update(m0, ~.+ s(expt1, by = exp1, bs = "bs", k = 10, m = c(2, 1)))

mexp2 <- update(mexp1$surf, ~.+ s(expt2, by = exp2, bs = "bs", k = 10, m = c(2, 1)))


# select model
mod <- mexp2$both

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
expeff1 <- GetExposureEff(mod, exp_var = "exp1")
plotExposureEffect(expeff1, pick = "dive")
abline(v = c(exp_T[1], exp_T[1] + 24 * 60))
expeff2 <- GetExposureEff(mod, exp_var = "exp2")
plotExposureEffect(expeff2, pick = "dive")
abline(v = c(exp_T[2], exp_T[2] + 24 * 60))

# plot exposure against baseline prediction 
limits <- c(exp_T[1] - 24*60, exp_T[1] + 24 * 3 * 60)
plot(dat$time, dat$dive, pch = 19, xlab = "Time of Dive Start", ylab = "Surface Duration", xlim = limits)
abline(v = c(exp_T[1], exp_T[1] + 24 * 60), col = "firebrick", lty = "dashed")
lines(dat$time[dat$time < exp_T[1]], pred$dive[dat$time < exp_T[1]], col = "steelblue", lwd = 1.5)
lines(dat$time[dat$time > exp_T[1] + 24 * 60], pred$dive[dat$time > exp_T[1] + 24 * 60], col = "steelblue", lwd = 1.5)
matlines(expeff1$time, expeff1$dive$pred[,1:20], col = "firebrick", lwd = 1.5, alpha = 0.7)
matlines(expeff1$time, expeff1$dive$base[,1:20], col = "steelblue")

limits <- c(exp_T[2] - 24*60, exp_T[2] + 24 * 1 * 60)
plot(dat$time, dat$dive, pch = 19, xlab = "Time of Dive Start", ylab = "Surface Duration", xlim = limits)
abline(v = c(exp_T[2], exp_T[2] + 24 * 60), col = "firebrick", lty = "dashed")
lines(dat$time[dat$time < exp_T[2]], pred$dive[dat$time < exp_T[2]], col = "steelblue", lwd = 1.5)
lines(dat$time[dat$time > exp_T[2] + 24 * 60], pred$dive[dat$time > exp_T[2] + 24 * 60], col = "steelblue", lwd = 1.5)
matlines(expeff2$time, expeff2$dive$pred[,1:20], col = "firebrick", lwd = 1.5, alpha = 0.7)
matlines(expeff2$time, expeff2$dive$base[,1:20], col = "steelblue")

# plot estimated intensities against truth
plot(tgr, dive_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Dive Intensity")
lines(mod$sm$ints, pred$diveI, lwd = 1.5, col = "firebrick")
plot(tgr, surf_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Surface Intensity")
lines(mod$sm$ints, pred$surfI, lwd = 1.5, col = "firebrick")
lines(tgr, surf_I(tgr, exp = FALSE), col = "steelblue")



