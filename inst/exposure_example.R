library(CTMCdive)

# Simulate data -----------------------------------------------------------

# total observation time 
T <- 24 * 60 * 7 

# exposure time 
exp_T <- T / 2

# time step 
dt <- 0.01

# time-varying intensities  
tgr <- seq(0, T, by = dt)
dive_I <- function(t, exp = TRUE, spike = TRUE) {
  eff <- 0.06 * (sin(2 * pi * t / T - pi / 2)) + 6 * 0.03
  a <- exp_T - 24 * 60 
  b <- exp_T + 12 * 60
  if(exp & !spike) eff <- ifelse(t > exp_T, eff - 0.15 * (t - a)/(exp_T - a) * (t - b)/(exp_T - b) * (t > a) * (t < b), eff)
  if(exp & spike) eff <- ifelse(t > exp_T, eff - 0.15 * (t >= exp_T) * (t <= exp_T + 100), eff)
  return(eff)
}
surf_I <- function(t) {
  #return(rep(0.05, length(t)))
  return(0.02 * (sin(2 * pi * t / T) + 1.2))
}

# kappa 
kappa <- list(dive = 3, surf = 3)

# plot truth 
plot(tgr, dive_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Dive Intensity")
lines(tgr, dive_I(tgr, exp = FALSE), col = "steelblue")
abline(v = c(exp_T, exp_T + 24 * 60), col = "firebrick", lty = "dashed")
plot(tgr, surf_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Surface Intensity")

# simulate data
seed <- sample(1:65555, size = 1)
set.seed(seed)
dat <- simulateCTMC2(dive_I, surf_I, T, dt, tstart = exp_T, kappa = kappa)

# add exposure data
dat$expt <- ifelse(dat$time >= exp_T & dat$time < exp_T + 24 * 60, dat$time - exp_T, 0)
dat$expf <- factor(ifelse(dat$time >= exp_T & dat$time < exp_T + 24 * 60, 1, 0), ordered = TRUE)
dat$spike <- factor(ifelse(dat$time >= exp_T & dat$time <= exp_T + 1e-5, 1, 0))

# plot data
plot(dat$time, dat$dive, pch = 19, xlab = "Time of Dive Start", ylab = "Dive Duration")
abline(v = c(exp_T, exp_T + 24 * 60), col = "firebrick", lty = "dashed")
plot(dat$time, dat$surf, pch = 19, xlab = "Time of Dive Start", ylab = "Surface Duration")
abline(v = c(exp_T, exp_T + 24 * 60), col = "firebrick", lty = "dashed")

# Fit Model ---------------------------------------------------------------

# setup model
forms <- list(surface ~ s(time, bs = "cs") + + s(time, by = expf, bs = "ts", m = 1) + s(expf, bs = "re"),
              dive ~ s(time, bs = "cs") + s(time, by = expf, bs = "ts", m = 1) + s(expf, bs = "re"))
# fit model
mod_sm <- FitCTMCdive(forms, dat, dt = 1, print = TRUE)

# setup model
forms <- list(surface ~ s(time, bs = "cs") + spike,
              dive ~ s(time, bs = "cs") + spike)
# fit model
mod_spike <- FitCTMCdive(forms, dat, dt = 1, print = TRUE)

# Response 
AIC(mod_sm, mod_spike)

# Best model
mod <- mod_spike

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

# test residuals for normality 
ks.test(rsurf, "pnorm")
ks.test(rdive, "pnorm")

# estimated exposure effect
exp_var <- "spike"
expeff <- GetExposureEff(mod, exp_var = exp_var)
plotExposureEffect(expeff, pick = "dive")
plotExposureEffect(expeff, pick = "surface")

# plot exposure against baseline prediction 
limits <- c(exp_T - 1 * 24 * 60, exp_T + 2 * 24 * 60)
plot(dat$time, dat$surf, pch = 19, xlab = "Time of Dive Start", ylab = "Surface Duration", xlim = limits)
abline(v = c(exp_T, exp_T + 24 * 60), col = "firebrick", lty = "dashed")
lines(dat$time[dat$time < exp_T], pred$surface[dat$time < exp_T], col = "steelblue", lwd = 1.5)
lines(dat$time[dat$time > exp_T + 24 * 60], pred$surface[dat$time > exp_T + 24 * 60], col = "steelblue", lwd = 1.5)
matlines(expeff$time, expeff$surf$pred[,1:20], col = "firebrick", lwd = 1.5, alpha = 0.7)
matlines(expeff$time, expeff$surf$base[,1:20], col = "steelblue")

# plot estimated intensities against truth
plot(tgr, dive_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Dive Intensity", ylim = c(0, 0.25))
lines(mod$sm$ints, pred$diveI, lwd = 1.5, col = "firebrick")
lines(tgr, dive_I(tgr, exp = FALSE), col = "steelblue")
plot(tgr, surf_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Surface Intensity")
lines(mod$sm$ints, pred$surfI, lwd = 1.5, col = "firebrick")




