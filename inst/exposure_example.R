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
dive_I <- function(t, exp = TRUE) {
  eff <- 0.05 * (cos(2 * pi * t / T) + 1.2)
  if(exp) eff <- ifelse(t > exp_T, eff + 0.1 * exp(-(t - exp_T) / 500), eff)
  return(eff)
}
surf_I <- function(t) {
  return(0.02 * (sin(2 * pi * t / T) + 1.2))
}
# plot truth 
plot(tgr, dive_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Dive Intensity")
lines(tgr, dive_I(tgr, exp = FALSE), col = "steelblue")
abline(v = c(exp_T, exp_T + 24 * 60), col = "firebrick", lty = "dashed")
plot(tgr, surf_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Surface Intensity")

# simulate data
set.seed(sample(1:65555, size = 1))
dat <- simulateCTMC(dive_I, surf_I, T, dt)

# add exposure data
dat$expt <- ifelse(dat$time >= exp_T & dat$time < exp_T + 24 * 60, dat$time - exp_T, 0)
dat$expf <- factor(ifelse(dat$time >= exp_T & dat$time < exp_T + 24 * 60, 1, 0))

# plot data
plot(dat$time, dat$dive, pch = 19, xlab = "Time of Dive Start", ylab = "Dive Duration")
abline(v = c(exp_T, exp_T + 24 * 60), col = "firebrick", lty = "dashed")
plot(dat$time, dat$surf, pch = 19, xlab = "Time of Dive Start", ylab = "Surface Duration")
abline(v = c(exp_T, exp_T + 24 * 60), col = "firebrick", lty = "dashed")

# Fit Model ---------------------------------------------------------------

# setup model
forms <- list(surface ~ s(time, bs = "cs"),
              dive ~ s(time, bs = "cs", by = expf))
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


# estimated exposure effect
tgrid <- seq(exp_T, exp_T + 24 * 60, 1)
predgrid <- data.frame(ID = 1, dive = 1, surface = 1, time = tgrid, expf = 1)
expeff <- GetExposureEff(mod, predgrid, exp = "expf")
plot(expeff, pick = "dive")
ends <- 5900

# plot estimated intensities against truth
plot(tgr, dive_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Dive Intensity")
lines(mod$sm$ints, pred$diveI, lwd = 1.5, col = "firebrick")
lines(tgr, dive_I(tgr, exp = FALSE), col = "steelblue")
abline(v = ends, lty = "dashed")
plot(tgr, surf_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Surface Intensity")
lines(mod$sm$ints, pred$surfI, lwd = 1.5, col = "firebrick")




