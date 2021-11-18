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
  eff <- 0.06 * (sin(2 * pi * t / T - pi / 2)) + 6 * 0.03
  a <- exp_T - 24 * 60 
  b <- exp_T + 12 * 60
  if(exp) eff <- ifelse(t > exp_T, eff - 0.15 * (t - a)/(exp_T - a) * (t - b)/(exp_T - b) * (t > a) * (t < b), eff)
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

# setup sims
set.seed(15839)
nsims <- 100
select <- rep(0, 2)
mods <- vector(mode = "list", length = nsims)
pb <- txtProgressBar(min = 1, max = nsims, initial = 0, style = 3)
for (sim in 1:nsims) {
  setTxtProgressBar(pb, sim)
  # simulate data
  dat <- simulateCTMC2(dive_I, surf_I, T, dt, tstart = exp_T, kappa = kappa)
  
  # add exposure data
  dat$expt <- ifelse(dat$time >= exp_T & dat$time < exp_T + 24 * 60, dat$time - exp_T, 0)
  dat$expf <- factor(ifelse(dat$time >= exp_T & dat$time < exp_T + 24 * 60, 1, 0), ordered = TRUE)
  
  # Basic model
  f0 <- list(surface ~ s(time, bs = "cs"),
             dive ~ s(time, bs = "cs"))
  m0 <- FitCTMCdive(f0, dat, dt = 1, print = FALSE)
  forms <- list(surface ~ s(time, bs = "cs"),
                dive ~ s(time, bs = "cs") + s(time, by = expf, bs = "ts", m = 1) + s(expf, bs = "re"))
  mexp <- FitCTMCdive(forms, dat, dt = 1, print = FALSE)
  mods[[sim]] <- list(m0 = m0, mexp = mexp)
  aic <- AIC(m0, mexp)
  if (aic["mexp", 2] <= aic["m0", 2] - 2) {
    expeff <- GetExposureEff(mexp, exp_var = "expf")
    sig <- (expeff$surf$ci[1,] * expeff$surf$ci[2,]) > 0
    if (any(sig)) {
      modno <- 2
    } else {
      modno <- 1 
    }
  } else {
    modno <- 1 
  }
  select[modno] <- select[modno] + 1  
  cat(select, "\n")
}


