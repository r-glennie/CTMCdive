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
surf_I <- function(t, exp = FALSE) {
  eff <- 0.02 * (sin(2 * pi * t / T) + 1.2) + 0.05
  a <- exp_T - 24 * 60 
  b <- exp_T + 12 * 60
  if(exp) eff <- ifelse(t > exp_T, eff - 0.05 * (t - a)/(exp_T - a) * (t - b)/(exp_T - b) * (t > a) * (t < b), eff)
  return(eff)
}

# kappa
kappa <- list(dive = 3, surf = 3)

# plot truth 
plot(tgr, dive_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Dive Intensity")
plot(tgr, surf_I(tgr), type = "l", lwd = 1.5, xlab = "Time", ylab = "Surface Intensity")
lines(tgr, surf_I(tgr, exp = FALSE), col = "steelblue")
abline(v = c(exp_T, exp_T + 24 * 60), col = "firebrick", lty = "dashed")

# setup sims
set.seed(15839)
nsims <- 100
select <- rep(0, 2)
mods <- vector(mode = "list", length = nsims)
pb <- txtProgressBar(min = 1, max = nsims, initial = 0, style = 3)
sim <- 1 
while (sim <= nsims) {
  err <- FALSE
  setTxtProgressBar(pb, sim)
  # simulate data
  dat <- simulateCTMC2(dive_I, surf_I, T, dt = 0.01, kappa = kappa)
  
  # add exposure data
  dat$expt <- ifelse(dat$time >= exp_T & dat$time < exp_T + 24 * 60, dat$time - exp_T, 0)
  dat$exp <- ifelse(dat$time >= exp_T & dat$time < exp_T + 24 * 60, 1, 0)
  dat$expf <- factor(dat$exp, ordered = TRUE)
  
  # Basic model
  f0 <- list(surface ~ s(time, bs = "ts"),
             dive ~ s(time, bs = "ts"))
  m0 <- try(FitCTMCdive(f0, dat, dt = 1, print = FALSE))
  if ("try-error" %in% class(m0)) {sim <- sim - 1; err <- TRUE}
  forms <- list(surface ~ s(time, bs = "ts") + expf + s(expt, by = expf, bs ="cs", k = 10, m = 1),
                dive ~ s(time, bs = "ts"))
  mexp <- try(FitCTMCdive(forms, dat, dt = 1, print = FALSE))
  if ("try-error" %in% class(mexp)) {sim <- sim - 1; err <- TRUE}
  mods[[sim]] <- list(m0 = m0, mexp = mexp)
  aic <- AIC(m0, mexp)
  if (anyNA(aic)) {
    sim <- sim - 1 
    next
  }
  if (aic["mexp", 2] < aic["m0", 2]- 2) {
    expeff <- try(GetExposureEff(mexp, exp_var = "expf"))
    if ("try-error" %in% class(expeff)) {
      sim <- sim - 1; err <- TRUE
    } else {
      sig <- (expeff$dive$ci[1,] * expeff$dive$ci[2,]) > 0
      if (any(sig)) {
        modno <- 2
      } else {
        modno <- 1 
      }
    }
  } else {
    modno <- 1 
  }
  if(!err) select[modno] <- select[modno] + 1  
  sim <- sim + 1
  cat(select, "\n")
}


