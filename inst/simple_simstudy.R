library(CTMCdive)

# Setup -------------------------------------------------------------------
set.seed(46644)
nsims <- 100

# total observation time 
T <- 24 * 60 * 7

# time step
dt <- 0.1

# constant intensities  
dive_I <- function(t) {return(rep(0.06, length(t)))}
surf_I <- function(t) {return(rep(0.04, length(t)))}

# set kappa
kappa <- list(dive = 5, surf = 5)

# setup model
forms <- list(surface ~ 1,
              dive ~ 1)

# space to save fitted models
mods <- vector(mode = "list", length = nsims)

# Simstudy ----------------------------------------------------------------

for (sim in 1:nsims) {
  cat(sim, " / ", nsims, "\r")
  
  # simulate data
  dat <- simulateCTMC2(dive_I, surf_I, T, dt, kappa = kappa, print = FALSE)
  
  # fit model
  mods[[sim]] <- FitCTMCdive(forms, dat, dt = 1, print = FALSE)
  
}


# Results -----------------------------------------------------------------

dive <- sapply(mods, FUN = function(x) {exp(x$res$dive$`Est.`)})
hist(dive)
summary(dive)

surf <- sapply(mods, FUN = function(x) {exp(x$res$surface$`Est.`)})
hist(surf)
summary(surf)

kappa_dive <- sapply(mods, FUN = function(x) {exp(x$rep$par.fixed["log_kappa_dive"])})
kappa_surf <- sapply(mods, FUN = function(x) {exp(x$rep$par.fixed["log_kappa_surf"])})

hist(kappa_dive)
hist(kappa_surf)
