# Simulate data
library(CTMCdive)

# Simulate data -----------------------------------------------------------

# constant hazards
dive_I <- function(t) {return(rep(0.06, length(t)))}
surf_I <- function(t) {return(rep(0.04, length(t)))}

# set kappa
kappa <- list(dive = 5, surf = 5)

# set individual-level random effects
#sd <- c(0.2, 0.5) # individual-level
sd <- c(1, 0.5, 0.7) # correlated random effects

# total observation time 
T <- 60 * 24 * 7

# time step 
dt <- 0.1

# simulate data
dat <- simulateCTMC2(dive_I, surf_I, T, dt, kappa = kappa, sd = sd)

# plot data
plot(dat$dive, pch = 19, xlab = "Time of Dive Start", ylab = "Dive Duration")
plot(dat$surf, pch = 19, xlab = "Time of Dive Start", ylab = "Surface Duration")
plot(dat$dive, dat$surf)

# Fit Model ---------------------------------------------------------------

# setup model
forms <- list(surface ~ 1,
              dive ~ 1)
# fit model
m <- FitCTMCdive(forms, dat, print = TRUE)
m_ind <- FitCTMCdive(forms, dat, print = TRUE, re = "ind")
m_corr <- FitCTMCdive(forms, dat, print = TRUE, re = "corr")

AIC(m, m_ind, m_corr)

# see results
mod <- m_corr
plot(mod)
exp(mod$res$surface[,1])
exp(mod$res$dive[,1])

# check residuals
pred <- predict(mod)
rdive <- pred$rdive
hist(rdive)
qqnorm(rdive);qqline(rdive)
ks.test(rdive, "pnorm")

rsurf <- pred$rsurf
hist(rsurf)
qqnorm(rsurf);qqline(rsurf)
ks.test(rsurf, "pnorm")


