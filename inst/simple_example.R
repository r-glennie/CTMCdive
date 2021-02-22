# Simulate data
library(CTMCdive)

# Simulate data -----------------------------------------------------------

# constant hazards
dive_I <- function(t) {return(rep(0.06, length(t)))}
surf_I <- function(t) {return(rep(0.04, length(t)))}

# total observation time 
T <- 60 * 24 * 7

# time step 
dt <- 1

# simulate data
dat <- simulateCTMC(dive_I, surf_I, T, dt)

# plot data
plot(dat$dive, pch = 19, xlab = "Time of Dive Start", ylab = "Dive Duration")
plot(dat$surf, pch = 19, xlab = "Time of Dive Start", ylab = "Surface Duration")

# Fit Model ---------------------------------------------------------------

# setup model
forms <- list(surface ~ 1,
              dive ~ 1)
# fit model
mod <- FitCTMCdive(forms, dat, print = TRUE)

# see results
mod
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


