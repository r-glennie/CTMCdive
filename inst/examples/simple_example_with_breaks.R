# Simulate data
library(CTMCdive)

# Simulate data -----------------------------------------------------------

# constant hazards
dive_I <- function(t) {return(rep(0.06, length(t)))}
surf_I <- function(t) {return(rep(0.04, length(t)))}

# set kappa
kappa <- list(dive = 5, surf = 5)

# total observation time 
T <- 60 * 24 * 7

# time step 
dt <- 0.1

# simulate data
dat <- simulateCTMC2(dive_I, surf_I, T, dt, kappa = kappa)

# make breaks
olddat <- dat
dat <- dat[-c(30:50, 100:130),]

# set breaks 
breaks <- data.frame(start = c(olddat$time[30], olddat$time[100]), 
                     end = c(olddat$time[51], olddat$time[131]))

# plot data
plot(dat$time, dat$dive, pch = 19, xlab = "Time of Dive Start", ylab = "Dive Duration")
abline(v = c(breaks$start, breaks$end), lty = "dotted", col = "steelblue")
plot(dat$time, dat$surf, pch = 19, xlab = "Time of Dive Start", ylab = "Surface Duration")
abline(v = c(breaks$start, breaks$end), lty = "dotted", col = "steelblue")

# Fit Model ---------------------------------------------------------------

# setup model
forms <- list(surface ~ 1,
              dive ~ 1)
# fit model
mod <- FitCTMCdive(forms, dat, breaks = breaks, print = TRUE)

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


