# Simulate data
library(CTMCdive)
library(expm)

## hacky simulation code
T <- 1000
dt <- 0.01
t <- seq(0, T, by = dt)
n <- length(t)
s1 <- cos(2 * pi * t / 1000)
s1 <- s1 - mean(s1)
diveI <- exp(log(0.05) - s1)
s2 <- 2.3 * s1
surfI <- exp(log(0.08) + s2)

#diveI <- rep(0.05, length(diveI))
#surfI <- rep(0.08, length(surfI))

#plot(t, diveI, type = "l")
#plot(t, surfI, type = "l")

s <- rep(0, n)
s[1] <- 1
dat <- data.frame(ID = 1, dive = 0, surface = 0, time = 1)
cur <- 1
diving <- TRUE
for (i in 2:n) {
  trm <- matrix(c(-surfI[i], diveI[i], surfI[i], -diveI[i]), nc = 2)
  tpm <- expm(trm*dt)
  s[i] <- sample(1:2, size = 1, prob = tpm[s[i - 1],])
  if (s[i] == 1) {
    dat$dive[cur] <- dat$dive[cur] + dt
  }
  if (s[i] == 2) {
    dat$surface[cur] <- dat$surface[cur] + dt
    diving <- FALSE
  }
  if (s[i] == 1 & !diving) {
    dat <- rbind(dat, data.frame(ID = 1, dive = 0, surface = 0, time = i*dt))
    cur <- cur + 1
    diving <- TRUE
  }

}


# setup model
forms <- list(surface ~ s(time, bs="cs"),
              dive ~ s(time, bs="cs"))

# fit model
simdat <- dat
mod <- FitCTMCdive(forms, simdat, model = "sd", print = TRUE)

mod

# predict durations
pred <- predict(mod)

# plot fit
plot(mod)
