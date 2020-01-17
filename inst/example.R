# Simulate data 
library(CTMCdive)
library(expm)

T <- 1000
dt <- 0.1
t <- seq(0, T, by = dt)
n <- length(t)
diveI <- 0.1 * (-cos(2 * pi * t / 1000) + 1)
surfI <- 0.2 * (-cos(2 * pi * t / 1000) + 1)

diveI <- rep(0.05, length(diveI))
surfI <- rep(0.08, length(surfI))

plot(t, diveI, type = "l")

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

forms <- list(surface ~ 1, 
              dive ~ 1)

# fit model 
mod <- FitCTMCdive(forms, dat, model = "iid", print = TRUE)

mod


#pred <- mod$sm$A_grid %*% mod$rep$par.random[1:9]
#pred <- pred + mod$rep$par.fixed[1]
#pred <- exp(pred)
#plot(seq(dat$time[1], dat$time[nrow(dat)] + dat$dive[nrow(dat)] + dat$surface[nrow(dat)], length = 1000), pred)
#lines(t, diveI, col = "red")
