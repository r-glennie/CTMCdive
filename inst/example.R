# Simulate data
library(CTMCdive)
library(mgcv)
#devtools::load_all()
library(expm)

## hacky simulation code
set.seed(10987)
T <- 10000
dt <- 0.5
t <- seq(0, T, by = dt)
n <- length(t)

exposure <- FALSE
if (exposure) {
  exp <- rep(0, n)
  exp[t > 5000 & t < 6000] <- 1 
  exp <- as.factor(exp)
  gm <- gam(int ~ s(t, by = exp), data = data.frame(int = 1, t = t, exp = exp), fit = FALSE)
  X <- gm$X
  betaexp <- c(-10, -10, -10*0.5, 0, 2 * 0.5, -1, -1, -0.5, 1) * c(0.3,0.36, 0.82,-0.5,0.93,0.21,0.2,-0.02,0.38)
  beta <- c(-3,0.3,0.36,0.82,-0.5,0.93,0.21,0.2,-0.02,0.38, betaexp)
  s1 <- X %*% beta 
  plot(t, s1, type = "l", xlab = "time", ylab = "dive intensity")
  abline(v = c(5000, 6000), col = "red", lty = "dotted")
  
  beta2 <- c(-2.53,-0.18,-0.16,1.36,1.37,0.82,0.66,0.78,-0.99,0.18,-0.18,-0.16,1.36,1.37,0.82,0.66,0.78,-0.99,0.18)
  s2 <- X %*% beta2
  plot(t, s2, type = "l")
} else {
  gm <- gam(int ~ s(t), data = data.frame(int = 1, t = t), fit = FALSE)
  X <- gm$X
  beta <- c(-3,0.3,0.36,0.82,-0.5,0.93,0.21,0.2,-0.02,0.38)
  s1 <- X %*% beta
  plot(t, s1, type = "l")
  
  beta2 <- c(-2.53,-0.18,-0.16,1.36,1.37,0.82,0.66,0.78,-0.99,0.18)
  s2 <- X %*% beta2
  plot(t, s2, type = "l", xlab = "time", ylab = "surface intensity")
}

diveI <- exp(s1)
surfI <- exp(s2)

plot(t, diveI, type = "l")
plot(t, surfI, type = "l")

s <- rep(0, n)
s[1] <- 1
dat <- data.frame(ID = 1, dive = 0, surface = 0, time = 1)
cur <- 1
diving <- TRUE
for (i in 2:n) {
  cat(i, " / ", n, "\r")
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

# plot data 
plot(dat$time + 0.5 * dat$dive, dat$dive, type = "b", pch = 19, xlab = "time", ylab = "dive duration")
abline(v = c(5000, 6000), col = "red", lty = "dotted")
plot(dat$time + dat$dive + 0.5 * dat$surface, dat$surface, type = "b", pch = 19, xlab ="time", ylab = "surface duration")
abline(v = c(5000, 6000), col = "red", lty = "dotted")

# setup model
forms <- list(surface ~ s(time, bs="cs"),
              dive ~ s(time, bs="cs"))
# fit model
mod <- FitCTMCdive(forms, dat, print = TRUE)

mod 

tgrid <- seq(5000, 6000, 1)
subdat <- dat[dat$time > 5000 & dat$time <= 6000,]
freq <- as.numeric(table(sapply(tgrid, FUN = function(x) which.min(abs(x - subdat$time)))))
predgrid <- subdat[rep(seq(nrow(subdat)), freq), ] 
predgrid$time <- tgrid
predgrid$exp <- 1 

predgrid <- ExpandCovs(dat, tgrid)
expeff <- GetExposureEff(mod, predgrid)
plot(expeff)

# predict durations
pred <- predict(mod)

# plot fit
plot(mod)

rand <- mod$rep$par.random
Xdive <- mod$sm$A_grid_dive
Xsurf <- mod$sm$A_grid_surf

preddive <- exp (mod$rep$par.fixed["par_dive"] + Xdive %*% rand[names(rand) == "s_dive"])
plot(mod$sm$ints, preddive, type = "l", ylim = c(min(diveI), max(diveI)), xlab = "time", ylab = "dive intensity")
lines(t, diveI, col = "blue")
predsurf <- exp (mod$rep$par.fixed["par_surf"] + Xdive %*% rand[names(rand) == "s_surf"])
plot(mod$sm$ints, predsurf, type = "l", xlab= "time", ylab="surface intensity")
lines(t, surfI, col = "blue")

## exposure model
dat$exp <- rep(0, nrow(dat))
dat$exp[dat$time > 5000 & dat$time < 6000] <- 1
dat$exp <- as.factor(dat$exp)

# setup model
forms <- list(surface ~ s(time, bs="cs"),
              dive ~ s(time, bs="cs", by = exp))
# fit model
modexp  <- FitCTMCdive(forms, dat, print = TRUE)

modexp

expeff <- GetExposureEff(modexp, predgrid)
plot(expeff)

tgrid <- seq(5000, 6000, 1)
predgrid <- data.frame(ID = 1, dive = 1, surface = 1, time = tgrid, exp = 1)
expeff <- GetExposureEff(modexp, predgrid)
plot(expeff)

plot(tgrid, expeff$dive$mean, type = "l", ylim = range(expeff$dive$ci))
lines(tgrid, expeff$dive$ci[1,], lty = "dashed")
lines(tgrid, expeff$dive$ci[2,], lty = "dashed")
abline(h = 0, col = "red", lty = "dotted")
wh <- min(which(expeff$dive$ci[1,] < 0 & expeff$dive$ci[2,] > 0))
abline(v = tgrid[wh], col = "blue", lty = "dotted")

AIC(mod, modexp)

# plot fit
plot(modexp)

rand <- modexp$rep$par.random
fixed <- modexp$rep$par.fixed

gmdive <- modexp$sm$gam_dive
Xdive <- predict(gmdive, data.frame(ID = 1, dive = 1, surface = 1, time = t, exp = exp), type = "lpmatrix")
pardive <- c(fixed[names(fixed) == "par_dive"], rand[names(rand) == "s_dive"])
preddive <- exp(Xdive %*% pardive)
plot(t, preddive, type = "l", ylim = c(min(diveI), max(diveI)), xlab = "time", ylab = "dive intensity")
lines(t, diveI, col = "blue")

Xbase <- predict(gmdive, data.frame(ID = 1, dive = 1, surface = 1, time = t, exp = 0), type = "lpmatrix")
basedive <- exp(Xbase  %*% pardive)
lines(t, basedive, col = "red")
V <- solve(modexp$rep$jointPrecision)
betas <- rmvn(1000, c(modexp$rep$par.fixed, modexp$rep$par.random), V)
betadive <- betas[, c(names(fixed) == "par_dive", names(rand) == "s_dive")]
simdive <- exp(Xdive %*% t(betadive))
simbase <- exp(Xbase %*% t(betadive))

quantbase <- apply(simbase, 1, quantile, prob = c(0.025, 0.975))
lines(t, quantbase[1,], col = "red", lty = "dashed")
lines(t, quantbase[2,], col = "red", lty = "dashed")


simexp <- simdive - simbase
simexp <- simexp[t > 5000 & t < 6000,]
meanexp <- rowMeans(simexp)
expt <- t[t>5000 & t < 6000]
quant <- apply(simexp, 1, quantile, prob = c(0.025, 0.975))
plot(expt, meanexp, type = "l", ylim = c(min(quant), max(quant)), xlab = "time", ylab = "exposure effect")
lines(expt, quant[1,], lty = "dotted")
lines(expt, quant[2,], lty = "dotted")
abline(h = 0, col = "red", lty = "dotted")

## plot exposure effect on dive 

GetExposureEff <- function(mod, predgrid, nsims = 1000) {
  basedat <- predgrid 
  basedat$exp <- 0 
  Xdive <- predict(mod$sm$gam_dive, predgrid, type = "lpmatrix")
  Xdivebase <- predict(mod$sm$gam_dive, basedat, type = "lpmatrix")
  Xsurf <- predict(mod$sm$gam_surf, predgrid, type = "lpmatrix")
  Xbasesurf <- predict(mod$sm$gam_surf, basedat, type = "lpmatrix")
  V <- solve(mod$rep$jointPrecision)
  betas <- rmvn(nsims, c(mod$rep$par.fixed, mod$rep$par.random), V)
  betadive <- t(betas[, c(names(fixed) == "par_dive", names(rand) == "s_dive")])
  betasurf <- t(betas[, c(names(fixed) == "par_surf", names(rand) == "s_surf")])
  simdive <- exp(Xdive %*% betadive)
  simbasedive <- exp(Xdivebase %*% betadive)
  simsurf <- exp(Xsurf %*% betasurf)
  simbasesurf <- exp(Xbasesurf %*% betasurf)
  simdiveexp <- simdive - simbasedive
  simsurfexp <- simsurf - simbasesurf
  simdiveexp <- simdiveexp[predgrid$exp == 1,]
  simsurfexp <- simsurfexp[predgrid$exp == 1,]
  dive <- surf <- NULL
  dive$mean <- rowMeans(simdiveexp)
  dive$ci <- apply(simdiveexp, 1, quantile, prob = c(0.025, 0.975))
  surf$mean <- rowMeans(simsurfexp)
  surf$ci <- apply(simsurfexp, 1, quantile, prob = c(0.025, 0.975))
  res <- list(dive = dive, surf = surf, predgrid = predgrid)
  class(res) <- "ExposureEffect"
  return(res) 
}

plot.ExposureEffect <- function(expeff, pick = "both") {
  if (pick == "dive" | pick == "both") {
    plot(expeff$predgrid$time, expeff$dive$mean, type = "l", ylim = range(expeff$dive$ci), xlab = "Time", ylab = "Dive intensity exposure effect")
    lines(expeff$predgrid$time, expeff$dive$ci[1,], lty = "dashed")
    lines(expeff$predgrid$time, expeff$dive$ci[2,], lty = "dashed")
    abline(h = 0, col = "red", lty = "dotted")
    wh <- which(expeff$dive$ci[1,] < 0 & expeff$dive$ci[2,] > 0)
    if (length(wh) > 0) {
      wh <- min(wh)
      abline(v = expeff$predgrid$time[wh], col = "blue", lty = "dotted")
      cat("Return to baseline for dive: ", expeff$predgrid$time[wh] - min(expeff$predgrid$time[expeff$predgrid$exp == 1]), "\n")
    } else {
      cat("No estimated return to baseline for dive\n")
    }
  } 
  if (pick == "surf" | pick == "both") {
    plot(expeff$predgrid$time, expeff$surf$mean, type = "l", ylim = range(expeff$surf$ci), xlab = "Time", ylab = "Surface intensity exposure effect")
    lines(expeff$predgrid$time, expeff$surf$ci[1,], lty = "dashed")
    lines(expeff$predgrid$time, expeff$surf$ci[2,], lty = "dashed")
    abline(h = 0, col = "red", lty = "dotted")
    wh <- which(expeff$surf$ci[1,] < 0 & expeff$surf$ci[2,] > 0)
    if (length(wh) > 0) {
      wh <- min(wh)
      abline(v = expeff$predgrid$time[wh], col = "blue", lty = "dotted")
      cat("Return to baseline for surface: ", expeff$predgrid$time[wh] - min(expeff$predgrid$time[expeff$predgrid$exp == 1]), "\n")
    } else {
      cat("No estimated return to baseline for surface\n")
    }
  }
}

## try double mod

forms <- list(surface ~ s(time, bs="cs", by = exp),
              dive ~ s(time, bs="cs", by = exp))
# fit model
modexp2  <- FitCTMCdive(forms, dat, print = TRUE)

AIC(mod, modexp, modexp2)
