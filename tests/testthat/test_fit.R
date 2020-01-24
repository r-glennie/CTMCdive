# test FitCTMCdive

context("FitCTMCdive")

# load data etc
library(CTMCdive)

load(system.file("simdat.RData", package="CTMCdive"))

test_that("point estimates are correct",{
  forms <- list(surface ~ s(time, bs="cs"),
                dive ~ s(time, bs="cs"))

  # fit model
  mod <- FitCTMCdive(forms, simdat, print = FALSE)

  expect_equal(mod$res$surface$Est, -1.841388492)
  expect_equal(mod$res$dive$Est, -3.073709637)
})






