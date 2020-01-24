# test FitCTMCdive

context("FitCTMCdive")

# load data etc
library(CTMCdive)

load(system.file("simdat.RData", package="CTMCdive"))

test_that("point estimates are correct",{
  forms <- list(surface ~ 1,
                dive ~ 1)

  # fit model
  mod <- FitCTMCdive(forms, simdat, model = "sd", print = FALSE)

  expect_equal(mod$res$surface$Est, -1.841388492)
  expect_equal(mod$res$dive$Est, -3.073709637)
})






