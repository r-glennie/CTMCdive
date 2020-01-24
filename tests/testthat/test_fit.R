# test FitCTMCdive

context("FitCTMCdive")

# load data etc
library(CTMCdive)

load(system.file("simdat.RData", package="CTMCdive"))

test_that("surface~s(time), dive~s(time)",{
  forms <- list(surface ~ s(time, bs="cs"),
                dive ~ s(time, bs="cs"))

  # fit model
  mod <- FitCTMCdive(forms, simdat, print = FALSE)

  expect_equal(mod$res$surface$Est, -1.841388492)
  expect_equal(mod$res$dive$Est, -3.073709637)
})


test_that("surface~1, dive~s(time)",{
  forms <- list(surface ~ 1,
                dive ~ s(time, bs="cs"))

  # fit model
  mod <- FitCTMCdive(forms, simdat, print = FALSE)

  expect_equal(mod$res$surface$Est, -2.90144483)
  expect_equal(mod$res$dive$Est, -3.073735243)
})


test_that("surface~s(time), dive~1",{
  forms <- list(surface ~ s(time, bs="cs"),
                dive ~ 1)

  # fit model
  mod <- FitCTMCdive(forms, simdat, print = FALSE)

  expect_equal(mod$res$surface$Est, -1.841374046)
  expect_equal(mod$res$dive$Est, -3.073729866)
})





