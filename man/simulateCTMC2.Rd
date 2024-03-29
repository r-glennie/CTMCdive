% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{simulateCTMC2}
\alias{simulateCTMC2}
\title{Simulate CTMCdive data}
\usage{
simulateCTMC2(
  dive_fn,
  surf_fn,
  T,
  dt,
  kappa = list(dive = 1, surf = 1),
  sd = NULL,
  print = TRUE
)
}
\arguments{
\item{dive_fn}{function that given a vector of times computes dive intensity at those times}

\item{surf_fn}{function that given a vector of times computes surface intensity at those times}

\item{dt}{time step for simulation}

\item{n}{number of dive-surfacings to simulate}
}
\value{
CTMC dataset
}
\description{
Simulate CTMCdive data
}
