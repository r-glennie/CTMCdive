% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctmc_dive.R
\name{predict.CTMCdive}
\alias{predict.CTMCdive}
\title{Predict mean duration from fitted CTMC model for dives and surfacing}
\usage{
\method{predict}{CTMCdive}(x, newdata = NULL, ...)
}
\arguments{
\item{x}{CTMCdive fitted model x}

\item{newdata}{new data frame to make predictions for}

\item{\dots}{unused (for S3 compatability)}
}
\value{
a list of surface and dive mean predictions
}
\description{
Predict mean duration from fitted CTMC model for dives and surfacing
}
