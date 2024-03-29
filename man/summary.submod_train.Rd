% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/submod_train.R
\name{summary.submod_train}
\alias{summary.submod_train}
\title{Subgroup Identification (Summary)}
\usage{
\method{summary}{submod_train}(object, round_est = 4, round_SE = 4, round_CI = 4, ...)
}
\arguments{
\item{object}{Trained submod_train model.}

\item{round_est}{Rounding for trt ests (default=4)}

\item{round_SE}{Rounding for trt SEs (default=4)}

\item{round_CI}{Rounding for trt CIs (default=4)}

\item{...}{Any additional parameters, not currently passed through.}
}
\value{
List of key outputs (1) Number of Identified Subgroups, and (2) Treatment effect estimates, 
SEs, and CIs for each subgroup/estimand
}
\description{
Summary for subgroup identification function.
}
