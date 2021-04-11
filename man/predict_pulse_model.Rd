% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pulse_diffusion_model.R
\name{predict_pulse_model}
\alias{predict_pulse_model}
\title{predict pulse model  (for internal use)}
\usage{
model$predict(pars=NULL, n=10000, ...)
}
\arguments{
\item{pars}{numeric vector; vector of parameters. If NULL, uses model$solution$pars.}

\item{n}{integer; number of decisions to simulate for each condition. If the number of conditions is equal to the length of the data, e.g. if using as_function with a continuous predictor, ignores \code{n} and simulates one decision per condition}

\item{method}{string; "euler" for euler-maruyama simulation or "fp" for Fokker-Planck method}

\item{...}{additional arguments passed to method (either \code{sim_pulse} or \code{pulse_fp_fpt})}
}
\value{
data.table with simulation conditions, decision (upper or lower boundary) and response time
}
\description{
Predict behavior with given pulse model parameters.
This function is only intended for use with a diffusion model object,
and should not be called directly outside of the diffusion model class.
}
\keyword{internal}