% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{bounds}
\alias{bounds}
\alias{hyperbolic_ratio_bound}
\alias{weibull_bound}
\title{Collapsing boundary functions}
\usage{
hyperbolic_ratio_bound(t, a, kappa = 0, tc = 0.25)

weibull_bound(t, a, a_prime = 0, kappa = 1, tc = 0.25)
}
\arguments{
\item{t}{vector; time points to evaluate boundary}

\item{a}{numeric; initial boundary}

\item{kappa}{numeric; slope of collapse}

\item{tc}{numeric; time constant of collapse}

\item{a_prime}{numeric; degree of collapse (weibull only)}
}
\value{
boundary at time t
}
\description{
evaluate diffusion model boundary using the hyperbolic ratio or weibull functions
}
