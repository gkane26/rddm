% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ddm_integral_fpt}
\alias{ddm_integral_fpt}
\title{Simulate drift diffusion model with fixed or collapsing boundary}
\usage{
ddm_integral_fpt(
  v,
  a,
  t0,
  z = 0.5,
  dc = 0,
  sv = 0,
  sz = 0,
  st0 = 0,
  aprime = 0,
  kappa = 0,
  tc = 0.25,
  s = 1,
  sv_points = 19,
  sz_points = 11,
  st0_points = 11,
  dt = 0.01,
  max_time = 10,
  bounds = 0L,
  n_threads = 1L
)
}
\arguments{
\item{v}{numeric; drift rate}

\item{a}{numeric; initial boundary}

\item{t0}{numeric; non-decision time}

\item{z}{numeric; starting point, 0 < z < 1, default = .5}

\item{dc}{numeric; drift criterion, the zero point of the drift rate (the drift rate v = v + dc); default = 0}

\item{sv}{numeric; standard deviation of variability in drift rate, sv >= 0, default = 0}

\item{sz}{numeric; variability in starting point. Uniform from [z-sz/2, z+sz/2], 0 < sz < z, default = 0}

\item{st0}{numeric; variability in non-decision time. Uniform from [t0-st0/2, t0+st0/2], 0 < st0 < t0, default = 0}

\item{aprime}{numeric; degree of collapse, default = 0}

\item{kappa}{numeric; slope of collapse, default = 0}

\item{tc}{numeric; time constant of collapse, default = .25}

\item{s}{numeric; standard deviation in wiener diffusion noise, default = 1}

\item{sv_points}{integer; number of points to approximate integral over drift rate variability, default = 19}

\item{sz_points}{integer; number of points to approximate integral over starting point variability, default = 11}

\item{st0_points}{integer; number of points to approximate integral over non-decision time variability, default = 11}

\item{dt}{numeric; time step of simulation, default = .01}

\item{max_time}{numeric; max time of simulation, default = 10}

\item{bounds}{int: 0 for fixed, 1 for hyperbolic ratio collapsing bounds, 2 for weibull collapsing bounds}

\item{n_threads}{integer; number of threads to run in parallel, default = 1}
}
\value{
data frame with two columns: response (1 for upper boundary, 0 for lower), and response time
}
\description{
Simulate drift diffusion model with fixed or collapsing boundary
}
