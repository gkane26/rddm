% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drift_diffusion_model.R
\name{simulate_diffusion_model}
\alias{simulate_diffusion_model}
\title{simulate diffusion model  (for internal use)}
\usage{
model$simulate(n, par_values, par_names=NULL, ...)
}
\arguments{
\item{n}{integer; number of decisions to simulate for each condition. If the number of conditions is equal to the length of the data, e.g. if using as_function with a continuous predictor, ignores \code{n} and simulates one decision per condition}

\item{pars}{numeric vector; vector of parameters}

\item{...}{additional arguments passed to \code{sim_ddm}}
}
\value{
data.table with simulation conditions, decision (upper or lower boundary) and response time
}
\description{
simulate DDM with given diffusion model parameters..
This function is only intended for use with a diffusion model object,
and should not be called directly outside of the diffusion model class.
Please use \code{sim_ddm} as a standalone function.
}
\keyword{internal}
