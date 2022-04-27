# rddm R package 

The goal of rddm is to provide a user-friendly syntax in R to work with diffusion models of decision-making.
rddm is centered around two models: the drift diffusion model (DDM; e.g., Ratcliff et al., 1978, [Bogacz et al., 2006](https://pubmed.ncbi.nlm.nih.gov/17014301/))
and the pulse diffusion model. 

- Drift Diffusion Model and it's variants:
  - Pure DDM: using only drift rate, boundary, non-decision time, and starting point)
  - Extended DDM: The pure DDM + trial-to-trial variabliity in the drift rate, non-decision time, and starting point
  - Collapsing Bounds: the bounds collapse over time according to a linear, hyperbolic ratio, or weibull function
  - Urgency Signal: a linear or logistic urgency signal
<br>
<br>
- Pulse Diffusion Model: an extension of the DDM that considers moment-by-moment, time-varying evidence, 
such as a train of pulses of sound or light. This model is adapted from [Brunton et. al., 2013](https://www.science.org/doi/10.1126/science.1233912).
As implemented here, the pulse diffusion model is the exact same weiner diffusion process as the DDM, 
but the drift is a function of the momentary stimulus or evidence towards the upper and lower boundary.
Important notes:
  - This model uses the same parameter convention as the DDM, with three exceptions.
  Trial-to-trial variabliity in the drift rate is replaced with moment-to-moment variability in the drift rate.
  Trial-to-trial variability in the starting point uses a normal distribution, whereas it is uniform in the DDM.
  - This is actually an OU-process, not just an weiner diffusion process. There is an additional parameter, labmda. 
  For lambda = 0, this model is equal to a weiner diffusion process, but for lambda != 0, it is an OU-process.
  - This model offers the same collapsing bounds and urgency paramters as the DDM, but it is not recommended to use 
  these parameters if you are also using allowing lambda to vary.
  
  
### Installation

`rddm` can be installed directly from github using the `remotes` package:
```
# first, check for remotes package
if (!require(remotes)) install.packages("remotes")

# install rddm
remotes::install_github("gkane26/rddm", build_vignettes=TRUE)
```

#### Note :: OpenMP is used to speed up diffusion model simulations.

For best performance, please make sure OpenMP support is enabled. On Mac, it is recommended to use the gcc compiler installed via [homebrew](https://brew.sh/). Install homebrew, then run in terminal: `brew install gcc`. To use `gcc` as your compiler for R packages, specify the following in your `~/.R/Makevars` file:
```
VER=11
CC=gcc-$(VER)
CXX=g++-$(VER)
CXX11=g++-$(VER)

FLIBS=-L/usr/local/opt/gcc/lib/gcc/$(VER)

SHLIB_OPENMP_CFLAGS=-fopenmp
SHLIB_OPENMP_CXXFLAGS=-fopenmp
```
Check your version of gcc by running `brew info gcc` in terminal.


### Using rddm

Please see the `getting-started` vignette. After installation, in your R session run: `vignette("getting-started", "rddm")`.

