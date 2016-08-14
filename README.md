This repository contains code for a manuscript about Hamiltonian Monte
Carlo algorithms for Bayesian models in ecology.

This study can be reproduced by installing necessary software and R
packages, and then running the script run_simulation.R. 

The models directory contains case studies, with data (as .RDS files), Stan
and JAGS implementations, and an R script to run them (using some functions
defined in startup.R). Within each model there is a folder containing
plots of results.

The algorithms directory contains an R script containing versions of the
static HMC and NUTs algorithms. These versions do not contain the full
suite of features used in Stan, but may be useful for learning more about
their behavior. The script mcmc.R contains the algorithms and helper
functions, and examples.R a short script showing how to use them.
