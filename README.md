This repository contains code for the paper:

Monnahan, C.C., J.T. Thorson, T.A. Branch (In press). Faster estimation of
Bayesian models in ecology using Hamiltonian Monte Carlo. Methods in
Ecology and Evolution.

This study can be reproduced by installing necessary software and R
packages, and then running the script run_simulation.R. See that file for
more information.

The models directory contains case studies, with data (as .RDS files), Stan
and JAGS implementations (as .stan and .jags files), and an R script to run
them (using some functions defined in startup.R). Within each model there
is a folder containing plots of results.

The algorithms directory contains an R script containing versions of the
static HMC and NUTS algorithms implemented in R. These versions do not
contain the full suite of features used in Stan, but may be useful for
learning more about their behavior. The script mcmc.R contains the
algorithms and helper functions, and examples.R a short script showing how
to use them.
