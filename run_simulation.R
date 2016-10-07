### ------------------------------------------------------------
## R script for the paper:
## Monnahan, C.C., Thorson, J.T., and Branch, T.A. (in press). "Faster
## estimation of Bayesian models in ecology using Hamiltonian Monte
## Carlo". Methods in Ecology and Evolution.

## Distinct models are contained in folders under the 'models' folder, with
## _nc appended in the case of both centered and non-centered (nc)
## models. Within each model folder is a run_model.R script, the model
## files (.stan and .jags), and the empirical data (data.RDS). For
## simulated models there is also a script to generate data. The fitted
## model objects are also saved (appended with the seed) in the 'fits'
## folder, and diagnostic and quick exploratory plots in the 'plots'
## folder.

## There are three "modes" for running this model. The first is to "verify"
## that the Stan and JAGS models are matching (to help debug). To do this
## the code runs both models with a high thinning rate so that samples are
## nearly independent and can be compared via QQ plots. Here we dont care
## about timing or efficiency. This takes a long time to run and only needs
## to be done once, so verify is set to FALSE below for each model, but
## should be run for any new model added here. The resulting samples from
## the Stan run are saved as sims.ind.RDS, and used as starting values in
## the empirical mode (below).

## The second mode is to compare efficiencies for real data. In this case
## there is no thinning and timing is important. Multiple repicates are
## run, with each software platform starting from the same place for each
## replicate/seed. Each model run in this mode gets saved to the results
## folder, with "empirical" appended.

## Last, there is simulated data for some (but not all) of the models. This
## behaves just like empirical data but for data generated inside.

## See README.md for more information.

### ------------------------------------------------------------
### Step 1: prepare working space; load libraries, functions, and global
### variables
## Set the reference directory from which everything else is relative.
main.dir <- 'C:/Users/Cole/gradmcmc/'
setwd(main.dir)
source("startup.R")
Nreps <- 20                 # number of replicates
Nout.ind <- 1000            # number of independent samples if verify==TRUE
set.seed(241)
seeds <- sample(1:1e5, size=Nreps)         #
## Old setting for testing static HMC vs NUTS. Leave as NULL since broken.
lambda.vec <- NULL
## Which metric to use for NUTS. Paper only used estimated diagonal (Stan
## default) but others could be used.
metric <- c('unit_e', 'diag_e', 'dense_e')[2]
## Suppress JAGS and Stan output to file via sink? Useful after debugging
## to clean up console and judge progress.
sink <- TRUE
## The paper was run with these key software versions
version$version.string                  # R version 3.2.3
packageVersion('rstan')                 # 2.11.1
packageVersion('R2jags')                # 0.5.7
packageVersion('rjags')                 # 4.4
### End of Step 1.
### ------------------------------------------------------------

### ------------------------------------------------------------
### Step 2: Run the models.
## Run multivariate normal, empirical and simulated
m <- 'mvnd'                             # model name
verify <- FALSE                         # whether to verify
delta <- 0.8                            # adapt_delta for Stan
Nout <- 20000                           # no. of samples out
Nthin <- 1                              # thin rate for emp/sim modes (leave at 1!)
Nthin.ind <- 100                        # thin rate for verify mode
## Settings for simulation mode. cor is a factor for independent (0) or
## from wishart (1) (see paper). Npar is how many parameters.
cor.vec <- c(0,1)
Npar.vec <- c(2, 5, 15, 25, 50)
source(paste0('models/',m,'/run_model.R'))

## Run MVN with varying correlations and a fixed Npar
m <- 'mvnc'
delta <- 0.8
Npar <- 5
Nout <- 20000; Nthin <- 1; Nthin.ind <- 100
cor.vec <- c(0, .25, .5, .75, .85, .9, .95)
Npar.vec <- c(2, 25, 50)
source(paste0('models/',m,'/run_model.R'))

## Run simulated growth tests, cross between centered/noncentered
Nout <- 20000; Nthin <- 1; Nthin.ind <- 500
Npar.vec <- c(2:10, floor(10^(seq(1.25, 2.5, by=.25))))
delta <- 0.8
m <- 'growth'
source(paste0('models/',m,'/run_model.R'))
delta <- 0.8
m <- 'growth_nc'
source(paste0('models/',m,'/run_model.R'))

## State space logistic. Note had to raise delta here.
Nout <- 20000; Nthin <- 1; Nthin.ind <- 500
delta <- 0.95
m <- 'ss_logistic'
source(paste0('models/',m,'/run_model.R'))
m <- 'ss_logistic_nc'
delta <- 0.95
source(paste0('models/',m,'/run_model.R'))

## Red kite example from Kery and Schaub; 8.4 w/ informative prior
Nout <- 20000; Nthin <- 1; Nthin.ind <- 100
delta <- 0.8
m <- 'redkite'
source(paste0('models/',m,'/run_model.R'))

## swallows -- Modified from example 14.5 from Korner-Nievergelt et al
Nout <- 20000; Nthin <- 1; Nthin.ind <- 100
delta <- .8
m <- 'swallows_nc'
source(paste0('models/',m,'/run_model.R'))
delta <- .8
m <- 'swallows'
source(paste0('models/',m,'/run_model.R'))

## Wildflower; from Bolker et al 2013
Nout <- 20000; Nthin <- 1; Nthin.ind <- 500
delta <- .8
m <- 'wildflower_nc'
source(paste0('models/',m,'/run_model.R'))
delta <- .8
m <- 'wildflower'
source(paste0('models/',m,'/run_model.R'))
### End of Step 2.
### ------------------------------------------------------------

### ------------------------------------------------------------
### Step 3: Load and prepare result data frames for plotting and tables
source('load_data.R')
source('make_plots.R')
### End of Step 3.
### ------------------------------------------------------------

