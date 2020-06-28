library(R.matlab)
library(fields)
library(colorRamps)
library(rstan)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
options(mc.cores = parallel::detectCores())
library(brms)

read.csv('D.csv')

prior_go           <- get_prior(y ~ x + (x|region/id_X),data=D)
prior_go$prior[1]  <- 'uniform(-10,10)'
prior_go$prior[2]  <- 'uniform(-3,0)'
prior_go$prior[8]  <- 'uniform(1E-10,10)'
prior_go$prior[10] <- 'uniform(1E-10,3)'
prior_go$prior[11] <- 'uniform(1E-10,10)'
prior_go$prior[13] <- 'uniform(1E-10,3)'
	
nchains <- 4
niter   <- 2000	
	
fit <- brm(y ~ x + (x|region/id_X),data=D,iter=niter, chains=nchains, prior=prior_go)

