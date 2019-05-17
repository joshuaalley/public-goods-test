# Joshua Alley
# Texas A&M University
# Simulate fake data from ML model and attempt to recover parameters

# Initial idea from: http://modernstatisticalworkflow.blogspot.com/2017/04/an-easy-way-to-simulate-fake-data-from.html 


# working directory set through projects

# load packages 
library(MASS)
library(tidyverse)
library(rstan)
library(shinystan)
library(bayesplot)


# set-up global STAN options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


# set seed
set.seed(12)

# Set-up simulatiom data
n.sim = 5000 # 2500 simulated observations


# pull data into a list and add some simulations
sim.data.noest <- list(N = n.sim, # 2000 obs
                       y = rcauchy(n.sim, location = 0, scale = .25), # cauchy outcome
                       state = rep(1:50, times = n.sim/50), # 50 states
                       S = 50, 
                       year = rep(1:200, times = n.sim/200), # 200 years
                       T = 200,
                       A = 200, # 200 alliances
                       Z = state.mem.mat[1001:6000, 1:200], # real data- hard to simulate 
                       X = mvrnorm(n.sim, mu = c(1, -1), Sigma = matrix(c(4,2,2,1),2,2)), # simulate state-level IVs from multivariate normal dist
                       M = 2, # two state-level variables
                       run_estimation = 0
)


# Run STAN model 
compiled.ml <- stan_model("data/ml-model-sim.stan")

# Run model to generate draws from posterior and parameter values
sim.out <- sampling(compiled.ml, data = sim.data.noest,
                    iter = 1000, warmup = 500, chains = 4)

# Check diagnostics
check_hmc_diagnostics(sim.out)



# Pull parameter values and simulated data from one draw of posterior (12th)
draw <- 12
y.sim <- extract(sim.out, pars = "y_sim")[[1]][draw, ] # outcome 
true.beta <- extract(sim.out, pars = "beta")[[1]][draw, ] # beta estimates
true.gamma <- extract(sim.out, pars = "gamma")[[1]][draw, ] # gamma estimates
true.theta <- extract(sim.out, pars = "theta")[[1]][draw] # lambda estimates
true.sigma <- extract(sim.out, pars = "sigma")[[1]][draw] # sigma estimates
true.sigma.state <- extract(sim.out, pars = "sigma_state")[[1]][draw] # sigma state estimates
true.sigma.year <- extract(sim.out, pars = "sigma_year")[[1]][draw] # sigma year estimates
true.sigma.all <- extract(sim.out, pars = "sigma_all")[[1]][draw] # sigma alliance estimates



# Fit model on draws from generated quantities block. 
# keep same data on IVs and structure as generated above
sim.data.est <- sim.data.noest # copy list 
sim.data.est$y <- y.sim # replace simulated y with y from STAN draw
sim.data.est$run_estimation <- 1 # estimate the likelihood


# run the model on this simulated data: attempt to recover parameters
sim.out.est <- sampling(compiled.ml, data = sim.data.est,
                        iter = 1000, warmup = 500, chains = 4)


# Check diagnostics
check_hmc_diagnostics(sim.out.est)


### Plot estimated parameters against "true values" from earlier simulated data
sim.est.sum <- extract(sim.out.est, pars = c("beta", "gamma", "theta", "sigma", 
                                             "sigma_state", "sigma_year", "sigma_all"), 
                       permuted = TRUE)

colnames(sim.est.sum$beta) <- c("beta1", "beta2")
colnames(sim.est.sum$gamma) <- paste0('gamma', 1:200)


# Use mcmc_areas to plot credible intervals and overlap
# Set the color scheme 
color_scheme_set("gray")

# Start with beta- second-level regression parameters
mcmc_areas(sim.est.sum$beta, pars = c("beta1"), prob = .9) +
  vline_at(true.beta[1], color = "black", size = 2) 
mcmc_areas(sim.est.sum$beta, pars = c("beta2"), prob = .9) +
  vline_at(true.beta[2], color = "black", size = 2) 


# Gamma parameters
mcmc_areas(sim.est.sum$gamma, pars = c("gamma1"), prob = .9) +
  vline_at(true.gamma[1], color = "black", size = 2) 
mcmc_areas(sim.est.sum$gamma, pars = c("gamma100"), prob = .9) +
  vline_at(true.gamma[100], color = "black", size = 2) 
mcmc_areas(sim.est.sum$gamma, pars = c("gamma175"), prob = .9) +
  vline_at(true.gamma[175], color = "black", size = 2) 



# Calculate accuracy of credible intervals: how many intervals contain true gamma? 
gamma.summary.sim <- summary(sim.out.est, pars = c("gamma"), probs = c(0.05, 0.95))$summary
gamma.summary.sim <- cbind.data.frame(true.gamma, gamma.summary.sim)

# create a dummy indicator of accurate coverage
gamma.summary.sim$accurate <- ifelse(gamma.summary.sim$true.gamma > gamma.summary.sim$`5%` & # greater than lower bound
                                       gamma.summary.sim$true.gamma < gamma.summary.sim$`95%`,
                                     1, 0) # smaller than upper bound
sum(gamma.summary.sim$accurate)


# Create dataframe of parameters 
sim.est.hyper <- cbind.data.frame(sim.est.sum$theta, sim.est.sum$sigma,
                                  sim.est.sum$sigma_state, sim.est.sum$sigma_year,
                                  sim.est.sum$sigma_all)
colnames(sim.est.hyper) <- c("theta", "sigma", "sigma_state", 
                             "sigma_year", "sigma_all")

# then theta- overall mean
mcmc_areas(sim.est.hyper, pars = c("theta"), prob = .9) +
  vline_at(true.theta[1], color = "black", size = 2) 
ggsave("appendix/theta-sim-res.pdf", height = 6, width = 8)

# check whether model recovers variance parameter and hyperparameters
# use ggplot b/c mcmc_areas doesn't take numeric vector input
# first-level variance sigma
mcmc_areas(sim.est.hyper, pars = c("sigma"), prob = .9) +
  vline_at(true.sigma[1], color = "black", size = 2) 

# state variance hyperparameter
mcmc_areas(sim.est.hyper, pars = c("sigma_state"), prob = .9) +
  vline_at(true.sigma.state[1], color = "black", size = 2) 

# year variance hyperparameter
mcmc_areas(sim.est.hyper, pars = c("sigma_year"), prob = .9) +
  vline_at(true.sigma.year[1], color = "black", size = 2) 


# alliance coefficient variance parameter
mcmc_areas(sim.est.hyper, pars = c("sigma_all"), prob = .9) +
  vline_at(true.sigma.all[1], color = "black", size = 2) 
ggsave("appendix/sall-sim-res.pdf", height = 6, width = 8)
