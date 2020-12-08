# Joshua Alley
# Set up R script for project testing the public goods theory of alliances


# set seed
set.seed(12)



# load packages 
library(arm)
library(MASS)
library(reshape2)
library(plm)
library(texreg)
library(stargazer)
library(tidyverse)
library(rstan)
library(shinystan)
library(bayesplot)
library(gridExtra)
library(hexbin)
library(conflicted)


# set-up global STAN options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


# manage conflicts
conflict_scout()

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("lead", "dplyr")
conflict_prefer("extract", "rstan")
conflict_prefer("monitor", "rstan")
conflict_prefer("traceplot", "rstan")
conflict_prefer("display", "xtable")
conflict_prefer("expand", "tidyr")
conflict_prefer("melt", "reshape2")
conflict_prefer("scale_discrete_manual", "ggplot2")
conflict_prefer("coefplot", "texreg")
conflict_prefer("between", "dplyr")
conflict_prefer("chol2inv", "Matrix")
conflict_prefer("rcond", "Matrix")
conflict_prefer("pack", "tidyr")
conflict_prefer("unpack", "tidyr")
conflict_prefer("position", "ggplot2")
conflict_prefer("matches", "dplyr")
conflict_prefer("is_null", "purrr")


# set up functions 
# Autocorrelation of ratio-simulation
calcrho <- function(rho, rho1, rho2) {
    rho*(1-rho1*rho2)/ sqrt((1-rho1^2)*(1-rho2^2))
}
