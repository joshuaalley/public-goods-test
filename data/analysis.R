# Joshua Alley
# Texas A&M University
# Empirical test of public goods theory of alliances


# Load packages
library(here)
library(arm)
library(reshape2)
library(MASS)
library(plm)
library(texreg)
library(interflex)
library(sampleSelection)
library(margins)
library(tidyverse)
library(rstan)
library(bayesplot)
library(shinystan)



# Set working directory to current folder 
setwd(here::here())
getwd()


# Set-up STAN guidelines
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


# Set seed 
set.seed(12)


# Load state-ally-year data 
state.ally.year <- read.csv("data/alliance-state-year.csv")
state.vars <- read.csv("data/state-vars.csv")



# summarize state-ally-year data with key indicators
state.ally.year$treaty.pres <- ifelse(state.ally.year$atopid > 0, 1, 0)

state.ally.sum <- state.ally.year %>%
  group_by(ccode, year) %>%
  summarize(
    treaty.count = n(),
    total.ally.expend = sum(ally.spend, na.rm = TRUE),
    avg.treaty.contrib = mean(alliance.contrib, na.rm = TRUE),
    avg.dem.prop = mean(avg.democ, na.rm = TRUE),
    avg.num.mem = mean(num.mem, na.rm = TRUE),
    treaty.pres = max(treaty.pres)
  )

# Sumarize key variables
summary(state.ally.sum$avg.treaty.contrib)
ggplot(state.ally.sum, aes(x = avg.treaty.contrib)) + geom_histogram()
summary(state.ally.sum$total.ally.expend)


# Create a full dataset
state.char.full <- left_join(state.vars, state.ally.sum)


# Fill missing values of alliance variables with zero
state.char.full[, 36: ncol(state.char.full)][is.na(state.char.full[, 36: ncol(state.char.full)])] <- 0




# Log-transform allied expenditures variable
state.char.full$ln.ally.expend <- log(state.char.full$total.ally.expend + 1)


# Difference allied spending: consistent with approach of Plumper and Neumayer
state.char.full <- state.char.full %>%
                  group_by(ccode) %>%
                 mutate(
                   lag.ally.expend = lag(ln.ally.expend),
                   diff.ally.expend = ln.ally.expend - lag.ally.expend
                 ) 


# Set up interaction with a separate dataset
state.char.inter <- as.data.frame(state.char.full) # interflex doesn't take tibble input
state.char.inter <- filter(state.char.inter, ln.gdp >= 18) # filter really tiny states out

# summarize interaction variables
summary(state.char.inter$ln.gdp)
ggplot(state.char.inter, aes(x = ln.gdp)) + geom_density()


summary(state.char.inter$total.ally.expend)
ggplot(state.char.inter, aes(x = total.ally.expend)) + geom_density()


summary(state.char.inter$ln.ally.expend)
ggplot(state.char.inter, aes(x = ln.ally.expend)) + geom_histogram()




### First test: aboslute size (GDP)
# Interact changes in allied spending and GDP
# Total allied spending: pooling regression
m1.pg.abs <- lm(ln.milex ~ diff.ally.expend + ln.gdp + diff.ally.expend:ln.gdp +
                     lag.ln.milex + avg.num.mem + avg.dem.prop + 
                     atwar + civilwar.part + polity  + 
                     lsthreat + cold.war,
                   data = state.char.inter
                )
summary(m1.pg.abs)
# Calculate marginal effects
margins(m1.pg.abs)
cplot(m1.pg.abs, x = "ln.gdp", dx = "diff.ally.expend", what = "effect",
      main = "Average Marginal Effect of Changes in Allied Spending")
abline(h = 0)



# FGLS 
m2.pg.abs <- pggls(ln.milex ~ diff.ally.expend + ln.gdp + diff.ally.expend:ln.gdp +
                       avg.dem.prop + lag.ln.milex +
                       atwar + civilwar.part + polity + ln.gdp + avg.num.mem +
                       lsthreat + cold.war,
                     data = state.char.inter, subset = (majpower == 0),
                     index = c("ccode", "year"),
                     effect = "individual", # unrestricted error covariance
                     model = "pooling")
summary(m2.pg.abs)


# binning estimator
bin.abs <- inter.binning(Y = "ln.milex", D = "diff.ally.expend", X = "ln.gdp", 
              Z = c("lag.ln.milex", "atwar", "civilwar.part", "polity", 
                    "lsthreat", "cold.war", "avg.num.mem", "avg.dem.prop"), 
              data = state.char.inter,
              na.rm = TRUE
)
bin.abs

# Kernel: 10+ minute run time 
kernel.abs <- inter.kernel(Y = "ln.milex", D = "diff.ally.expend", X = "ln.gdp", 
             Z = c("lag.ln.milex", "atwar", "civilwar.part", "polity", 
                   "lsthreat", "cold.war", "avg.num.mem", "avg.dem.prop"), 
             data = state.char.inter, 
             na.rm = TRUE,
             nboots = 200, parallel = TRUE, cores = 4
)
kernel.abs


# Check for non-random selection into alliances and compare allied states
# Still a statistically significant interaction
heckit.ally.spend <- heckit(selection = treaty.pres ~ lag.ln.milex + 
                              ln.gdp + polity + atwar + lsthreat,
                            outcome = ln.milex ~ diff.ally.expend + ln.gdp + 
                              diff.ally.expend:ln.gdp +
                              lag.ln.milex +
                              atwar + civilwar.part + polity + ln.gdp +
                              lsthreat + cold.war,
                            data = state.char.inter,
                            method = "ml")
summary(heckit.ally.spend)





### Second test: relative size expressed as contribution to alliance
# estimate interactions
# filter out cases with no alliances
inter.data.rel <- filter(state.char.full, avg.treaty.contrib > 0)
inter.data.rel <- as.data.frame(inter.data.rel)

# Total allied spending: pooling regression
m1.pg.rel <- lm(ln.milex ~ diff.ally.expend + avg.treaty.contrib + diff.ally.expend:avg.treaty.contrib +
                     lag.ln.milex + avg.num.mem + avg.dem.prop + 
                     atwar + civilwar.part + polity  + 
                     lsthreat + cold.war,
                   data = inter.data.rel)
summary(m1.pg.rel)
# Calculate marginal effects
margins(m1.pg.rel)
cplot(m1.pg.rel, x = "avg.treaty.contrib", dx = "diff.ally.expend", what = "effect",
      main = "Average Marginal Effect of Changes in Allied Spending")
abline(h = 0)

# FGLS 
m2.pg.rel <- pggls(ln.milex ~ diff.ally.expend + avg.treaty.contrib + diff.ally.expend:avg.treaty.contrib +
                       avg.dem.prop + lag.ln.milex +
                       atwar + civilwar.part + polity + ln.gdp + avg.num.mem +
                       lsthreat + cold.war,
                     data = inter.data.rel, subset = (majpower == 0),
                     index = c("ccode", "year"),
                     effect = "individual", # unrestricted error covariance
                     model = "pooling")
summary(m2.pg.rel)


# binning estimator
bin.rel <- inter.binning(Y = "ln.milex", D = "diff.ally.expend", X = "avg.treaty.contrib", 
              Z = c("lag.ln.milex", "atwar", "civilwar.part", "polity", 
                    "lsthreat", "cold.war", "avg.num.mem", "avg.dem.prop"), 
              data = inter.data.rel,
              na.rm = TRUE
)
bin.rel 

# Kernel: 10+ minute run time 
kernel.rel <- inter.kernel(Y = "ln.milex", D = "diff.ally.expend", X = "avg.treaty.contrib", 
             Z = c("lag.ln.milex", "atwar", "civilwar.part", "polity", 
                   "lsthreat", "cold.war", "avg.num.mem", "avg.dem.prop"), 
             data = inter.data.rel, 
             na.rm = TRUE,
             nboots = 200, parallel = TRUE, cores = 4
)
kernel.rel

# Check for non-random selection into alliances and compare allied states
# less evidence of interaction
heckit.rel <- heckit(selection = treaty.pres ~ lag.ln.milex + 
                              ln.gdp + polity + atwar + lsthreat,
                            outcome = ln.milex ~ diff.ally.expend + avg.treaty.contrib +
                              diff.ally.expend:avg.treaty.contrib +
                              lag.ln.milex +
                              atwar + civilwar.part + polity + ln.gdp +
                              lsthreat + cold.war,
                            data = state.char.full,
                            method = "ml")
summary(heckit.rel)




### Second approach: Estimate response of each state to every alliance they are in
# Requires estimating the model in STAN


# Load state-year matrix of alliance participation:
atop.cow.year <- read.csv("data/atop-cow-year.csv")

# Create a dataset of state-year alliance membership:
atop.cow.year <- group_by(atop.cow.year, atopid, ccode, year)
state.mem <- atop.cow.year %>% select(atopid, ccode, year)
state.mem <-  mutate(state.mem, member = 1)
state.mem <- distinct(state.mem, atopid, ccode, year, .keep_all = TRUE)

# This matrix has a binary indicator of which alliances states are a member of in a given year
state.mem <- spread(state.mem, key = atopid, value = member, fill = 0)

# Remove the zero or no alliance category
state.mem <- subset(state.mem, select = -(3))


# Add state membership in alliances to this data
reg.state.data <-  state.vars %>%
                   select(ccode, year, ln.milex, lag.ln.milex,
                          atwar, civilwar.part, rival.milex, ln.gdp, polity, 
                          cold.war, disputes, majpower) %>%
                  left_join(state.mem)


# Create a matrix of state membership in alliances (Z in STAN model)
reg.state.data <- reg.state.data[complete.cases(reg.state.data), ]
state.mem.mat <- as.matrix(reg.state.data[, 13: ncol(reg.state.data)])

# Rescale variables
reg.state.data[, 5:11] <- lapply(reg.state.data[, 5:11], 
                                 function(x) rescale(x, binary.inputs = "0/1"))

reg.state.mat <- as.matrix(reg.state.data[, 4:12])

# Set-up data for STAN
# create a state index variable
reg.state.data$state.id <- reg.state.data %>% group_indices(ccode)
# Create a year index variable 
reg.state.data$year.id <- reg.state.data %>% group_indices(year)



# Define the data list 
stan.data <- list(N = nrow(reg.state.data), y = reg.state.data[, 3],
                  state = reg.state.data$state.id, S = length(unique(reg.state.data$state.id)),
                  year = reg.state.data$year.id, T = length(unique(reg.state.data$year.id)),
                  A = ncol(state.mem.mat),
                  Z = state.mem.mat, 
                  X = reg.state.mat, M = ncol(reg.state.mat)
)

# Compile the model code
model.1 <- stan_model(file = "data/ml-model-stan.stan")

# Variational Bayes- use to check coefficients and posterior predictions on model
ml.model.vb <- vb(model.1, data = stan.data, seed = 12)
launch_shinystan(ml.model.vb)

# posterior predictive check from variational Bayes- did not converge
# so treat these predictions with caution
y = reg.state.data[, 3]
vb.model.sum <- extract(ml.model.vb)
ppc_dens_overlay(y, vb.model.sum$y_pred[1:100, ])



# Regular STAN
system.time(
  ml.model <- sampling(model.1, data = stan.data, 
                       iter = 2000, warmup = 1000, chains = 4
  )
)

# diagnose full model
launch_shinystan(ml.model)

check_hmc_diagnostics(ml.model)


# Extract coefficients from the model
ml.model.sum <- extract(ml.model, permuted = TRUE)

# Posterior predictive distributions relative to observed data
yrep.full <- ml.model.sum$y_pred

# plot posterior predictive denisty of first 100 simulations
ppc_dens_overlay(y, yrep.full[1:100, ])


# Summarize lamdba 
lambda.summary <- summary(ml.model, pars = c("lambda"), probs = c(0.05, 0.95))$summary
lambda.summary <- cbind.data.frame(as.numeric(colnames(state.mem.mat)), lambda.summary)
colnames(lambda.summary) <- c("atopid", "lambda.mean", "lambda.se.mean",
                              "lambda.sd", "lambda.5", "lambda.95",
                              "lambda.neff", "lambda.rhat")
lambda.summary$lambda.positive <- ifelse((lambda.summary$lambda.5 > 0 & lambda.summary$lambda.95 > 0), 1, 0)
sum(lambda.summary$lambda.positive)
lambda.summary$lambda.negative <- ifelse((lambda.summary$lambda.5 < 0 & lambda.summary$lambda.95 < 0), 1, 0)
sum(lambda.summary$lambda.negative)

# Plot posterior means of alliance coefficients
ggplot(lambda.summary, aes(x = lambda.mean)) +
  geom_density() +
  ggtitle("Posterior Means of Alliance Coefficients")

ggplot(lambda.summary, aes(x = lambda.mean)) +
  geom_histogram(bins = 60) +
  ggtitle("Posterior Means of Alliance Coefficients")



# Plot points with error bars by ATOPID
ggplot(lambda.summary, aes(x = atopid, y = lambda.mean)) +
  geom_errorbar(aes(ymin = lambda.5, 
                    ymax = lambda.95,
                    width=.01), position = position_dodge(0.1)) +
  geom_point(position = position_dodge(0.1))


# plot non-zero treaties
lambda.summary %>%
  filter(lambda.positive == 1 | lambda.negative == 1) %>% 
  ggplot(aes(x = atopid, y = lambda.mean)) +
  geom_errorbar(aes(ymin = lambda.5, 
                    ymax = lambda.95,
                    width=.01), position = position_dodge(0.1)) +
  geom_point(position = position_dodge(0.1))


# Load ATOP data for comparison
atop <- read.csv("data/atop-additions.csv")

# Join alliance coefficients with ATOP data
alliance.coefs <- left_join(atop, lambda.summary)


# Plot unconditional military support lambdas: lots of nulls
alliance.coefs %>%
  filter(uncond.milsup == 1) %>% 
  ggplot(aes(x = atopid, y = lambda.mean)) +
  geom_errorbar(aes(ymin = lambda.5, 
                    ymax = lambda.95,
                    width=.01), position = position_dodge(0.1)) +
  geom_point(position = position_dodge(0.1))

# Plot positive and negative, colored by unconditional military support
alliance.coefs %>%
filter(lambda.positive == 1 | lambda.negative == 1) %>% 
  ggplot(aes(x = atopid, y = lambda.mean, color = uncond.milsup)) +
  geom_errorbar(aes(ymin = lambda.5, 
                    ymax = lambda.95,
                    width=.01), position = position_dodge(0.1)) +
  geom_point(position = position_dodge(0.1))

# Plot lambdas against latent strength
ggplot(alliance.coefs, aes(y = lambda.mean, x = latent.str.mean)) + geom_point()

# non-negative Coefficients with error bars, colored by latent strength 
alliance.coefs %>%
  filter(lambda.positive == 1 | lambda.negative == 1) %>% 
  ggplot(aes(x = atopid, y = lambda.mean, color = latent.str.mean)) +
  geom_errorbar(aes(ymin = lambda.5, 
                    ymax = lambda.95,
                    width=.01), position = position_dodge(0.1)) +
  geom_point(position = position_dodge(0.1))
