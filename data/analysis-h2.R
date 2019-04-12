# Joshua Alley
# Texas A&M University
# Empirical test of public goods theory of alliances
# Examine Hypothesis 2 



# Load packages
library(tidyverse)
library(rstan)
library(bayesplot)
library(shinystan)


# Set-up STAN guidelines
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


# Set seed 
set.seed(12)



### Second approach: Multilevel model- heterogenous effects across treaties
# Estimate a unique impact of changes in relative contribution for each alliance
# Requires estimating the model in STAN


# Load state-year matrix of alliance participation:
atop.cow.year <- read.csv("data/atop-cow-year.csv")
# Merge state contributions to each alliance
# and keep only alliances that promise military support
atop.cow.year <- select(state.ally.year, atopid, ccode, year, alliance.contrib) %>%
  left_join(atop.cow.year) %>%
  group_by(atopid, ccode, year) %>%
  filter(defense == 1 | offense == 1)

# Create a dataset of state-year alliance membership:
state.mem <- atop.cow.year %>% select(atopid, ccode, year, alliance.contrib)
state.mem <- distinct(state.mem, atopid, ccode, year, .keep_all = TRUE)


# This matrix has each state's contribution to every alliance in a given year
# If a state is not a member of the alliance, corresponding matrix element = 0
state.mem <- spread(state.mem, key = atopid, value = alliance.contrib, fill = 0)

# Remove the zero or no alliance category
state.mem <- subset(state.mem, select = -(3))


# Add state membership in alliances to this data
reg.state.data <-  state.vars %>%
  select(ccode, year, growth.milex,
         atwar, civilwar.part, rival.milex, ln.gdp, polity, 
         cold.war, disputes, majpower) %>%
  left_join(state.mem)

# fill in missing alliance data with zeros
reg.state.data[, 12:ncol(reg.state.data)][is.na(reg.state.data[, 12:ncol(reg.state.data)])] <- 0

# Create a matrix of state membership in alliances (Z in STAN model)
reg.state.data <- reg.state.data[complete.cases(reg.state.data), ]
state.mem.mat <- as.matrix(reg.state.data[, 12: ncol(reg.state.data)])

# Rescale state regression variables
reg.state.data[, 5:11] <- lapply(reg.state.data[, 5:11], 
                                 function(x) rescale(x, binary.inputs = "0/1"))

reg.state.mat <- as.matrix(reg.state.data[, 4:11])

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

# Variational Bayes- use to check model will run and give reasonable predictions
ml.model.vb <- vb(model.1, data = stan.data, seed = 12)
# Does not converge

# posterior predictive check from variational Bayes- did not converge
# so treat these predictions with caution
y = reg.state.data[, 3]
vb.model.sum <- extract(ml.model.vb)
ppc_dens_overlay(y, vb.model.sum$y_pred[1:100, ])

# Clear variational Bayes results from environment
rm(ml.model.vb)
rm(vb.model.sum)


# Run model with full Bayes
system.time(
  ml.model <- sampling(model.1, data = stan.data, 
                       iter = 2000, warmup = 1000, chains = 4
  )
)

# diagnose full model
launch_shinystan(ml.model)

check_hmc_diagnostics(ml.model)


# plot energy 
post.ml <- nuts_params(ml.model) # extract posterior draws
color_scheme_set("red")
mcmc_nuts_energy(post.ml)
ggsave("appendix/energy-plot.pdf", height = 6, width = 8) 


# plot r-hats
ml.rhat <- rhat(ml.model)
mcmc_rhat_hist(ml.rhat)
ggsave("appendix/rhat-plot.pdf", height = 6, width = 8) 



# Extract coefficients from the model
ml.model.sum <- extract(ml.model, pars = c("beta", "gamma", 
                                           "sigma", "sigma_year", "sigma_state",
                                           "y_pred"),permuted = TRUE)

# Posterior predictive distributions relative to observed data
yrep <- ml.model.sum$y_pred[1:100, ]

# plot posterior predictive denisty of first 100 simulations
ppc_dens_overlay(y, yrep)


# Summarize gamma
gamma.summary <- summary(ml.model, pars = c("gamma"), probs = c(0.05, 0.95))$summary
gamma.summary <- cbind.data.frame(as.numeric(colnames(state.mem.mat)), gamma.summary)
colnames(gamma.summary) <- c("atopid", "gamma.mean", "gamma.se.mean",
                             "gamma.sd", "gamma.5", "gamma.95",
                             "gamma.neff", "gamma.rhat")

# tabulate number of positive and negative estimates
gamma.summary$gamma.positive <- ifelse((gamma.summary$gamma.5 > 0 & gamma.summary$gamma.95 > 0), 1, 0)
sum(gamma.summary$gamma.positive) # 16 treaties: increasing contribution to alliance leads to increased spending
gamma.summary$gamma.negative <- ifelse((gamma.summary$gamma.5 < 0 & gamma.summary$gamma.95 < 0), 1, 0)
sum(gamma.summary$gamma.negative) # 13 treaties: increasing contribution to alliance leads to decreased spending


# Ignore uncertainty in estimates: are posterior means positive or negative? 
gamma.summary$positive.lmean <- ifelse(gamma.summary$gamma.mean > 0, 1, 0)
sum(gamma.summary$positive.lmean) # 147 treaties
gamma.summary$negative.lmean <- ifelse(gamma.summary$gamma.mean < 0, 1, 0)
sum(gamma.summary$negative.lmean) # 138 treaties

# Plot posterior means of alliance coefficients
ggplot(gamma.summary, aes(x = gamma.mean)) +
  geom_density() + theme_classic() +
  ggtitle("Posterior Means of Alliance Coefficients")

ggplot(gamma.summary, aes(x = gamma.mean)) +
  geom_histogram(bins = 50) + theme_classic() +
  labs(x = "Posterior Mean") +
  ggtitle("Distribution of Alliance Coefficient Posterior Means")
ggsave("manuscript/alliance-coefs-hist.pdf", height = 6, width = 8)


# Plot points with error bars by ATOPID
ggplot(gamma.summary, aes(x = atopid, y = gamma.mean)) +
  geom_errorbar(aes(ymin = gamma.5, 
                    ymax = gamma.95,
                    width=.01), position = position_dodge(0.1)) +
  geom_point(position = position_dodge(0.1)) + geom_hline(yintercept = 0) +
  theme_classic()


# plot non-zero treaties
gamma.summary %>%
  filter(gamma.positive == 1 | gamma.negative == 1) %>% 
  ggplot(aes(x = atopid, y = gamma.mean)) +
  geom_errorbar(aes(ymin = gamma.5, 
                    ymax = gamma.95,
                    width=.01), position = position_dodge(0.1)) +
  geom_point(position = position_dodge(0.1)) + geom_hline(yintercept = 0) + 
  theme_classic()


# Load ATOP data for comparison
atop <- read.csv("data/atop-additions.csv")

# Join alliance coefficients with ATOP data
alliance.coefs <- left_join(atop, gamma.summary)


# Plot by start year of alliance
ggplot(alliance.coefs, aes(x = begyr, y = gamma.mean)) +
  geom_errorbar(aes(ymin = gamma.5, 
                    ymax = gamma.95,
                    width=.01), position = position_dodge(0.01)) +
  geom_point(position = position_dodge(0.01)) + geom_hline(yintercept = 0) +
  labs(x = "Start Year of Alliance", y = "Coefficient for Alliance Contribution") +
  theme_classic() 
ggsave("manuscript/alliance-coefs-year.pdf", height = 6, width = 8)


# Positive and negative only
alliance.coefs %>%
  filter(gamma.positive == 1 | gamma.negative == 1) %>% 
  ggplot(aes(x = begyr, y = gamma.mean)) +
  geom_errorbar(aes(ymin = gamma.5, 
                    ymax = gamma.95,
                    width=.01), position = position_dodge(0.1)) +
  geom_point(position = position_dodge(0.1)) + geom_hline(yintercept = 0) +
  labs(x = "Start Year of Alliance", y = "Coefficient for Alliance Contribution") +
  theme_classic() 


# 15 / 272 defense pacts have a positive association between contribution and changes in spending
table(alliance.coefs$gamma.positive, alliance.coefs$defense)
# 12 / 272 defense pacts have a negative association beween contribution and changes in spending 
table(alliance.coefs$gamma.negative, alliance.coefs$defense)


# For non-zero alliances 
alliance.coefs$atopid <- reorder(alliance.coefs$atopid, alliance.coefs$gamma.mean)
alliance.coefs %>%
  filter(gamma.positive == 1 | gamma.negative == 1) %>% 
  ggplot(mapping = aes(x = atopid, y = gamma.mean)) + 
  geom_col() +
  scale_fill_brewer(palette = "Greys") +
  #  geom_text(aes(label = round(gamma.mean, digits = 3)), nudge_y = 0.075, size = 4) +
  labs(x = "ATOPid", y = "Posterior Mean of Alliance Parameter") +
  coord_flip() + theme_classic() 
ggsave("manuscript/nonzero-alliance-coefs.pdf", height = 6, width = 8)


# Plot gammas against latent strength
ggplot(alliance.coefs, aes(y = gamma.mean, x = latent.str.mean)) + 
  geom_point()  + theme_classic()

# non-negative Coefficients with error bars, colored by latent scope 
alliance.coefs %>%
  filter(gamma.positive == 1 | gamma.negative == 1) %>% 
  ggplot(aes(x = atopid, y = gamma.mean, color = latent.str.mean)) +
  geom_errorbar(aes(ymin = gamma.5, 
                    ymax = gamma.95,
                    width=.01), position = position_dodge(0.1)) +
  geom_point(position = position_dodge(0.1)) + geom_hline(yintercept = 0) +
  theme_classic()



# non-negative coefficients with error bars, colored by offense
alliance.coefs %>%
  filter(gamma.positive == 1 | gamma.negative == 1) %>% 
  ggplot(aes(x = atopid, y = gamma.mean, color = factor(offense))) +
  geom_errorbar(aes(ymin = gamma.5, 
                    ymax = gamma.95,
                    width=.01), position = position_dodge(0.1)) +
  geom_point(position = position_dodge(0.1)) + geom_hline(yintercept = 0) +
  theme_classic()




### 
# Robustness check: estimate STAN model on states only in alliances
# Filter out obs where states are not in at least one alliance
reg.data.all <- reg.state.data %>%
  select(-c(state.id, year.id)) %>% 
  mutate(allied.cap = rowSums(.[12: ncol(reg.data.all)])) %>% 
  filter(allied.cap != 0)

state.mem.all <- as.matrix(reg.data.all[, 12: (ncol(reg.data.all) - 1)])

# Define regression data 
reg.mat.all <- as.matrix(reg.data.all[, 4:11])

# Add state and year indicators
reg.data.all$state.id <- group_indices(reg.data.all, ccode)
reg.data.all$year.id <- group_indices(reg.data.all, year)

# Define the data list 
stan.data.all <- list(N = nrow(reg.data.all), y = reg.data.all[, 3],
                      state = reg.data.all$state.id, S = length(unique(reg.data.all$state.id)),
                      year = reg.data.all$year.id, T = length(unique(reg.data.all$year.id)),
                      A = ncol(state.mem.mat),
                      Z = state.mem.all, 
                      X = reg.mat.all, M = ncol(reg.mat.all)
)


# Run model with full Bayes
system.time(
  ml.model.all <- sampling(model.1, data = stan.data.all, 
                           iter = 2000, warmup = 1000, chains = 4
  )
)

# diagnose full model
check_hmc_diagnostics(ml.model.all)



# Extract coefficients from the model
sum.ml.all <- extract(ml.model.all, pars = c("gamma"), permuted = TRUE)


# Summarize gamma
gamma.summary.all <- summary(ml.model.all, pars = c("gamma"), probs = c(0.05, 0.95))$summary
gamma.summary.all <- cbind.data.frame(as.numeric(colnames(state.mem.all)), gamma.summary.all)
colnames(gamma.summary.all) <- c("atopid", "gamma.mean", "gamma.se.mean",
                                 "gamma.sd", "gamma.5", "gamma.95",
                                 "gamma.neff", "gamma.rhat")

# tabulate number of positive and negative estimates
gamma.summary.all$gamma.positive <- ifelse((gamma.summary.all$gamma.5 > 0 & gamma.summary.all$gamma.95 > 0), 1, 0)
sum(gamma.summary.all$gamma.positive) # 18 treaties: increasing contribution to alliance leads to increased spending
gamma.summary.all$gamma.negative <- ifelse((gamma.summary.all$gamma.5 < 0 & gamma.summary.all$gamma.95 < 0), 1, 0)
sum(gamma.summary.all$gamma.negative) # 15 treaties: increasing contribution to alliance leads to decreased spending

# Prediction success 
18 / 285

# Ignore uncertainty in estimates: are posterior means positive or negative? 
gamma.summary.all$positive.lmean <- ifelse(gamma.summary.all$gamma.mean > 0, 1, 0)
sum(gamma.summary.all$positive.lmean) # 138 treaties
gamma.summary.all$negative.lmean <- ifelse(gamma.summary.all$gamma.mean < 0, 1, 0)
sum(gamma.summary.all$negative.lmean) # 147 treaties

# Plot posterior means of alliance coefficients
ggplot(gamma.summary.all, aes(x = gamma.mean)) +
  geom_density() + theme_classic() +
  ggtitle("Posterior Means of Alliance Coefficients")

ggplot(gamma.summary.all, aes(x = gamma.mean)) +
  geom_histogram(bins = 50) + theme_classic() +
  labs(x = "Posterior Mean") +
  ggtitle("Distribution of Alliance Coefficient Posterior Means")
ggsave("appendix/all-sample-gamma")


# Create a data frame with gamma pars from both samples 
gamma.comp <- cbind.data.frame(gamma.summary$gamma.mean, gamma.summary.all$gamma.mean)
colnames(gamma.comp) <- c("full", "allies")
gamma.comp <- melt(gamma.comp)

# plot 
ggplot(gamma.comp, aes(x = value, fill = variable)) + 
  geom_density(alpha = 0.25) +
  scale_fill_manual(name = "Sample", values=c("#999999", "#000000")) +
  ggtitle("Posterior Distributions of Treaty Scope: Major and Non-Major Powers") +
  theme_classic()
ggsave("appendix/sample-comp-gamma.pdf")



