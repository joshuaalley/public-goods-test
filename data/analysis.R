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
library(stargazer)


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
    total.ally.expend = sum(ally.spend[defense == 1 | offense == 1], na.rm = TRUE),
    avg.treaty.contrib = mean(alliance.contrib[defense == 1 | offense == 1], na.rm = TRUE),
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
state.char.full[, 42: ncol(state.char.full)][is.na(state.char.full[, 42: ncol(state.char.full)])] <- 0




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
summary(state.char.inter$change.gdp)
ggplot(state.char.inter, aes(x = change.gdp)) + geom_density()


summary(state.char.inter$total.ally.expend)
ggplot(state.char.inter, aes(x = total.ally.expend)) + geom_density()


summary(state.char.inter$ln.ally.expend)
ggplot(state.char.inter, aes(x = ln.ally.expend)) + geom_histogram()



# tests for panel unit roots
y.panel <- state.char.inter %>% select(ccode, ln.milex, year) %>%
           mutate(ln.milex = approx(year, ln.milex, year)$y) %>% # interpolate missing data points
            spread(ccode, ln.milex) %>% select(-(year))
y.panel <- as.matrix(y.panel)
# formal test is not working at the moment- reasons are unclear but some points surrounded by missing values after interpolation
# choi.test <- pCADFtest(y.panel, max.lag.y = 4, crosscorr = .05, type = "trend")

# plot on ln shows clear evidence of non-stationarity
ggplot(state.char.inter, aes(x = year, y = ln.milex, colour = as.factor(ccode))) + 
  geom_line()

# Growth is more mean-reverting
ggplot(state.char.inter, aes(x = year, y = asinh(growth.milex), colour = as.factor(ccode))) + 
  geom_line()



### First test: aboslute size (GDP)
# Interact changes in allied spending and GDP
# Total allied spending: pooling regression
m1.pg.abs <- rlm(growth.milex ~ diff.ally.expend + ln.gdp + diff.ally.expend:ln.gdp +
                     avg.num.mem + avg.dem.prop + 
                     atwar + civilwar.part + polity  +
                     lsthreat + cold.war,
                   data = state.char.inter
                )
summary(m1.pg.abs)
plotreg(m1.pg.abs)
stargazer(m1.pg.abs)

# Calculate marginal effects
margins(m1.pg.abs)
cplot(m1.pg.abs, x = "ln.gdp", dx = "diff.ally.expend", what = "effect",
      main = "Marginal Effect of Changes in Allied Spending on Growth in Military Spending",
      xlab = "ln(GDP)", ylab = "Average M.E. of Changes in Allied Spending")
abline(h = 0)

# Switch the directions 
cplot(m1.pg.abs, x = "diff.ally.expend", dx = "ln.gdp", what = "effect",
      main = "Marginal Effect GDP on Growth in Military Spending",
      xlab = "Change in Allied Capability", ylab = "Average M.E. of ln(GDP)")
abline(h = 0)



# OLS
m2.pg.abs <- lm(growth.milex ~ diff.ally.expend + ln.gdp + diff.ally.expend:ln.gdp +
                       avg.num.mem + avg.dem.prop + 
                       atwar + civilwar.part + polity + ln.gdp + 
                       lsthreat + cold.war,
                     data = state.char.inter)
summary(m2.pg.abs)



# binning estimator
bin.abs <- inter.binning(Y = "growth.milex", D = "diff.ally.expend", X = "ln.gdp", 
              Z = c("lag.ln.milex", "atwar", "civilwar.part", "polity", 
                    "lsthreat", "cold.war", "avg.num.mem", "avg.dem.prop"), 
              data = state.char.inter,
              na.rm = TRUE,
              Ylabel = "Growth in Military Spending",
              Dlabel = "Allied Spending",
              Xlabel = "ln(GDP)", theme.bw = TRUE
)
bin.abs
ggsave("appendix/inter-bin-abs.pdf", height = 6, width = 8)


# Kernel: 10+ minute run time 
kernel.abs <- inter.kernel(Y = "growth.milex", D = "diff.ally.expend", X = "ln.gdp", 
             Z = c("lag.ln.milex", "atwar", "civilwar.part", "polity", 
                   "lsthreat", "cold.war", "avg.num.mem", "avg.dem.prop"), 
             data = state.char.inter, 
             na.rm = TRUE,
             nboots = 200, parallel = TRUE, cores = 4,
             Ylabel = "Growth in Military Spending",
             Dlabel = "Allied Spending",
             Xlabel = "ln(GDP)", theme.bw = TRUE
)
kernel.abs
ggsave("appendix/inter-kernel-abs.pdf", height = 6, width = 8)


# Check for non-random selection into alliances and compare allied states
# Still a statistically significant interaction
heckit.ally.spend <- heckit(selection = treaty.pres ~ lag.ln.milex + 
                              ln.gdp + polity + atwar + lsthreat,
                            outcome = growth.milex ~ diff.ally.expend + ln.gdp + 
                              diff.ally.expend:ln.gdp +
                              atwar + civilwar.part + polity + ln.gdp +
                              lsthreat + cold.war,
                            data = state.char.inter,
                            method = "ml")
summary(heckit.ally.spend)

# Create a table for the appendix
stargazer(heckit.ally.spend)



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






### Transform outcome with Inverse hyperbolic sine to reduce the extremely heavy tails 
# Total allied spending: pooling regression
m1.abs.ihs <- rlm(asinh(growth.milex) ~ diff.ally.expend + ln.gdp + diff.ally.expend:ln.gdp +
                   avg.num.mem + avg.dem.prop + 
                   atwar + civilwar.part + polity  +
                   lsthreat + cold.war,
                 data = state.char.inter
)
summary(m1.abs.ihs)
plotreg(m1.abs.ihs)
stargazer(m1.abs.ihs)

# Calculate marginal effects
margins(m1.abs.ihs)
cplot(m1.abs.ihs, x = "ln.gdp", dx = "diff.ally.expend", what = "effect",
      main = "Marginal Effect of Changes in Allied Spending on Growth in Military Spending",
      xlab = "ln(GDP)", ylab = "Average M.E. of Changes in Allied Spending")
abline(h = 0)

# Switch the directions 
cplot(m1.abs.ihs, x = "diff.ally.expend", dx = "ln.gdp", what = "effect",
      main = "Marginal Effect GDP on Growth in Military Spending",
      xlab = "Change in Allied Capability", ylab = "Average M.E. of ln(GDP)")
abline(h = 0)


# Results are very similar- check distribution of weights across observations
plot(m1.pg.abs$residuals, m1.pg.abs$w) # original model
abline(h = 0)
plot(m1.abs.ihs$residuals, m1.abs.ihs$w) # transformed DV


# OLS estiamtion 
m2.abs.ihs <- lm(asinh(growth.milex) ~ diff.ally.expend + ln.gdp + diff.ally.expend:ln.gdp +
                    avg.num.mem + avg.dem.prop + 
                    atwar + civilwar.part + polity  +
                    lsthreat + cold.war,
                  data = state.char.inter
)
summary(m2.abs.ihs)

cplot(m2.abs.ihs, x = "ln.gdp", dx = "diff.ally.expend", what = "effect",
      main = "Marginal Effect of Changes in Allied Spending on Growth in Military Spending",
      xlab = "ln(GDP)", ylab = "Average M.E. of Changes in Allied Spending")
abline(h = 0)


# Take this and ols estiamtes in one appendix table
stargazer(m2.pg.abs, m2.abs.ihs, m1.abs.ihs)

# plot residuals from transformed regression 
qqnorm(m2.abs.ihs$residuals, main = "Normal Q-Q Plot: Regression Residuals")
qqline(m2.abs.ihs$residuals)
# Export plot
dev.copy(pdf,'appendix/res-qq-plot.pdf')
dev.off()


# Marginal effects plots all in one place
par(mfrow = c(2, 2))

# Main model in paper: robust regression 
cplot(m1.pg.abs, x = "ln.gdp", dx = "diff.ally.expend", what = "effect",
      main = "Robust Reg",
      xlab = "ln(GDP)", ylab = "Average M.E. of Allied Cap")

# rreg- transformed DV
cplot(m1.abs.ihs, x = "ln.gdp", dx = "diff.ally.expend", what = "effect",
      main = "Robust Reg: Transformed DV",
      xlab = "ln(GDP)", ylab = "Average M.E. of Allied Cap")

# OLS 
cplot(m2.pg.abs, x = "ln.gdp", dx = "diff.ally.expend", what = "effect",
      main = "Linear Regression",
      xlab = "ln(GDP)", ylab = "Average M.E. of Allied Cap")

# OLS with transformed DV 
cplot(m2.abs.ihs, x = "ln.gdp", dx = "diff.ally.expend", what = "effect",
      main = "Linear Regression: Transformed DV",
      xlab = "ln(GDP)", ylab = "Average M.E. of Allied Cap")
# Export plot
dev.copy(pdf,'appendix/me-plots.pdf')
dev.off()
dev.off()


### Additional single-level test: relative size expressed as contribution to alliance
# estimate interactions
# filter out cases with no alliances
inter.data.rel <- filter(state.char.full, treaty.pres == 1)
inter.data.rel <- as.data.frame(inter.data.rel)

# Total allied spending: pooling regression
m1.pg.rel <- rlm(growth.milex ~ diff.ally.expend + avg.treaty.contrib + diff.ally.expend:avg.treaty.contrib +
                   lag.ln.milex + avg.num.mem + avg.dem.prop + 
                   atwar + civilwar.part + polity  + 
                   lsthreat + cold.war,
                 data = inter.data.rel)
summary(m1.pg.rel)
# Calculate marginal effects
margins(m1.pg.rel)
cplot(m1.pg.rel, x = "avg.treaty.contrib", dx = "diff.ally.expend", what = "effect",
      main = "Marginal Effect of Changes in Allied Spending on Military Spending",
      xlab = "Average Alliance Contribution", ylab = "Average M.E. of Changes in Allied Spending")
abline(h = 0)

# OLS
m2.pg.rel <- lm(growth.milex ~ diff.ally.expend + avg.treaty.contrib + diff.ally.expend:avg.treaty.contrib +
                     avg.dem.prop + lag.ln.milex +
                     atwar + civilwar.part + polity + ln.gdp + avg.num.mem +
                     lsthreat + cold.war,
                   data = inter.data.rel
                )
summary(m2.pg.rel)


# binning estimator
bin.rel <- inter.binning(Y = "growth.milex", D = "diff.ally.expend", X = "avg.treaty.contrib", 
                         Z = c("lag.ln.milex", "atwar", "civilwar.part", "polity", 
                               "lsthreat", "cold.war", "avg.num.mem", "avg.dem.prop"), 
                         data = inter.data.rel,
                         na.rm = TRUE
)
bin.rel 

# Kernel: 10+ minute run time 
kernel.rel <- inter.kernel(Y = "growth.milex", D = "diff.ally.expend", X = "avg.treaty.contrib", 
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
                     outcome = growth.milex ~ diff.ally.expend + avg.treaty.contrib +
                       diff.ally.expend:avg.treaty.contrib +
                       lag.ln.milex +
                       atwar + civilwar.part + polity + ln.gdp +
                       lsthreat + cold.war,
                     data = state.char.full,
                     method = "ml")
summary(heckit.rel)

