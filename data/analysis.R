# Joshua Alley
# Texas A&M University
# Empirical test of public goods theory of alliances
# Examine Hypothesis 1 and 2


### Multilevel model- heterogenous effects across treaties
# Estimate a unique impact of changes in relative contribution for each alliance
# Requires estimating the model in STAN


# Load state-year matrix of alliance participation:
state.ally.year <- read.csv("data/alliance-state-year.csv")
state.vars <- read.csv("data/state-vars.csv")
atop.cow.year <- read.csv("data/atop-cow-year.csv")

# Merge state contributions to each alliance
# and keep only alliances that promise military support
atop.cow.year <- select(state.ally.year, atopid, ccode, year, contrib.gdp) %>%
  left_join(atop.cow.year) %>%
  group_by(atopid, ccode, year) %>%
  filter(defense == 1 | offense == 1) %>%
   mutate(bilat = ifelse(num.mem == 2, 1, 0))


# look at distribution of economic weight
summary(atop.cow.year$contrib.gdp)
ggplot(atop.cow.year, aes(x = contrib.gdp)) + geom_histogram()

# check missing by years
summary(subset(atop.cow.year, year < 1919, select = contrib.gdp)) # bilateral
summary(subset(atop.cow.year, year >= 1919, select = contrib.gdp)) # multilateral

# filter for post-45 only
atop.cow.year <- filter(atop.cow.year, year >= 1919)

# comparison bilateral and multilateral
summary(subset(atop.cow.year, bilat == 1, select = contrib.gdp)) # bilateral
summary(subset(atop.cow.year, bilat == 0, select = contrib.gdp)) # multilateral
ggplot(atop.cow.year, aes(x = as.factor(bilat), y = contrib.gdp)) + 
  geom_boxplot() + theme_bw()






# classify states as small or large in their alliances
# -1 if below 1st quartile, 1 if above in bilateral
# -1 if above 3rd quartile, 1 if above in multilateral
atop.cow.year$econ.size <- ifelse((atop.cow.year$contrib.gdp < 0.5000 & 
                                     atop.cow.year$bilat == 1) | # bilateral
                                  (atop.cow.year$contrib.gdp <= 0.08 & 
                                     atop.cow.year$bilat == 0), # multilateral
                                  -1, 0)
atop.cow.year$econ.size[(atop.cow.year$contrib.gdp >= 0.5000 & 
                             atop.cow.year$bilat == 1) | # bilateral
                          (atop.cow.year$contrib.gdp > 0.08 & 
                             atop.cow.year$bilat == 0)] <- 1
table(atop.cow.year$econ.size)


# weighted size: negative for small, positive for large
atop.cow.year$econ.size.w <- ifelse((atop.cow.year$contrib.gdp < 0.5 & 
                                     atop.cow.year$bilat == 1) | # bilateral
                                    (atop.cow.year$contrib.gdp <= 0.0707 & 
                                       atop.cow.year$bilat == 0), # multilateral
                                    (atop.cow.year$contrib.gdp - 1), atop.cow.year$contrib.gdp)
ggplot(atop.cow.year, aes(x = econ.size.w)) + geom_histogram()
summary(subset(atop.cow.year, bilat == 1, select = econ.size.w)) # bilateral
summary(subset(atop.cow.year, bilat == 0, select = econ.size.w)) # multilateral


# validity check: NATO coding
# US, UK and France all clear the 3rd quartile. 
nato <- filter(atop.cow.year, atopid == 3180) %>%
              select(contrib.gdp, econ.size, econ.size.w)


# Create a dataset of state-year alliance membership:
state.mem <- atop.cow.year %>% select(atopid, ccode, year, econ.size)
state.mem <- distinct(state.mem, atopid, ccode, year, .keep_all = TRUE)


# This matrix has each state's contribution to every alliance in a given year
# If a state is not a member of the alliance, corresponding matrix element = 0
state.mem <- spread(state.mem, key = atopid, value = econ.size, fill = 0)

# Flag alliances with some GDP data
apply(state.mem, 2, function(x) max(x, na.rm = TRUE) == 0 & 
                                min(x, na.rm = TRUE) == 0)

# remove alliances with no GDP data
state.mem <- state.mem[, apply(state.mem, 2, function(x) !(max(x, na.rm = TRUE) == 0 & 
                                 min(x, na.rm = TRUE) == 0))
                       ]

# Add state membership in alliances to this data
reg.state.data <- state.vars %>%
  select(ccode, year, growth.milex,
         atwar, civilwar.part, rival.milex, ln.gdp, polity, 
         cold.war, disputes, majpower) %>%
  filter(year >= 1919) %>%
  left_join(state.mem)

# fill in missing alliance data with zeros
reg.state.data[, 12:ncol(reg.state.data)][is.na(reg.state.data[, 12:ncol(reg.state.data)])] <- 0

# Create a matrix of state membership in alliances (Z in STAN model)
reg.state.data <- reg.state.data[complete.cases(reg.state.data), ]

state.mem.mat <- as.matrix(reg.state.data[, 12: ncol(reg.state.data)])
ncol(state.mem.mat)

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
                  X = reg.state.mat, M = ncol(reg.state.mat),
                  Z = state.mem.mat
)

# Compile the model code
model.1 <- stan_model(file = "data/ml-model-stan.stan")

# Variational Bayes- use to check model will run before going full STAN 
ml.model.vb <- vb(model.1, data = stan.data, seed = 12)
# Does not converge: diagnostics suggest problems

# Clear variational Bayes results from environment
rm(ml.model.vb)


# Run model with full Bayes
system.time(
  ml.model <- sampling(model.1, data = stan.data, 
                       iter = 2100, warmup = 1000, chains = 4
  )
)

# diagnose full model
# launch_shinystan(ml.model) 
check_hmc_diagnostics(ml.model)


# plot energy 
post.ml <- nuts_params(ml.model) # extract posterior draws
color_scheme_set("gray")
mcmc_nuts_energy(post.ml)
ggsave("appendix/energy-plot.png", height = 6, width = 8) 


# plot r-hats
ml.rhat <- rhat(ml.model)
mcmc_rhat_hist(ml.rhat)
ggsave("appendix/rhat-plot.png", height = 6, width = 8) 



# Extract coefficients from the model
ml.model.sum <- extract(ml.model, pars = c("beta", "gamma", "theta",
                                           "alpha",
                                           "sigma", "sigma_year", 
                                           "sigma_state", "sigma_all"),
                        permuted = TRUE)


# Summarize gamma
gamma.summary <- summary(ml.model, pars = c("gamma"), probs = c(0.05, 0.95))$summary
gamma.summary <- cbind.data.frame(as.numeric(colnames(state.mem.mat)), gamma.summary)
colnames(gamma.summary) <- c("atopid", "gamma.mean", "gamma.se.mean",
                             "gamma.sd", "gamma.5", "gamma.95",
                             "gamma.neff", "gamma.rhat")



# tabulate number of positive and negative estimates
# based on 90% credible intervals
gamma.summary$gamma.positive <- ifelse((gamma.summary$gamma.5 > 0 & gamma.summary$gamma.95 > 0), 1, 0)
sum(gamma.summary$gamma.positive) # 0 
gamma.summary$gamma.negative <- ifelse((gamma.summary$gamma.5 < 0 & gamma.summary$gamma.95 < 0), 1, 0)
sum(gamma.summary$gamma.negative) # 0


# Ignore uncertainty in estimates: are posterior means positive or negative? 
gamma.summary$positive.lmean <- ifelse(gamma.summary$gamma.mean > 0, 1, 0)
sum(gamma.summary$positive.lmean) # 13 treaties
gamma.summary$negative.lmean <- ifelse(gamma.summary$gamma.mean < 0, 1, 0)
sum(gamma.summary$negative.lmean) # 191 treaties


# Plot posterior means of alliance coefficients
ggplot(gamma.summary, aes(x = gamma.mean)) +
  geom_histogram(bins = 50) + theme_classic() +
  labs(x = "Posterior Mean", y = "Number of Alliances") +
  ggtitle("Distribution of Alliance Coefficient Posterior Means")



# Load ATOP data for comparison
atop <- read.csv("data/atop-additions.csv")

# Join alliance coefficients with ATOP data
alliance.coefs <- left_join(atop, gamma.summary) %>%
  filter(begyr >= 1899)


# Plot by start year of alliance
mean(ml.model.sum$theta)
ggplot(alliance.coefs, aes(x = begyr, y = gamma.mean)) +
  geom_errorbar(aes(ymin = gamma.5, 
                    ymax = gamma.95,
                    width=.01), position = position_dodge(0.01)) +
  geom_point(position = position_dodge(0.01)) + 
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = mean(ml.model.sum$theta), linetype = "dashed") +
  labs(x = "Start Year of Alliance", y = "Coefficient for Alliance Contribution") +
  theme_bw()
ggsave("manuscript/alliance-coefs-year.png", height = 6, width = 8)



# Calculate positive and negative posterior probability
positive.check <- function(x){
  mean(x > 0)
}
gamma.probs <- apply(ml.model.sum$gamma, 2, positive.check)
gamma.probs <- cbind.data.frame(gamma.probs, gamma.summary$atopid,
                                gamma.summary$gamma.mean)
colnames(gamma.probs) <- c("pos.post.prob", "atopid", "gamma.mean")


# Plot posterior probabilities
gamma.probs$atopid <- reorder(gamma.probs$atopid, gamma.probs$pos.post.prob)
gamma.probs$over.50 <- gamma.probs$pos.post.prob - .50

# For all alliances: plot relative posterior probability 
ggplot(gamma.probs, aes(x = atopid, y = over.50)) + 
  geom_col(color = "grey", fill = "black") +
  scale_y_continuous(breaks = seq(from = -.5, to = .5, .25),
                     labels = c("100% Negative", "75% negative", "Even", 
                                "75% Positive", "100% Positive")) +
  labs(x = "Alliance", y = "Posterior Probability") +
  theme_classic() +
  theme(axis.text.x = element_blank(), # remove atopid labels
        axis.ticks.x = element_blank()) +
  ggtitle("Posterior Probability of Alliance Coefficients")


# non-zero given 90% cutoff
gamma.probs$non.zero <- ifelse(gamma.probs$pos.post.prob >= .90 | 
                                gamma.probs$pos.post.prob <= .10, 1, 0)
sum(gamma.probs$non.zero)
# positive and negative
gamma.probs$nz.pos <- ifelse(gamma.probs$pos.post.prob >= .90 & 
                               gamma.probs$non.zero == 1, 1, 0)
sum(gamma.probs$nz.pos) # 0
gamma.probs$nz.neg <- ifelse(gamma.probs$pos.post.prob <= .10 & 
                               gamma.probs$non.zero == 1, 1, 0)
sum(gamma.probs$nz.neg) # 0


# Look at distribution of hyperparameters
# Variance hyperparameter
plot(density(ml.model.sum$sigma_all))
summary(ml.model.sum$sigma_all)
# mean hyperparameter
plot(density(ml.model.sum$theta))
summary(ml.model.sum$theta)


### estimate total effect of alliance participation on small and large members
total.all <- state.mem.mat %*% t(ml.model.sum$gamma)

total.all <- t(as.data.frame(
              apply(total.all, 1, function(x) quantile(x, 
                                    probs = c(.05, .95)))
))

total.all <- cbind.data.frame(reg.state.data$ccode,
                              reg.state.data$year,
                              total.all)

# gamma summary has the rough data needed to give estimates
glimpse(gamma.summary)

# predicted impact on large members
pred.imp.lg <- select(gamma.summary, atopid,
                      gamma.mean, gamma.5, gamma.95) %>%
                 mutate(mean.pred = sinh(gamma.mean), # undo transformation
                        pred.5 = sinh(gamma.5),
                        pred.95 = sinh(gamma.95))
glimpse(pred.imp.lg)  

# predicted impact on small members
pred.imp.sm <- select(gamma.summary, atopid,
                      gamma.mean, gamma.5, gamma.95) %>%
  mutate(mean.pred = -1 * sinh(gamma.mean), # undo transformation
         pred.5 = -1 * sinh(gamma.5), # and prediction w/ -1 in state.mem.mat
         pred.95 = -1 * sinh(gamma.95))
glimpse(pred.imp.sm)  

# bind rows
pred.imp.all <- bind_rows("large" = pred.imp.lg,
                          "small" = pred.imp.sm,
                          .id = "size")

# pull key alliance
pred.imp.key <- filter(pred.imp.all,
                  atopid == 3180 | # NATO
                  atopid == 3240 | # US-RoK
                  atopid == 3205 | # Arab League
                  atopid == 3150 | # OAS
                  atopid == 3285 # warsaw pact  
                   ) %>% 
                select(atopid, size,
                  mean.pred, pred.5, pred.95
                )
 
# plot 
ggplot(pred.imp.key, aes(y = factor(atopid), x = mean.pred,
                         group = factor(size),
                         color = factor(size)
                         )) +
  geom_point(position = position_dodge(0.3),
             size = 1.5) + # mean prediction
  geom_errorbarh( # error bars 
   aes(xmax = pred.95, xmin = pred.5),
   position = position_dodge(0.3),
   height = .1, size = 1.1
  ) +
  scale_color_manual(name="Economic Size", # manual colors
                       values = c("large" = "gray1",
                                  "small" = "gray45"),
                      breaks=c("large", "small"),
                      labels=c("Large", "Small")) +
  scale_y_discrete(name = "Alliance", # label alliances
                   breaks = c(3150, 3180, 3205, 3240, 3285),
                   labels = c("OAS", "NATO", "Arab League",
                   "US-Korea", "Warsaw Pact")) +
  labs(x = "Predicted Spending Growth") +
  theme_bw()
ggsave("manuscript/pred-sum.png", height = 6, width = 8)




### simulate impact of increasing share of allied GDP, given max positive gamma
# create relevant dataframe of coefficients
coef.sim <- cbind(ml.model.sum$alpha, ml.model.sum$beta, ml.model.sum$gamma[, 75])

# Create hypothetical dataset 
all.data.lshare <- numeric(ncol(reg.state.mat) + 2)
names(all.data.lshare) <- c("cons", colnames(reg.state.mat), "econ.share")

# summarize contrib/share of allied GDP
summary(atop.cow.year$contrib.gdp)

# Set values of variables for simulation 
all.data.lshare["cons"] <- 1 
all.data.lshare["atwar"] <- 0
all.data.lshare["civilwar.part"] <- 0 
all.data.lshare["rival.milex"] <- median(reg.state.data$rival.milex)
all.data.lshare["ln.gdp"] <- median(reg.state.data$ln.gdp)
all.data.lshare["polity"] <-median(reg.state.data$polity)
all.data.lshare["cold.war"] <- 0
all.data.lshare["disputes"] <- 0
all.data.lshare["majpower"] <- 0
all.data.lshare["econ.share"] <- -1 # set share to minimum

# Matrix multiplication for prediction: low 
pred.lshare <- coef.sim%*%all.data.lshare



# move share to highest value = 1 
all.data.hshare <- all.data.lshare
all.data.hshare["econ.share"] <- 1

# Matrix multiplication for prediction: high 
pred.hshare <- coef.sim%*%all.data.hshare


# Difference in predicted values
pred.diff <- pred.hshare - pred.lshare


# Compare predictions"
pred.data <- cbind.data.frame(pred.lshare, pred.hshare, pred.diff)
colnames(pred.data) <- c("Low Share", "High Share", "Difference")

# plot 90% credible intervals 
color_scheme_set("gray")
mcmc_intervals(pred.data, prob = .9) +
  labs(x = "Predicted Percentage Change in Military Spending", y = "Share of Allied GDP") +
  theme_bw()
ggsave("appendix/pred-change-share.png", height = 6, width = 8)


# Remove fit model from workspace
saveRDS(ml.model, "data/ml-model.rds")
rm(ml.model)
