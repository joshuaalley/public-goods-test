# Joshua Alley
# Texas A&M University
# Analysis with weighted alliance


# load key package
library(hexbin)

# Uses positive and negative GDP share values from the analysis-h1 script

# Create a dataset of state-year alliance membership:
state.mem.w <- atop.cow.year %>% select(atopid, ccode, year, econ.size.w)
state.mem.w <- distinct(state.mem.w, atopid, ccode, year, .keep_all = TRUE)


# This matrix has each state's contribution to every alliance in a given year
# If a state is not a member of the alliance, corresponding matrix element = 0
state.mem.w <- spread(state.mem.w, key = atopid, value = econ.size.w, fill = 0)


# Flag alliances with some GDP data
apply(state.mem.w, 2, function(x) max(x, na.rm = TRUE) == 0 & 
        min(x, na.rm = TRUE) == 0)

# remove alliances with no GDP data
state.mem.w <- state.mem.w[, apply(state.mem.w, 2, function(x) !(max(x, na.rm = TRUE) == 0 & 
                                                             min(x, na.rm = TRUE) == 0))
                       ]


# Add state membership in alliances to this data
reg.state.data.w <- state.vars %>%
  select(ccode, year, growth.milex,
         atwar, civilwar.part, rival.milex, ln.gdp, polity, 
         cold.war, disputes, majpower) %>%
  filter(year >= 1919) %>%
  left_join(state.mem.w)

# fill in missing alliance data with zeros
reg.state.data.w[, 12:ncol(reg.state.data.w)][is.na(reg.state.data.w[, 12:ncol(reg.state.data.w)])] <- 0

# Create a matrix of state membership in alliances (Z in STAN model)
reg.state.data.w <- reg.state.data.w[complete.cases(reg.state.data.w), ]

state.mem.matw <- as.matrix(reg.state.data.w[, 12: ncol(reg.state.data.w)])


# Rescale state regression variables
reg.state.data.w[, 5:11] <- lapply(reg.state.data.w[, 5:11], 
                                 function(x) rescale(x, binary.inputs = "0/1"))

reg.state.matw <- as.matrix(reg.state.data.w[, 4:11])

# Set-up data for STAN
# create a state index variable
reg.state.data.w$state.id <- reg.state.data.w %>% group_indices(ccode)
# Create a year index variable 
reg.state.data.w$year.id <- reg.state.data.w %>% group_indices(year)



# Define the data list 
stan.data.w <- list(N = nrow(reg.state.data.w), y = reg.state.data.w[, 3],
                  state = reg.state.data.w$state.id, S = length(unique(reg.state.data.w$state.id)),
                  year = reg.state.data.w$year.id, T = length(unique(reg.state.data.w$year.id)),
                  A = ncol(state.mem.matw),
                  X = reg.state.matw, M = ncol(reg.state.matw),
                  Z = state.mem.matw
)

# Compile the model code
model.1 <- stan_model(file = "data/ml-model-stan.stan")

# Run model with full Bayes
system.time(
  ml.model.w <- sampling(model.1, data = stan.data.w, 
                       iter = 2100, warmup = 1000, chains = 4
  )
)

# diagnose full model
# launch_shinystan(ml.model) 
check_hmc_diagnostics(ml.model.w)




# Extract coefficients from the model
model.sum.w <- extract(ml.model.w, pars = c("gamma", "theta",
                                            "sigma_all"),
                        permuted = TRUE)


# Summarize gamma
gamma.summary.w <- summary(ml.model.w, pars = c("gamma"), probs = c(0.05, 0.95))$summary
gamma.summary.w <- cbind.data.frame(as.numeric(colnames(state.mem.matw)), gamma.summary.w)
colnames(gamma.summary.w) <- c("atopid", "gamma.mean", "gamma.se.mean",
                             "gamma.sd", "gamma.5", "gamma.95",
                             "gamma.neff", "gamma.rhat")



# tabulate number of positive and negative estimates
# based on 90% credible intervals
gamma.summary.w$gamma.positive <- ifelse((gamma.summary.w$gamma.5 > 0 & gamma.summary.w$gamma.95 > 0), 1, 0)
sum(gamma.summary.w$gamma.positive) 
gamma.summary.w$gamma.negative <- ifelse((gamma.summary.w$gamma.5 < 0 & gamma.summary.w$gamma.95 < 0), 1, 0)
sum(gamma.summary.w$gamma.negative) 

# Ignore uncertainty in estimates: are posterior means positive or negative? 
gamma.summary.w$positive.lmean <- ifelse(gamma.summary.w$gamma.mean > 0, 1, 0)
sum(gamma.summary.w$positive.lmean)
gamma.summary.w$negative.lmean <- ifelse(gamma.summary.w$gamma.mean < 0, 1, 0)
sum(gamma.summary.w$negative.lmean) 


# Plot posterior means of alliance coefficients
ggplot(gamma.summary.w, aes(x = gamma.mean)) +
  geom_histogram(bins = 50) + theme_classic() +
  labs(x = "Posterior Mean", y = "Number of Alliances") +
  ggtitle("Distribution of Alliance Coefficient Posterior Means: Weighted Z")


# Plot points with error bars by ATOPID
ggplot(gamma.summary.w, aes(x = atopid, y = gamma.mean)) +
  geom_errorbar(aes(ymin = gamma.5, 
                    ymax = gamma.95,
                    width=.01), position = position_dodge(0.1)) +
  geom_point(position = position_dodge(0.1)) + geom_hline(yintercept = 0) +
  theme_classic() + coord_flip()



# Load ATOP data for comparison
atop <- read.csv("data/atop-additions.csv")

# Join alliance coefficients with ATOP data
alliance.coefs.w <- left_join(atop, gamma.summary.w) %>%
  filter(begyr >= 1898)


# Plot by start year of alliance
mean(model.sum.w$theta)
ggplot(alliance.coefs.w, aes(x = begyr, y = gamma.mean)) +
  geom_errorbar(aes(ymin = gamma.5, 
                    ymax = gamma.95,
                    width=.01), position = position_dodge(0.01)) +
  geom_point(position = position_dodge(0.01)) + 
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = mean(ml.model.sum$theta), linetype = "dashed") +
  labs(x = "Start Year of Alliance", y = "Coefficient for Alliance Contribution") +
  theme_classic()
ggsave("appendix/alliance-coefs-year-w.pdf", height = 6, width = 8)



# Calculate positive and negative posterior probability
positive.check <- function(x){
  mean(x > 0)
}
gamma.probs.w <- apply(model.sum.w$gamma, 2, positive.check)
gamma.probs.w <- cbind.data.frame(gamma.probs.w, gamma.summary$atopid,
                                gamma.summary$gamma.mean)
colnames(gamma.probs.w) <- c("pos.post.prob", "atopid", "gamma.mean")


# Plot posterior probabilities
gamma.probs.w$atopid <- reorder(gamma.probs.w$atopid, gamma.probs.w$pos.post.prob)
gamma.probs.w$over.50 <- gamma.probs.w$pos.post.prob - .50

# For all alliances: plot relative posterior probability 
ggplot(gamma.probs.w, aes(x = atopid, y = over.50)) + 
  geom_col(color = "grey", fill = "black") +
  scale_y_continuous(breaks = seq(from = -.5, to = .5, .25),
                     labels = c("100% Negative", "75% negative", "Even", 
                                "75% Positive", "100% Positive")) +
  labs(x = "Alliance", y = "Posterior Probability") +
  theme_classic() +
  theme(axis.text.x = element_blank(), # remove atopid labels
        axis.ticks.x = element_blank()) +
  ggtitle("Posterior Probability of Alliance Coefficients: Weighted Z Values")


# non-zero given 90% cutoff
gamma.probs.w$non.zero <- ifelse(gamma.probs.w$pos.post.prob >= .90 | 
                                 gamma.probs.w$pos.post.prob <= .10, 1, 0)
sum(gamma.probs.w$non.zero)
# positive and negative
gamma.probs.w$nz.pos <- ifelse(gamma.probs.w$pos.post.prob >= .90 & 
                               gamma.probs.w$non.zero == 1, 1, 0)
sum(gamma.probs.w$nz.pos) 
gamma.probs.w$nz.neg <- ifelse(gamma.probs.w$pos.post.prob <= .10 & 
                               gamma.probs.w$non.zero == 1, 1, 0)
sum(gamma.probs.w$nz.neg) 


# Look at distribution of hyperparameters
# Variance hyperparameter
plot(density(model.sum.w$sigma_all))
summary(model.sum.w$sigma_all)
# mean hyperparameter
plot(density(model.sum.w$theta))
summary(model.sum.w$theta)




# Predicted military spending change for all individual alliances
# plot against economic weight
a <- ncol(state.mem.matw)
growth.pred <- rep(NA, a)
growth.pred <- list()

# Loop over matrix columns
for(i in 1:a){
  growth.pred[[i]] <- state.mem.matw[, i][state.mem.matw[, i] != 0] # filter out zeros 
  growth.pred[[i]] <- growth.pred[[i]]%*%t(model.sum.w$gamma[, i]) # multiply by gamma
  growth.pred[[i]] <- as.data.frame(growth.pred[[i]])
}

names(growth.pred) <- c(colnames(state.mem.matw)) # label each matrix with ATOPID


# Capture means and SD in lists of dataframes
growth.pred.mean <- lapply(growth.pred, function(x) as.data.frame(apply(x, 1, mean)))
growth.pred.sd <- lapply(growth.pred, function(x) as.data.frame(apply(x, 1, sd)))


# combine means and sds in a dataframe 
growth.pred.res <- cbind(do.call(rbind, growth.pred.mean), unlist(growth.pred.sd))
growth.pred.res$atopid <- as.numeric(substr(rownames(growth.pred.res), 1, 4))
colnames(growth.pred.res) <- c("mean.pred", "sd.pred", "atopid")
growth.pred.res$mean.pred <- sinh(growth.pred.res$mean.pred) # reverse IHS transformation
growth.pred.res$sd.pred <- sinh(growth.pred.res$sd.pred) # reverse IHS transformation
growth.pred.res$nz.weights <- state.mem.matw[state.mem.matw != 0]


# Create a dataframe with maximum predicted change, positive or negative 
growth.pred.res.max <- growth.pred.res %>% 
  select(atopid, nz.weights, mean.pred) %>% 
  group_by(atopid) %>%
  summarise_each(funs(.[which.max(abs(.))]))  



# plot: hard to read with all the data points
ggplot(growth.pred.res, aes(x = nz.weights, y = mean.pred)) +
  geom_hline(yintercept = 0) +
  geom_point(position = position_jitter(width = 0.1), alpha = .25) + 
  geom_smooth(method = "lm") + 
  labs(x = "Economic Weight", y = "Mean Predicted Military Spending Growth from Alliance") +
  ggtitle("Predicted Military Spending Growth and Treaty Depth") +
  theme_bw() 


# Another way to attack the clear overplotting problem
growth.weight.plot <- ggplot(growth.pred.res, aes(x = nz.weights, y = mean.pred)) +
  geom_hline(yintercept = 0) +
  stat_bin_hex(colour="white", na.rm=TRUE) +
  scale_fill_gradientn(colours=c("#999999","#333333"), 
                       name = "Frequency", 
                       na.value=NA) +
  labs(x = "Economic Weight",
       y = "Mean Predicted Spending Growth from Alliance") +
  ggtitle("Predicted Military Spending Growth and Treaty Depth") +
  theme_bw() 
growth.weight.plot

# take largest absolute value for each alliance
ggplot(growth.pred.res.max, aes(x = nz.weights, y = mean.pred)) +
  geom_hline(yintercept = 0) +
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(x = "Economic Weight", y = "Largest Mean Predicted Military Spending Growth from Alliance") +
  theme_classic() 





# Remove fit model from workspace
saveRDS(ml.model.w, "data/ml-model-w.rds")
rm(ml.model.w)
