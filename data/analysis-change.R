# Joshua Alley
# Analysis of Spending changes

# add changes in spending 
reg.state.data.change <- left_join(reg.state.data, 
                                   select(state.vars,
                                          ccode, year,
                                          change.ln.milex, 
                                          lag.ln.milex))

reg.state.mat.change <- as.matrix(cbind(reg.state.data.change[219],
                                reg.state.data[, 4:11]))


# Define the data list 
stan.data.change <- list(N = nrow(reg.state.data), y = reg.state.data.change$change.ln.milex,
                  state = reg.state.data$state.id, S = length(unique(reg.state.data$state.id)),
                  year = reg.state.data$year.id, T = length(unique(reg.state.data$year.id)),
                  A = ncol(state.mem.mat),
                  X = reg.state.mat.change, M = ncol(reg.state.mat.change),
                  Z = state.mem.mat
)


# Compile the model code
model.change <- stan_model(file = "data/ml-model-noasinh.stan")

# Variational Bayes- use to check model will run before going full STAN 
vb.change <- vb(model.change, data = stan.data.change, seed = 12)
# Does not converge: diagnostics suggest problems

# Clear variational Bayes results from environment
rm(vb.change)


# Run model with full Bayes
system.time(
  ml.model.change <- sampling(model.change, data = stan.data.change, 
                       iter = 2800, warmup = 1400, chains = 4
  )
)

# diagnose full model
# launch_shinystan(ml.model.change) 
check_hmc_diagnostics(ml.model.change)


# Extract coefficients from the model
model.sum.change <- extract(ml.model.change, pars = c("gamma", "theta",
                                            "sigma_all"),
                       permuted = TRUE)


# Summarize gamma
gamma.summary.c <- summary(ml.model.change, pars = c("gamma"), probs = c(0.05, 0.95))$summary
gamma.summary.c <- cbind.data.frame(as.numeric(colnames(state.mem.mat)), gamma.summary.c)
colnames(gamma.summary.c) <- c("atopid", "gamma.mean", "gamma.se.mean",
                               "gamma.sd", "gamma.5", "gamma.95",
                               "gamma.neff", "gamma.rhat")


# Ignore uncertainty in estimates: are posterior means positive or negative? 
gamma.summary.c$positive.lmean <- ifelse(gamma.summary.c$gamma.mean > 0, 1, 0)
sum(gamma.summary.c$positive.lmean)
gamma.summary.c$negative.lmean <- ifelse(gamma.summary.c$gamma.mean < 0, 1, 0)
sum(gamma.summary.c$negative.lmean) 


# Plot posterior means of alliance coefficients
ggplot(gamma.summary.c, aes(x = gamma.mean)) +
  geom_histogram(bins = 50) + theme_classic() +
  labs(x = "Posterior Mean", y = "Number of Alliances") +
  ggtitle("Distribution of Alliance Coefficient Posterior Means: Change in Spending")


# Plot points with error bars by ATOPID
ggplot(gamma.summary.c, aes(x = atopid, y = gamma.mean)) +
  geom_errorbar(aes(ymin = gamma.5, 
                    ymax = gamma.95,
                    width=.01), position = position_dodge(0.1)) +
  geom_point(position = position_dodge(0.1)) + geom_hline(yintercept = 0) +
  theme_classic() + coord_flip()


# Join alliance coefficients with ATOP data
alliance.coefs.c <- left_join(atop, gamma.summary.c) %>%
  filter(begyr >= 1898)


# Plot by start year of alliance
mean(model.sum.change$theta)
ggplot(alliance.coefs.c, aes(x = begyr, y = gamma.mean)) +
  geom_errorbar(aes(ymin = gamma.5, 
                    ymax = gamma.95,
                    width=.01), position = position_dodge(0.01)) +
  geom_point(position = position_dodge(0.01)) + 
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = mean(ml.model.sum$theta), linetype = "dashed") +
  labs(x = "Start Year of Alliance", y = "Economic Weight Parameter") +
  theme_classic()
ggsave("appendix/alliance-coefs-year-change.png", height = 6, width = 8)



# save model and remove from workspace
saveRDS(ml.model.change, "data/ml-model-change.rds")
rm(ml.model.change)
