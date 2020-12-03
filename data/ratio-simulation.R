# Joshua Alley
# Texas A&M University
# Simulation of correlations with IV in denominator of the DV


# 20 units
sim.data.list <- vector(mode = "list", length = 20)

for(i in 1:length(sim.data.list)){
# Simulate IV with autocorrelation
# purely random process with mean 0 and standard deviation 1.5
error.x <- rnorm(100, mean = 1, sd = .5)

# seed
sim.auto.x <- rnorm(1, mean = 5, sd = .75)

# the process
for (j in 2:length(error.x)) {
  sim.auto.x[j] <- .90*sim.auto.x[j-1] + error.x[j]
}

# simulate ratio
error.dv <- rnorm(100, mean = .1, sd = .05)

# start
sim.auto.dv <- rnorm(1, mean = .2, sd = .05)

# the process
for (k in 2:length(error.dv)){
  sim.auto.dv[k] <- .85*sim.auto.dv[k-1] + error.dv[k]
}

# check correlation from common time trend
print(cor(sim.auto.x, sim.auto.dv))


# create ratio
sim.ratio <- sim.auto.dv / sim.auto.x


# create a dataframe:
sim.ratio.data <- cbind.data.frame(sim.ratio, sim.auto.dv, sim.auto.x)
sim.ratio.data$obs <- as.numeric(rownames(sim.ratio.data))

# assign data to each element of list
sim.data.list[[i]] <- sim.ratio.data

}


# pull simulated data into a list
sim.df.full <- bind_rows(sim.data.list, .id = "unit") %>%
              group_by(obs) # grouping var for later cross-sections
summary(sim.df.full$sim.ratio)
ggplot(sim.df.full, aes(x = sim.ratio)) +
  geom_histogram()



### Use regression and correlation to caculate associations 

# calculate correlations by observation 
sim.ratio.cor <- summarize(sim.df.full,
                           ratio = cor(sim.ratio, sim.auto.x),
                           non.ratio = cor(sim.auto.dv, sim.auto.x)
                          ) 
summary(sim.ratio.cor$ratio)
summary(sim.ratio.cor$non.ratio)


# create data for plotting
sum.cor.long <- pivot_longer(sim.ratio.cor, 
                             -obs, 
                             names_to = "outcome", 
                             values_to = "correlation")

# plot correlations
ggplot(sum.cor.long, aes(x = correlation, fill = outcome)) +
  geom_density(alpha = .6)


# Linear regression results
sim.ratio.lm <- sim.df.full %>%
                 drop_na() %>%
                 summarize(
                          ratio.est = summary(lm(sim.ratio ~ sim.auto.x))[["coefficients"]][2, 1],
                          non.ratio.est = summary(lm(sim.auto.dv ~ sim.auto.x))[["coefficients"]][2, 1],
                          
                          ratio.tstat = summary(lm(sim.ratio ~ sim.auto.x))[["coefficients"]][2, 3],
                          non.ratio.tstat = summary(lm(sim.auto.dv ~ sim.auto.x))[["coefficients"]][2, 3]
                          )

# Plot t-stats
sum.lm.long <- pivot_longer(sim.ratio.lm, 
                              -c(obs, ratio.est, non.ratio.est), 
                                 names_to = "outcome", 
                                 values_to = "tstat") %>%
                          mutate(
                          stat.sig = ifelse(tstat > 2.09 | tstat < -2.09, 1, 0),
                          stat.sig.factor = as.factor(stat.sig)
                          )

# plot densities of t-stats
ggplot(sum.lm.long, aes(x = tstat, fill = outcome)) +
  geom_density(alpha = .6)


# look at stat.sig conclusions
statsig.nocor <- ggplot(sum.lm.long, aes(x = stat.sig.factor)) +
  facet_wrap(~outcome,
             labeller=labeller(outcome = c(non.ratio.tstat = "Non-Ratio",
                                           ratio.tstat = "Ratio"))
  ) +
  geom_bar() +
  theme_bw() +
  labs(
    x = "Statistical Significance",
    y = "Count",
    title = "Uncorrelated Variables") +
  scale_x_discrete(labels=c("0" = "No",
                            "1" = "Yes"))
statsig.nocor

# plot test-statistics with more detail
tstat.nocor <- ggplot(sum.lm.long, aes(y = tstat, x = as.factor(outcome))) +
                geom_hline(yintercept = 2.09) +
                geom_hline(yintercept = -2.09) +
                geom_boxplot(outlier.shape = NA) +
                geom_jitter() +
                theme_bw() +
                 labs(
                   x = "Outcome Measure",
                   y = "Test Statistic Estimate",
                   title = "Uncorrelated Variables"
                   ) +
                scale_x_discrete(labels=c("non.ratio.tstat" = "Non-Ratio",
                                   "ratio.tstat" = "Ratio"))
tstat.nocor

# tabulate results
table(sum.lm.long$outcome, sum.lm.long$stat.sig.factor)
chisq.test(x = sum.lm.long$outcome, y = sum.lm.long$stat.sig.factor)



### Same simulation, but with a correlation between the two outcomes


# look at the correlation between the outcomes to inform the simulation 
# full data filtered to match the sample
state.data19 <- filter(state.vars, year >= 1919)

# milex and gdp are correlated 
cor.test(state.data19$ln.milex, state.data19$ln.gdp)
cor.test(state.data19$change.ln.milex, state.data19$ln.gdp)
cor.test(state.data19$growth.milex, state.data19$ln.gdp)

ggplot(state.data19, aes(x = gdp, y = milex)) + geom_point()
state.data19 %>%
  filter(ln.gdp > 15) %>% # cut microstates for clarity
  ggplot(aes(x = ln.gdp, y = ln.milex)) + geom_point()


# 20 units
sim.data.cor <- vector(mode = "list", length = 20)

for(i in 1:length(sim.data.cor)){
  
  # Simulate correlated predictor and outcome
  burn.in <- 300
  n <- 100 
  rho <- cor(state.data19$ln.milex, state.data19$ln.gdp, 
             use = "pairwise.complete.obs") 
  rho1 <- 0.9
  rho2 <- 0.85
  q12 <- calcrho(rho, rho1, rho2)
  eps <- mvrnorm(n+burn.in, mu=c(5, .75), Sigma=cbind(c(1, q12),c(q12, 1)))
  
  sim.auto.x <- arima.sim(list(ar=rho1), n, innov = eps[burn.in+1:n, 1], start.innov = eps[1:burn.in,1])
  sim.auto.dv <- arima.sim(list(ar=rho2), n, innov = eps[burn.in+1:n, 2], start.innov = eps[1:burn.in,2])
  print(cor(sim.auto.x, sim.auto.dv))
  
  
  # create ratio
  sim.ratio <- sim.auto.dv / sim.auto.x
  
  # create a dataframe:
  sim.ratio.data <- cbind.data.frame(sim.ratio, sim.auto.dv, sim.auto.x, sim.growth)
  sim.ratio.data$obs <- as.numeric(rownames(sim.ratio.data))
  
  # assign data to each element of list
  sim.data.cor[[i]] <- sim.ratio.data
  
}


# pull simulated data into a list
sim.df.cor <- bind_rows(sim.data.cor, .id = "unit") %>%
  group_by(obs) # switch grouping var for later cross-sections

# replace any negatives with zero: actual value for some states
summary(sim.df.cor$sim.ratio)
sim.df.cor$sim.ratio[sim.df.cor$sim.ratio < 0] <- 0
summary(sim.df.cor$sim.ratio)

ggplot(sim.df.cor, aes(x = sim.ratio)) +
  geom_histogram()


# Calculate regressions on cross-sections of observations
# Linear regression results
sim.lm.cor <- sim.df.cor %>%           
           drop_na() %>%
           summarize(
              ratio.est = summary(lm(sim.ratio ~ sim.auto.x))[["coefficients"]][2, 1],
              non.ratio.est = summary(lm(sim.auto.dv ~ sim.auto.x))[["coefficients"]][2, 1],
    
              ratio.tstat = summary(lm(sim.ratio ~ sim.auto.x))[["coefficients"]][2, 3],
              non.ratio.tstat = summary(lm(sim.auto.dv ~ sim.auto.x))[["coefficients"]][2, 3]
            )

# Plot t-stats
sum.lm.cor <- pivot_longer(sim.lm.cor, 
                            -c(obs, ratio.est, non.ratio.est), 
                            names_to = "outcome", 
                            values_to = "tstat") %>%
  mutate(
    stat.sig = ifelse(tstat > 2.09 | tstat < -2.09, 1, 0),
    stat.sig.factor = as.factor(stat.sig)
  )

# plot densities of t-stats
ggplot(sum.lm.cor, aes(x = tstat, fill = outcome)) +
  geom_density(alpha = .6)

# look at stat.sig conclusions
statsig.cor <- ggplot(sum.lm.cor, aes(x = stat.sig.factor)) +
                facet_wrap(~outcome,
                labeller=labeller(outcome = c(non.ratio.tstat = "Non-Ratio",
                                           ratio.tstat = "Ratio"))
                  ) +
                 geom_bar() +
                theme_bw() +
                labs(
                x = "Statistical Significance",
                y = "Count",
                title = "Correlated Variables") +
               scale_x_discrete(labels=c("0" = "No",
                            "1" = "Yes"))
statsig.cor

table(sum.lm.cor$outcome, sum.lm.cor$stat.sig.factor)
chisq.test(x = sum.lm.cor$outcome, y = sum.lm.cor$stat.sig.factor)

# plot test-statistics with more detail
tstat.cor <- ggplot(sum.lm.cor, aes(y = tstat, x = as.factor(outcome))) +
               geom_hline(yintercept = 2.09) +
               geom_hline(yintercept = -2.09) +
               geom_boxplot(outlier.shape = NA) +
               geom_jitter() +
               theme_bw() +
                labs(
                 x = "Outcome Measure",
                 y = "Test Statistic Estimate",
                 title = "Correlated Variables"
                 ) +
               scale_x_discrete(labels=c(
                        "non.ratio.tstat" = "Non-Ratio",
                        "ratio.tstat" = "Ratio"))
tstat.cor



### Use growth as the outcome
# pull simulated data into a list and calculate growth
sim.df.growth <- bind_rows(sim.data.cor, .id = "unit") %>%
  group_by(unit) %>%
  mutate(
    growth.dv = sim.auto.dv / lag(sim.auto.dv),
    growth.x = sim.auto.x / lag(sim.auto.x)
  ) %>%
  group_by(obs) # switch grouping var for later cross-sections

summary(sim.df.growth$growth.dv)
summary(sim.df.growth$growth.x)


# Calculate regressions on cross-sections of observations
# Linear regression results
sim.lm.growth <- sim.df.growth %>%           
  drop_na() %>%
  summarize(
    ratio.est = summary(lm(sim.ratio ~ sim.auto.x))[["coefficients"]][2, 1],
    growth.est = summary(lm(growth.dv ~ sim.auto.x))[["coefficients"]][2, 1],
    
    ratio.tstat = summary(lm(sim.ratio ~ sim.auto.x))[["coefficients"]][2, 3],
    growth.tstat = summary(lm(growth.dv ~ sim.auto.x))[["coefficients"]][2, 3]
  )

# Plot t-stats
sum.lm.growth <- pivot_longer(sim.lm.growth, 
                           -c(obs, ratio.est, growth.est), 
                           names_to = "outcome", 
                           values_to = "tstat") %>%
  mutate(
    stat.sig = ifelse(tstat > 2.09 | tstat < -2.09, 1, 0),
    stat.sig.factor = as.factor(stat.sig)
  )

# plot densities of t-stats
ggplot(sum.lm.growth, aes(x = tstat, fill = outcome)) +
  geom_density(alpha = .6)

# look at stat.sig conclusions
statsig.growth <- ggplot(sum.lm.growth, aes(x = stat.sig.factor)) +
  facet_wrap(~outcome,
             labeller=labeller(outcome = c(growth.tstat = "Percentage Change",
                                           ratio.tstat = "Ratio"))
  ) +
  geom_bar() +
  theme_bw() +
  labs(
    x = "Statistical Significance",
    y = "Count",
    title = "Percentage Change vs Ratio") +
  scale_x_discrete(labels=c("0" = "No",
                            "1" = "Yes")) 
statsig.growth

table(sum.lm.growth$outcome, sum.lm.growth$stat.sig.factor)
chisq.test(x = sum.lm.growth$outcome, y = sum.lm.growth$stat.sig.factor)

# plot test-statistics with more detail
tstat.growth <- ggplot(sum.lm.growth, aes(y = tstat, x = as.factor(outcome))) +
  geom_hline(yintercept = 2.09) +
  geom_hline(yintercept = -2.09) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  theme_bw() +
  labs(
    x = "Outcome Measure",
    y = "Test Statistic Estimate",
    title = "Percentage Change vs Ratio"
  ) +
  scale_x_discrete(labels=c(
    "growth.tstat" = "Percentage Change",
    "ratio.tstat" = "Ratio"))
tstat.growth







# Combine plots and add to the appendix
grid.arrange(tstat.nocor, tstat.cor, tstat.growth, ncol = 3)
grid.arrange(statsig.nocor, statsig.cor, statsig.growth, ncol = 3)

sim.inferences <- arrangeGrob(statsig.nocor, statsig.cor, statsig.growth,
                               ncol = 3)
ggsave("appendix/sim-inferences.png", sim.inferences, height = 6, width = 8.4) #save file
