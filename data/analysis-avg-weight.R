# Joshua Alley
# Texas A&M University
# Empirical test of public goods theory of alliances
# Examine corelation between economic weight and military spending 


# Load state-ally-year data 
state.ally.year <- read.csv("data/alliance-state-year.csv")
state.vars <- read.csv("data/state-vars.csv")


# summarize state-ally-year data with key indicators
state.ally.year$treaty.pres <- ifelse(state.ally.year$atopid > 0, 1, 0)

state.ally.sum <- state.ally.year %>%
  group_by(ccode, year) %>%
  summarize(
    obs = n(),
    total.ally.expend = sum(ally.spend[defense == 1 | offense == 1], na.rm = TRUE),
    avg.treaty.weight = mean(alliance.contrib[defense == 1 | offense == 1], na.rm = TRUE),
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
state.char.full[, 42: ncol(state.char.full)][is.na(state.char.full[, 44: ncol(state.char.full)])] <- 0




# Log-transform allied expenditures variable
state.char.full$ln.ally.expend <- log(state.char.full$total.ally.expend + 1)

# IHS transforamtion of growth variable
state.char.full$ihs.growth.milex <- asinh(state.char.full$growth.milex)

# Difference allied spending: consistent with approach of Plumper and Neumayer
state.char.full <- state.char.full %>%
                  group_by(ccode) %>%
                 mutate(
                   lag.ally.expend = lag(ln.ally.expend),
                   diff.ally.expend = ln.ally.expend - lag.ally.expend
                 ) 


# Normalized GDP by year
state.char.full <- state.char.full %>% 
  group_by(year) %>%
  mutate(
    ln.gdp.norm = ln.gdp / max(ln.gdp, na.rm = TRUE)
  ) %>%
  filter(year >= 1919) %>%
  group_by()

# Compare the two variables 
ggplot(state.char.full, aes(x = ln.gdp)) + geom_histogram()
ggplot(state.char.full, aes(x = ln.gdp.norm)) + geom_histogram()




# plot on ln shows clear evidence of non-stationarity
ggplot(state.char.full, aes(x = year, y = ln.milex, group = as.factor(ccode))) + 
  geom_line(alpha = .5)

# Growth is more mean-reverting
ggplot(state.char.full, aes(x = year, y = asinh(growth.milex), group = as.factor(ccode))) + 
  geom_line(alpha = .5)


# filter out states:
state.char.sample <- filter(state.char.full, ln.gdp >= 18 & treaty.pres > 0)


### First test: 
m1.pg.abs <- rlm(growth.milex ~ avg.treaty.weight + ln.gdp.norm + 
                     avg.num.mem + avg.dem.prop + 
                     atwar + civilwar.part + polity  +
                     lsthreat + cold.war,
                   data = state.char.sample
                )
summary(m1.pg.abs)
plotreg(m1.pg.abs)

# OLS
m2.pg.abs <- lm(growth.milex ~ avg.treaty.weight + ln.gdp.norm + 
                       avg.num.mem + avg.dem.prop + 
                       atwar + civilwar.part + polity + 
                       lsthreat + cold.war,
                     data = state.char.sample)
summary(m2.pg.abs)




### Transform outcome with Inverse hyperbolic sine to reduce the extremely heavy tails 
m1.abs.ihs <- rlm(ihs.growth.milex ~ avg.treaty.weight + ln.gdp.norm + 
                   avg.num.mem + avg.dem.prop + 
                   atwar + civilwar.part + polity  +
                   lsthreat + cold.war,
                 data = state.char.sample
)
summary(m1.abs.ihs)
plotreg(m1.abs.ihs)


# Results are very similar- check distribution of weights across observations
plot(m1.pg.abs$residuals, m1.pg.abs$w) # original model
abline(h = 0)
plot(m1.abs.ihs$residuals, m1.abs.ihs$w) # transformed DV


# OLS estimation: positive driven by unusual obs
m2.abs.ihs <- lm(ihs.growth.milex ~ avg.treaty.weight + ln.gdp.norm + 
                    avg.num.mem + avg.dem.prop + 
                    atwar + civilwar.part + polity  +
                    lsthreat + cold.war,
                  data = state.char.sample
)
summary(m2.abs.ihs)


# Take this and ols estiamtes in one appendix table
stargazer(m1.pg.abs, m2.pg.abs, m2.abs.ihs, m1.abs.ihs,
          style = "all2",
          dep.var.labels = c("\\% Change Milex.", "IHS(\\% Change Milex.)"),
          covariate.labels = c("Avg. Economic Weight",
                               "ln(GDP)", "Avg. Alliance Size",
                               "Avg. Allied Democracy", "International War",
                               "Civil War Participant", "Regime Type",
                               "External Threat", "Cold War"),
          ci=TRUE, 
          star.char = c("", "", ""),
          notes = "95\\% Confidence Intervals in Parentheses.", 
          notes.append = FALSE,
          label = c("tab:avg-weight-res"))


### calculate likely effect of increasing average weight
summary(state.char.sample$avg.treaty.weight)

# Estimation sample
est.data.noint <- select(state.char.sample, c(growth.milex,  avg.treaty.weight, 
                          ln.gdp.norm, avg.num.mem, avg.dem.prop, 
                           atwar, civilwar.part, polity,
                           lsthreat, cold.war)) %>%
                   drop_na()

# Bootstrap
n.bs <- 5000
bs.sim <- matrix(NA, nrow = n.bs, ncol = length(coef(m1.pg.abs)))
for (i in 1:n.bs) {
  bs.samp <- sample(1:nrow(est.data.noint), nrow(est.data.noint), replace = T)
  bs.data <- est.data.noint[bs.samp, ]
  bs.est <- rlm(growth.milex ~ avg.treaty.weight + ln.gdp.norm + 
                  avg.num.mem + avg.dem.prop + 
                  atwar + civilwar.part + polity  +
                  lsthreat + cold.war,
                data = bs.data, maxit = 100000)
  bs.sim[i, ] <- coef(bs.est)
}
sim.res <- bs.sim


# Create vectors to store hypothetical values
x.low.weight <- numeric(ncol(est.data.noint))
names(x.low.weight) <- c("Intercept", colnames(select(est.data.noint, -c(growth.milex))))

# Set values for no pcj case
x.low.weight["Intercept"] <- 1
x.low.weight["avg.treaty.weight"] <- 0.05 # first quartile
x.low.weight["ln.gdp.norm"] <- median(est.data.noint$ln.gdp.norm)
x.low.weight["avg.num.mem"] <- median(est.data.noint$avg.num.mem) 
x.low.weight["avg.dem.prop"] <- median(est.data.noint$avg.dem.prop) 
x.low.weight["atwar"] <- 0 # no war
x.low.weight["civilwar.part"] <- 0 # no civil war part
x.low.weight["polity"] <- median(est.data.noint$polity)
x.low.weight["lsthreat"] <- median(est.data.noint$lsthreat)
x.low.weight["cold.war"] <- 1 # cold-war

x.low.weight

# Create vector with pcj
x.high.weight <- x.low.weight
x.high.weight["avg.treaty.weight"] <- .39 # third quartile

# Simulate quantities of interest
eval.low.weight <- sim.res%*%x.low.weight
eval.high.weight <- sim.res%*%x.high.weight

# difference in two scenarios
eval.weight.diff <- eval.high.weight - eval.low.weight
# calculate 90% credible interval
weight.diff.int <- quantile(eval.weight.diff, c(.05, .5, .95))
weight.diff.int 
plot(density(eval.weight.diff))


