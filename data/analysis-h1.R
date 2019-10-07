# Joshua Alley
# Texas A&M University
# Empirical test of public goods theory of alliances
# Examine Hypothesis 1


# Load packages
library(here)
library(arm)
library(reshape2)
library(MASS)
library(plm)
library(texreg)
library(interflex)
library(margins)
library(stargazer)
library(robustlmm)
library(tidyverse)



# Set working directory to current folder 
setwd(here::here())
getwd()


# set seed
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
  group_by()

# Compare the two variables 
ggplot(state.char.full, aes(x = ln.gdp)) + geom_histogram()
ggplot(state.char.full, aes(x = ln.gdp.norm)) + geom_histogram()



# Set up interaction with a separate dataset
state.char.inter <- as.data.frame(state.char.full) # interflex doesn't take tibble input
state.char.inter <- filter(state.char.inter, ln.gdp >= 18) # filter really tiny states out

# summarize interaction variables
summary(state.char.inter$ln.gdp)
ggplot(state.char.inter, aes(x = ln.gdp)) + geom_histogram()

# normalize GDP
summary(state.char.inter$ln.gdp.norm)
ggplot(state.char.inter, aes(x = ln.gdp.norm)) + geom_histogram()


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
ggplot(state.char.inter, aes(x = year, y = ln.milex, group = as.factor(ccode))) + 
  geom_line(alpha = .5)

# Growth is more mean-reverting
ggplot(state.char.inter, aes(x = year, y = asinh(growth.milex), group = as.factor(ccode))) + 
  geom_line(alpha = .5)



### First test: aboslute size (GDP)
# Interact changes in allied spending and GDP
# Total allied spending: pooling regression
m1.pg.abs <- rlm(growth.milex ~ diff.ally.expend + ln.gdp.norm + diff.ally.expend:ln.gdp.norm +
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
cplot(m1.pg.abs, x = "ln.gdp.norm", dx = "diff.ally.expend", what = "effect",
      main = "Allied Spending and Growth in Military Spending",
      xlab = "ln(GDP): Normalized", ylab = "Average M.E. of Changes in Allied Spending")
abline(h = 0)
# Export plot
dev.copy(pdf,'manuscript/abs-margins-plot.pdf')
dev.off()

# Switch the directions 
cplot(m1.pg.abs, x = "diff.ally.expend", dx = "ln.gdp.norm", what = "effect",
      main = "Marginal Effect GDP on Growth in Military Spending",
      xlab = "Change in Allied Capability", ylab = "Average M.E. of ln(GDP)")
abline(h = 0)



# OLS
m2.pg.abs <- lm(growth.milex ~ diff.ally.expend + ln.gdp.norm + diff.ally.expend:ln.gdp.norm +
                       avg.num.mem + avg.dem.prop + 
                       atwar + civilwar.part + polity + 
                       lsthreat + cold.war,
                     data = state.char.inter)
summary(m2.pg.abs)



# binning estimator
bin.abs <- inter.binning(Y = "ihs.growth.milex", D = "diff.ally.expend", X = "ln.gdp.norm", 
              Z = c("lag.ln.milex", "atwar", "civilwar.part", "polity", 
                    "lsthreat", "cold.war", "avg.num.mem", "avg.dem.prop"), 
              data = state.char.inter,
              na.rm = TRUE,
              Ylabel = "Growth in Military Spending",
              Dlabel = "Allied Spending",
              Xlabel = "ln(GDP): Normalized", theme.bw = TRUE
)
bin.abs
ggsave("appendix/inter-bin-abs.pdf", height = 6, width = 8)


# Kernel: 10+ minute run time 
# Subset to values above .7: insane uncertainty otherwise
kernel.abs <- inter.kernel(Y = "growth.milex", D = "diff.ally.expend", X = "ln.gdp.norm", 
             Z = c("lag.ln.milex", "atwar", "civilwar.part", "polity", 
                   "lsthreat", "cold.war", "avg.num.mem", "avg.dem.prop"), 
             data = subset(state.char.inter, state.char.inter$ln.gdp.norm > .7), 
             na.rm = TRUE,
             nboots = 200, parallel = TRUE, cores = 4,
             Ylabel = "Growth in Military Spending",
             Dlabel = "Allied Spending",
             Xlabel = "ln(GDP): Normalized", theme.bw = TRUE
)
kernel.abs
ggsave("appendix/inter-kernel-abs.pdf", height = 6, width = 8)



### Transform outcome with Inverse hyperbolic sine to reduce the extremely heavy tails 
# Total allied spending: pooling regression
m1.abs.ihs <- rlm(ihs.growth.milex ~ diff.ally.expend + ln.gdp.norm + diff.ally.expend:ln.gdp.norm +
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
cplot(m1.abs.ihs, x = "ln.gdp.norm", dx = "diff.ally.expend", what = "effect",
      main = "Marginal Effect of Changes in Allied Spending on Growth in Military Spending",
      xlab = "ln(GDP): Normalized", ylab = "Average M.E. of Changes in Allied Spending")
abline(h = 0)

# Switch the directions 
cplot(m1.abs.ihs, x = "diff.ally.expend", dx = "ln.gdp.norm", what = "effect",
      main = "Marginal Effect GDP on Growth in Military Spending",
      xlab = "Change in Allied Capability", ylab = "Average M.E. of ln(GDP)")
abline(h = 0)


# Results are very similar- check distribution of weights across observations
plot(m1.pg.abs$residuals, m1.pg.abs$w) # original model
abline(h = 0)
plot(m1.abs.ihs$residuals, m1.abs.ihs$w) # transformed DV


# OLS estiamtion 
m2.abs.ihs <- lm(asinh(growth.milex) ~ diff.ally.expend + ln.gdp.norm + diff.ally.expend:ln.gdp.norm +
                    avg.num.mem + avg.dem.prop + 
                    atwar + civilwar.part + polity  +
                    lsthreat + cold.war,
                  data = state.char.inter
)
summary(m2.abs.ihs)

cplot(m2.abs.ihs, x = "ln.gdp.norm", dx = "diff.ally.expend", what = "effect",
      main = "Marginal Effect of Changes in Allied Spending on Growth in Military Spending",
      xlab = "ln(GDP): Normalized", ylab = "Average M.E. of Changes in Allied Spending")
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
cplot(m1.pg.abs, x = "ln.gdp.norm", dx = "diff.ally.expend", what = "effect",
      main = "Robust Reg",
      xlab = "ln(GDP): Normalized", ylab = "Average M.E. of Allied Cap")

# rreg- transformed DV
cplot(m1.abs.ihs, x = "ln.gdp.norm", dx = "diff.ally.expend", what = "effect",
      main = "Robust Reg: IHS",
      xlab = "ln(GDP): Normalized", ylab = "Average M.E. of Allied Cap")

# OLS 
cplot(m2.pg.abs, x = "ln.gdp.norm", dx = "diff.ally.expend", what = "effect",
      main = "Linear Regression",
      xlab = "ln(GDP): Normalized", ylab = "Average M.E. of Allied Cap")

# OLS with transformed DV 
cplot(m2.abs.ihs, x = "ln.gdp.norm", dx = "diff.ally.expend", what = "effect",
      main = "Linear Regression: IHS",
      xlab = "ln(GDP): Normalized", ylab = "Average M.E. of Allied Cap")
# Export plot
dev.copy(pdf,'appendix/me-plots.pdf')
dev.off()
dev.off()




# Robust regression with random effects: I may need to RSTAN or brms this thing.
# This is freezing both machines at the moment: comment out. 
#rreg.re.abs <- rlmer(growth.milex ~ diff.ally.expend + ln.gdp.norm + diff.ally.expend:ln.gdp.norm +
#  avg.num.mem + avg.dem.prop + 
#  atwar + civilwar.part + polity  +
#  lsthreat + cold.war + (1|ccode) + (1|year),
#data = state.char.inter
#)
#summary(rreg.re.abs)

# No cplot or stargazer to report these results
#rreg.re.sum <- summary(rreg.re.abs)
#stargazer(rreg.re.sum$coefficients)


### Additional single-level test: relative size expressed as contribution to alliance
# estimate interactions
# filter out cases with no alliances
inter.data.rel <- filter(state.char.full, treaty.pres == 1)
inter.data.rel <- as.data.frame(inter.data.rel)

# Look at robust regression among alliance members
m1.pg.all <- rlm(growth.milex ~ diff.ally.expend + ln.gdp.norm + diff.ally.expend:ln.gdp.norm +
                   lag.ln.milex + avg.num.mem + avg.dem.prop + 
                   atwar + civilwar.part + polity  + 
                   lsthreat + cold.war,
                 data = inter.data.rel)
summary(m1.pg.all)

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
m2.pg.rel <- lm(ihs.growth.milex ~ diff.ally.expend + avg.treaty.contrib + diff.ally.expend:avg.treaty.contrib +
                     avg.dem.prop + lag.ln.milex +
                     atwar + civilwar.part + polity + ln.gdp.norm + avg.num.mem +
                     lsthreat + cold.war,
                   data = inter.data.rel
                )
summary(m2.pg.rel)
plot(rreg.re.abs)

# binning estimator
bin.rel <- inter.binning(Y = "ihs.growth.milex", D = "diff.ally.expend", X = "avg.treaty.contrib", 
                         Z = c("lag.ln.milex", "atwar", "civilwar.part", "polity", 
                               "lsthreat", "cold.war", "avg.num.mem", "avg.dem.prop"), 
                         data = inter.data.rel,
                         na.rm = TRUE
)
bin.rel 

# Kernel: 10+ minute run time 
kernel.rel <- inter.kernel(Y = "ihs.growth.milex", D = "diff.ally.expend", X = "avg.treaty.contrib", 
                           Z = c("lag.ln.milex", "atwar", "civilwar.part", "polity", 
                                 "lsthreat", "cold.war", "avg.num.mem", "avg.dem.prop"), 
                           data = inter.data.rel, 
                           na.rm = TRUE,
                           nboots = 200, parallel = TRUE, cores = 4
)
kernel.rel


