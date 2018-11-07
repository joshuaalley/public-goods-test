# Joshua Alley
# Texas A&M University
# Empirical test of public goods theory of alliances


# Load packages
library(here)
library(MASS)
library(plm)
library(texreg)
library(interflex)
library(sampleSelection)
library(margins)
library(tidyverse)
library(rstan)
library(shinystan)



# Set working directory to current folder 
setwd(here::here())
getwd()



# Load state-ally-year data 
state.ally.year <- read.csv("data/alliance-state-year.csv")
state.vars <- read.csv("data/state-vars.csv")



# summarize state-ally-year data with key indicators

state.ally.year$treaty.pres <- ifelse(state.ally.year$atopid > 0, 1, 0)

state.ally.sum <- state.ally.year %>%
  group_by(ccode, year) %>%
  summarize(
    treary.count = n(),
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


state.char.full <- state.char.full[complete.cases(state.char.full$ccode), ]
state.char.full <- unique(state.char.full) 



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
# Interact total allied spenidng and GDP
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




