# Joshua Alley

# execute all files in public-goods-test

# load packages and set seeds, and manage conflicts
source("data/setup-script.R")
# simulate problems with ratio outcomes 
source("data/ratio-simulation.R")
# primary analyses
source("data/analysis.R")
# analyses of weighted membership and changes
source("data/analysis-z-weight.R")
source("data/analysis-change.R")
# single-level regression results
source("data/analysis-avg-weight.R")
# simulate ML model 
source("data/simulation-check.R")