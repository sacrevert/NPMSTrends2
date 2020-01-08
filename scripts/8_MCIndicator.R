##############################################
#### 8. Combine trends into an indicator #####
##############################################
# O.L. Pescott, olipes@ceh.ac.uk
# Dec. 2019/Jan. 2020
############
# Data from script 5, scripts/5_processDataForTrendsFunction.R
rm(list=ls())
library(plyr)
library(coda) # required for summmary to function properly on coda.samples object
######################################
load(file = "outputs/sppRuns/grasslandResults_06012020_FULL.Rdata") # loads object sppModels
# Each species in sppModels list has three items:
# 1. x = original data; 2. tableOut = summary table of parameters extracted ('mC', 'mPsi', 'mu', 'annOcc', 'avgOcc')
# 3. coda.samples = actual mcmc samples (1000 iterations)
# Within 2, parameter 'mC' = annual mean cover for taxon including zeros.
# inspect <- sppModels[[1]]

### Function ###
# columns for mean and time-series se for parameter = pattern (pattern can only take one element at the moment)
getMuTse <- function(x, pattern = 'mC') { out <- summary(x[[3]])
                          muTse <- out$stat[,c(1,4)]
                          return(muTse[grep(rownames(muTse), pattern = pattern),]) }
# test
#test <- getMuTse(sppModels[[1]])
#

# Extract all summary stats for mC
mcStats <- lapply(sppModels, function(i) getMuTse(x = i))

### Function ###
# Create indicator, involving
# i. simulate N draws for a species using index mean and s.e.
# ii. normalise per species time series to baseline (year 1 = 2015)
# iii. take geometric mean across species at each time point
# (plot in separate function)

simDist <- function(a, b) rnorm(n = 1000, mean = a, sd = b)
#apply(test, 1, function(y) simDist(y['Mean'],y['Time-series SE'])) # works in isolation
simsList <- list()
simsList <- lapply(mcStats, function(x) apply(x, 1, function(y) simDist(y['Mean'],y['Time-series SE']))) # works in isolation )
simsLNorm <- lapply(simsList, function(x) (x/mean(x[,1]))*100 ) # normalise each species indices to year 1 species mean, scale to 100

# Calculate MSI across species for each simulation
geoMeans <- list() # replicate MSIs
for (i in 1:dim(simsLNorm[[1]])[1]){ # for each sim
  geoMeans[[i]] <- colMeans(do.call(rbind, lapply(simsLNorm, function(y) log(y[i,]))))
}
geoMeansDf <- do.call(rbind, geoMeans)
# exponentiate geometric means of simulations
meansDf <- exp(na.exclude(geoMeansDf)) # 8 NaNs, ignore for the moment, as only 8 out of 1000
avgMSI <- t(apply(meansDf, 2, function(y) c(mean(y), sd(y)/sqrt(1000)))) # means and standard errors across the 1000 (well, 992) MSIs (multi-species indicators)
avgMSI <- round(avgMSI, digits = 2)

### Plot multi-species indicator
library(ggplot2)
###
avgMSIDF <-data.frame(avgMSI)
names(avgMSIDF)[1:2] <- c("MSI", "s.e.")
avgMSIDF$Year <- 2015:2019
avgMSIDF$low <- avgMSIDF$MSI - avgMSIDF$s.e.
avgMSIDF$high <- avgMSIDF$MSI + avgMSIDF$s.e.

ggplot(avgMSIDF, aes_string(x = "Year", y = "MSI")) +
  geom_ribbon(data = avgMSIDF,
              aes_string(group = 1, ymin = "low", ymax = "high"),
              alpha = 0.2) +
  geom_line(size = 1, col = "black") +
  ylab("MSI") +
  xlab("Year") +
  scale_y_continuous(limits = c(75, 125)) +
  theme(plot.title = element_text(lineheight = .8, face = "bold"),
        legend.position = 'bottom')
