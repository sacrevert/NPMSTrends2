################################################################################
#### 5b. Create function to make model runs across multiple species easier.#####
################################################################################
# O.L. Pescott, olipes@ceh.ac.uk
# 09.01.2019
rm(list=ls())
# v0 -- original version
# v1 -- amnend function to remove all plots with no known states at time of processing
######################################
list.of.packages <- c("R2jags", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(R2jags) ## obviously requires that JAGS (preferably 4.3.0) is installed
library(dplyr)

# files for development
load(file = "outputs/grasslandsEg_26 03 2019.Rdata") # load saved files from 4_extractData.R

## Helper function (deletes plots without any known states)
#x %>% group_by(plot_id) %>% filter(all(!is.na(dominUnify)))
ilt<-function(x){ # inverse logit function
  exp(x)/(1+exp(x))
}
logit <- function(x){ # logit function
  log(x/(1 - x))
}
######################################

###### Future improvments (noted 10.01.2019)
# Have function where is npms15_18plots and npms15_18spp are not in the environment,
# then files loaded and some processing from 4_extractData done?
# Also need to think about species that are present in aggregates and separately
# e.g. Thymus polytrichus and Ulex gallii -- currently I am creating separate trends for these
######

# domin things
domins <- read.csv(file = "data/dominScores.csv", header = T, stringsAsFactors = F)
#domins <- read.csv("W:/PYWELL_SHARED/Pywell Projects/BRC/_BRC_projects/NPMS/Analyses/2018 08 - Per species trend analyses/r_proj/NPMStrends/data/dominScores.csv", header = T, stringsAsFactors = F)
# t.per <- c(0.001,0.01,0.03,0.05,0.1,0.25,0.33,0.50,0.75,0.95,0.99) these are the cutpoints used in Pescott et al. 2016

## Required data 1
#N = N # the total number of plots
#Y = Y # the total number of years
#n.Plot.pos = n.Plot.pos # the total number of positive cover observations across all plot visits (= samples in NPMS database terms)
#cpos.Cens = cpos.Cens # indicator (is censored?) -- can be T/F or 1/0, doesn't matter as long as is.logical(cpos.Cens) == T
## note that for the NPMS all non-zero observations are censored; this might not be true though if we combined with CS data or similar with more "precise" cover data
#cpos.Latent = cpos.Latent # just NA values for latent observations (unknown value underlying censored observation)
#lims = lims # the limits to the intervals used for censoring (again, these should be of one type only, unless data collected under different schemes are combined) -- see Pescott et al. (2016)

##########################################################
#### Make function that can be applied across species ####
##########################################################
# Reduced list of sample data, excluding species with no data #
spForMods <- sppDatList[!names(sppDatList) %in% excludedSpp]
#
## Function for preparing a species data for JAGS model and running model
runModels_v1 <- function(i) {
  x <- spForMods[[i]]
  x <- (x %>% group_by(plot_id) %>% filter(any(!is.na(dominUnify))) %>% as.data.frame()) # added in v1: keep plots with at least on non NA only
  x$year <- format(x$date.x, "%Y")
  x <- x[order(x$year, x$plot_id), ]
  uniPlots <- unique(x$plot_id) # unique plot IDs - 31/12/2018 = 933
  # create unique plot index (useful for cross-referencing later on?)
  plotIndex <- data.frame(plot = uniPlots, index = 1:length(uniPlots))
  N <- length(uniPlots)
  # check that all years are in the range of 2015:System time
  if ( max(unique(format(x$date.x, "%Y"))) > format(Sys.time(), "%Y") ) {
    print(paste("Error. Years violation: ", names(spForMods)[i]))
  } else {
    print(paste("Years ok: ", names(spForMods)[i]))
  }
  Y <- length(unique(format(x$date.x, "%Y")))
  # All Domin data, including zeros
  yOrig <- x$dominUnify
  n.Plot.pos <- length(x$dominUnify[x$dominUnify !='0' & !is.na(x$dominUnify)]) # 300
  cpos.Cens <- rep(1, n.Plot.pos)
  cpos.Latent <- rep(NA, n.Plot.pos)
  t <- c(1e-4, 0.01,
         0.01, 0.03,
         0.03, 0.05,
         0.05, 0.1,
         0.1, 0.25,
         0.25, 0.33,
         0.33, 0.5,
         0.5, 0.75,
         0.75, 0.95, 
         0.95, 0.9999)
  tdf <- as.data.frame(matrix(t, nrow = 10, ncol = 2, byrow = TRUE))
  colnames(tdf) <- c('L','U') # 'L'ower and 'U'pper bounds of categories
  tdf$int <- c(1,2,3,4,5,6,7,8,9,10)
  spPos <- x[x$dominUnify !='0' & !is.na(x$dominUnify),]
  # add row number to allow resorting after merge (merge with sort = F does not actually do what we want, could also used plyer::join)
  spPos$indexPos <- 1:nrow(spPos) 
  spPos <- merge(spPos, tdf, by.x = "dominUnify", by.y = "int", all.x = T, all.y = F)
  spPos <- spPos[order(spPos$indexPos),]
  lims <- spPos[,c("L","U")] # has the lower/upper cutpoints for all intervals
  spPos <- merge(spPos, plotIndex, by.x = "plot_id", by.y = "plot", all.x = T, all.y = F)
  spPos <- spPos[order(spPos$indexPos),]
  plot <- spPos$index
  spPos$year <- as.factor(spPos$year)
  levels(spPos$year) <- 1:length(unique(format(x$date.x, "%Y")))
  year <- spPos$year # an indicator linking a percentage cover observation to its parent year
  V2 <- nrow(x) # the total number of visits (samples), irrespective of whether there is a positive cover for a species or not
  #plotZ <- match(Achi_mill_PAN$plot_id, uniPlots)  # an indicator linking a visit to its parent (spatially unique) plot
  x <- merge(x, plotIndex, by.x = "plot_id", by.y = "plot", all.x = T, all.y = F)
  x <- x[order(x$year, x$plot_id), ]
  plotZ <- x$index
  x$year <- as.factor(x$year)
  levels(x$year) <- 1:length(unique(format(x$date.x, "%Y")))
  yearZ <- x$year # an indicator linking a visit to its parent year
  y <- y1 <-  numeric()
  for (j in 1:nrow(x)) {
    if(is.na(x$dominUnify[j])){ 
      y <- NA
          } else if (x$dominUnify[j]==0) { 
      y <- 0
          } else {
      y <- 1
          }
    y1 <- c(y1,y)
  }
  y <- y1 # visit-level detection history (binary)
  # Prepare data for JAGS
  Data <- list(N = N,
               Y = Y,
               n.Plot.pos = n.Plot.pos,
               cpos.Cens = cpos.Cens, # indicator (is censored?)
               cpos.Latent = cpos.Latent, # NA values for latent observations
               lims = lims,
               plot = plot,
               year = year,
               V2 = V2,
               plotZ = plotZ,
               yearZ = yearZ,
               y = y,
               yOrig = yOrig)
  zinit <- matrix(1, nrow = N, ncol = Y)
  inits.fn <- function() list(z = zinit,
                              tau.C = runif(1,1,5),
                              mu.C = rbeta(1,1,1),
                              mean.p = rbeta(1,1,1),
                              #gamma0 = rnorm(1,0,1),
                              gamma1 = runif(1,-5,5),
                              # cpos.Latent is approx. mid-points of the categories, used as initial values (dominUnif)
                              # intervals start from 1 (midpoint for the "zeroth" category not needed for these latent values for positive data)
                              cpos.Latent = c(0.001,0.025,0.04,0.075,0.175,0.29,0.375,0.625,0.85,0.975)[spPos$dominUnify] )
  #for ref only
  cPos.Init <- c(0.001,0.025,0.04,0.075,0.175,0.29,0.375,0.625,0.85,0.975)[spPos$dominUnify]
  ### MAKE SURE YOU HAVE THE RIGHT MODEL SCRIPT ###
  jagsModel <- rjags::jags.model(file= 'scripts/JAGS/JAGS_mod3.0.txt', data = Data, inits = inits.fn, n.chains = 3, n.adapt= 1000)
  ### MAKE SURE YOU HAVE THE RIGHT MODEL SCRIPT ###
  # Specify parameters for which posterior samples are saved
  #para.names <- c('psi')
  para.names <- c('mu.C', 'tau.C', 'gamma0', 'gamma1', 'annOcc', 'cPosAn')
  # Continue the MCMC runs with sampling
  samples <- rjags::coda.samples(jagsModel, variable.names = para.names, n.iter = 1500)
  ## Inspect results
  out <- summary(samples)
  mu_sd <- out$stat[,1:2] #make columns for mean and sd
  q <- out$quantile[,c(1,3,5)] #make columns for median and CI
  tableOut <- as.data.frame(cbind(mu_sd,q)) #make table
  ##################################
  ## add DIC or similar as an output
  ##################################
  return(list(x,tableOut,samples))
}


# 11.01.2019
### Low occupancy and low cover seems to result in annual average occupancy being always around 50%
### confirmed by simulations with avg low cover (0.1) and various levels of avg occupancy
### only with high cover values that occupancy can be reliably estimated?

# Function runModels_v1() will be applied across the list created in 4_extractData.R (sppDatList), minus excluded species
sppModels <- list()
#sppModels <- lapply(seq_along(spForMods[1]), function(i) runModels_v1(i)) # test
sppModels <- lapply(seq_along(spForMods), function(i) runModels_v1(i))
names(sppModels) <- names(spForMods)
save(sppModels, file = "outputs/sppRuns/grasslands_26032019_v2_gamma0only.Rdata") # this is the run across all models that used the reduced set of plots (excluding plots with only NAs)

names(sppModels) <- names(spForMods[1])
#mean(test[grep(rownames(test), pattern = "cPosAn") & regexpr(text = rownames(test), pattern = "(\\d+)\\D*\\z", perl = T),3])
test <- sppModels[[2]] # summarised model results now second item in each species' list
test[grep(rownames(test), pattern = 'annOcc'),]
mean(test[c(12:1002),1])
hist(test[c(12:1002),1], breaks = 1000)
mean(test[c(1003:1993),1])
hist(test[c(1003:1993),1], breaks = 1000)

######################################
## JAGS model
######################################
sink('scripts/JAGS/JAGS_mod3.0.txt')
cat("
model{
for (j in 1:Y){ # years
    # positive cover
    a[j] <- mu[j] * phi[j]
    b[j] <- (1 - mu[j]) * phi[j]
    cPos[j] ~ dbeta(a[j], b[j]) T(1e-16,0.9999999999999999) # random draw of cPos for that year, based on data-informed beta distribution
  for (i in 1:N){ # plots
      # ZI cover
      C[i,j] <- z[i,j] * cPos[j] # combine random draw of cPos for year j with estimated plot-specific occupancy
      
      # occupancy
      z[i,j] ~ dbern(psi[i,j]) # true PA state of a plot within a year depends on occupancy probability psi
      logit(psi[i,j]) <- m[j] # year random intercept on occupancy prob (same as Sparta)
      #psi[i,j] ~ dbeta(1,1)
      } # end of plots loop
} # end of years loop

## Derived values from state model
for (j in 1:Y){ # number of years
    annOcc[j] <- (sum(z[,j]))/N # annual occupancy
    mC[j] <- mean(C[,j]) # annual mean cover, including zeros
    mPsi[j] <- mean(psi[,j]) # annual mean occupancy prob
}
avgOcc <- mean(annOcc[]) # average annual occupancy

## Plot positive covers
for(k in 1:n.Plot.pos){ 
    cposCens[k] ~ dinterval(cposLatent[k], lims[k,])
    cposLatent[k] ~ dbeta(a[year[k]], b[year[k]]) T(1e-16,0.9999999999999999) # recorded cover when present follows beta distribution
}

## Observation model for all plot visits ([within-year] detection within plots)
for (a in 1:V2){
    x[a] ~ dbern(py[a]) # detectability influences detection
    py[a] <- z[plotZ[a], yearZ[a]] * p.dec[a] # true state x detectability
    p.dec[a] <- min(max(1e-16, p.Dec[a]), 0.9999999999999999) # trick to stop numerical problems
    
    # detectability model
    logit(p.Dec[a]) <- gamma0 + gamma1 * yOrig[a]
    yOrig[a] ~ dunif(0,10)
}

## Priors!
# Cover
phi[1] ~ dt(0, 0.01, 1)T(0,)
mu[1] ~ dbeta(1, 1)
for (j in 2:Y){
  mInt[j] ~ dnorm(logit(mu[j-1]), 4)
  mu[j] <- ilogit(mInt[j]) # mean follows low-variation random walk on logit-normal scale (sd = 0.5, so tau = 4)
  phi[j] ~ dt(0, 0.01, 1)T(0,) # phi's are independent half-Cauchy priors
}

# Detection model
mean.p ~ dbeta(1,1) # broad intercept on prob scale
gamma0 <- logit(mean.p) # transformed # note that this requires tweak to initial values
gamma1 ~ dnorm(0, 4) # informative prior to help with quasi-complete separation in logistic regression

# Intercept for state model for occupancy
mean.m ~ dbeta(1,1)
for(j in 1:Y){
  m[j] ~ dnorm(logit(mean.m), tau.m)
}  
tau.m <- 1/pow(sd.m, 2)
sd.m ~ dt(0, 0.001, 1)T(0,)

} ## END MODEL
", fill = TRUE)
sink()
