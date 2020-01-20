################################################################################
#### 5b. Create function to make model runs across multiple species easier.#####
################################################################################
# O.L. Pescott, olipes@ceh.ac.uk
# 09.01.2019
rm(list=ls())
# v0 -- original version
# v1 -- amend function to remove all plots with no known states at time of processing
######################################
list.of.packages <- c("R2jags", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(R2jags) ## obviously requires that JAGS (preferably 4.3.0) is installed
library(dplyr)

# files for development
load(file = "outputs/grasslandsEg_2020-01-06.Rdata") # load saved files from 4_extractData.R

## Helper functions
ilt <- function(x){ # inverse logit function
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
# head(spForMods[[1]]); names(spForMods[1])
# x <- spForMods[[1]]
## Function for preparing species datasets for the JAGS model and then running model
runModels_v1 <- function(i, file = 'scripts/JAGS/JAGS_mod4.1.txt') {
  x <- spForMods[[i]]
  #x <- (x %>% group_by(plot_id) %>% filter(any(!is.na(dominUnify))) %>% as.data.frame()) # added in v1: keep plots with at least one non NA only
  if ( is.factor(x$date.x) ) { x$date.x <- as.Date(as.character(x$date.x), format = "%d/%m/%Y") }
  x$grazing <- as.numeric(factor(x$grazing, levels(x$grazing)[c(2,3,1)]))
  #levels(x$cutting) <- "cut/mow" # not needed yet!
  x$title <- factor(x$title, levels(x$title)[c(3,1,2)])
  x$year <- format(x$date.x, "%Y")
  x <- x[order(x$year, x$plot_id), ]
  uniPlots <- unique(x$plot_id) # unique plot IDs - 31/12/2018 = 933
  # create unique plot index (useful for cross-referencing later on?)
  plotIndex <- data.frame(plot = uniPlots, index = 1:length(uniPlots))
  N <- length(uniPlots)
  # check that all years are in the range of 2015:System time
  if ( max(unique(format(x$date.x, "%Y"))) > format(Sys.time(), "%Y") ) {
    print(paste("Error. Year violation: ", names(spForMods)[i]))
  } else {
    print(paste("Years ok: ", names(spForMods)[i]))
  }
  Y <- length(unique(format(x$date.x, "%Y")))
  # All Domin data, including zeros
  yOrig <- x$dominUnify
  x2 <- x$grazing
  x3 <- x$title
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
         0.95, 0.9999,
         1e-4,0.9999) #this row is for presences with unknown cover (given the value 11 in the Domin scale)
  tdf <- as.data.frame(matrix(t, nrow = 11, ncol = 2, byrow = TRUE))
  colnames(tdf) <- c('L','U') # 'L'ower and 'U'pper bounds of categories
  tdf$int <- c(1,2,3,4,5,6,7,8,9,10,11)
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
  yOrig[yOrig=11] <- NA # remove the value 11 (not appropriate as covar in detection model)
  Data <- list(N = N,
               Y = Y,
               n.Plot.pos = n.Plot.pos,
               cposCens = cpos.Cens, # indicator (is censored?)
               cposLatent = cpos.Latent, # NA values for latent observations
               lims = lims,
               #plot = plot,
               year = year,
               V2 = V2,
               plotZ = plotZ,
               yearZ = yearZ,
               y = y,
               x2 = x2,
               x3 = x3,
               yOrig = yOrig,
               g2Levs = 3,
               g3Levs = 3,
               pi = c(0.33,0.33,0.33))
  zinit <- matrix(1, nrow = N, ncol = Y)
  inits.fn <- function() list(z = zinit,
                              sd.m = runif(1,0,10),
                              sd.g2 = runif(1,0,10),
                              sd.g3 = runif(1,0,10),
                              mean.m = runif(1,0,1),
                              mean.p = runif(1,0,1),
                              mean.g2 = runif(1,0,1),
                              mean.g3 = runif(1,0,1),
                              phi = rep(runif(1,1,10), Y),
                              g1 = runif(1,-5,5),
                              cposLatent = c(0.001,0.025,0.04,0.075,0.175,0.29,0.375,0.625,0.85,0.975,0.5)[spPos$dominUnify] 
                              )
  #for ref only
  cPos.Init <- c(0.001,0.025,0.04,0.075,0.175,0.29,0.375,0.625,0.85,0.975,0.5)[spPos$dominUnify]
  ### MAKE SURE YOU HAVE GIVEN THE RIGHT MODEL SCRIPT TO THE FUNCTION ###
  jagsModel <- rjags::jags.model(file = file, data = Data, inits = inits.fn, n.chains = 3, n.adapt= 500)
  ### MAKE SURE YOU HAVE THE RIGHT MODEL SCRIPT ###
  # Specify parameters for which posterior samples are saved
  #para.names <- c('psi')
  para.names <- c('mC', 'mPsi', 'mu', 'annOcc', 'avgOcc')
  # Continue the MCMC runs with sampling
  samples <- rjags::coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)
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
sppModels <- lapply(seq_along(spForMods), function(i, file = 'scripts/JAGS/JAGS_mod4.1.txt') runModels_v1(i))
names(sppModels) <- names(spForMods)
save(sppModels, file = "outputs/sppRuns/grasslandResults_06012020_FULL.Rdata") # this is the run across all models that used the reduced set of plots (excluding plots with only NAs)

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
sink('scripts/JAGS/JAGS_mod4.1.txt')
cat("
model{
for (j in 1:Y){ # years
    # positive cover
    a[j] <- mu[j] * phi[j]
    b[j] <- (1 - mu[j]) * phi[j]
    cPos[j] ~ dbeta(a[j], b[j])T(1e-4,0.9999) # random draw of cPos for that year, based on data-informed beta distribution
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
    cposLatent[k] ~ dbeta(a[year[k]], b[year[k]])T(1e-4,0.9999) # recorded cover when present follows beta distribution
}

## Observation model for all plot visits ([within-year] detection within plots)
for (a in 1:V2){
    y[a] ~ dbern(py[a]) # detectability influences detection
    py[a] <- z[plotZ[a], yearZ[a]] * p.dec[a] # true state x detectability
    p.dec[a] <- min(max(1e-4, p.Dec[a]), 0.9999) # trick to stop numerical problems
    
    # detectability model: g1 is cover effect, g2 is grazing effect (with missing values), g3 is survey level (no missing values)
    ##############################################################
    ## Currently breaks here when you add in categorical predictor g3
    ##############################################################
    logit(p.Dec[a]) <- g0 + g1*yOrig[a] + g2[x2[a]] + g3[x3[a]]
    yOrig[a] ~ dunif(0,10) # this could be improved using the Wilson/Irvine approach (i.e. imputing new values where missing, 
                           # rather than relying on noise-inducing uniform draws)
    x2[a] ~ dcat(pi) # missing values for grazing info (or assume 0): any equally likely (could supplement with veg height info)
    #x3[a] ~ dcat(pi) # missing values for survey level. Only needed for PPP as there are no missing values for this covar in dataset
}

### Priors ###
## Intercept for state model for occupancy
mean.m ~ dbeta(1,1)T(1e-4,0.9999)
for(j in 1:Y){
  m[j] ~ dnorm(logit(mean.m), tau.m)
}  
tau.m <- 1/pow(sd.mA, 2)
# Might want to check the following position (sd.m + 0.1) through PPP
sd.mA <- sd.m #+ 0.1 # see Kruschke page 486 (don't want shrinkage to be too strong when data sparse)
sd.m ~ dt(0, 0.1, 1)T(0,)

## Cover (with random walk prior)
#phi[1] ~ dt(0, 0.01, 1)T(0,)
phi[1] ~ dpar(1.5, 0.1) # alpha, c
mu[1] ~ dbeta(1, 1)T(1e-4,0.9999)
for (j in 2:Y){
  mInt[j] ~ dnorm(logit(mu[j-1]), 4)
  mu[j] <- ilogit(mInt[j]) # mean follows low-variation random walk on logit-normal scale (sd = 0.5, so tau = 4)
  #phi[j] ~ dt(0, 0.1, 1)T(0,) # phi's are independent half-Cauchy priors (make tau = 0.1 given likely values of a/b)
  phi[j] ~ dpar(1.5, 0.1) # alpha, c
}

## Detection model
mean.p ~ dbeta(1,1)T(1e-4,0.9999) # broad intercept on prob scale
g0 <- logit(mean.p) # transformed # note that this requires tweak to initial values
g1 ~ dunif(-5, 5)

# Hierarchical prior for levels of grazing
for (j in 1:g2Levs) { g2[j] ~ dnorm(g2mu, g2tau) }
mean.g2 ~ dbeta(1,1)T(1e-4,0.9999)
g2mu <- logit(mean.g2)
g2tau <- 1/pow(sd.g2A, 2)
#sd.g2A <- sd.g2 + 0.1 # see Kruschke page 486 (don't want shrinkage to be too strong when data sparse)
sd.g2A <- sd.g2
sd.g2 ~ dt(0, 0.01, 1)T(0,)

# Hierarchical prior for survey levels
for (j in 1:g3Levs) { g3[j] ~ dnorm(g3mu, g3tau) }
mean.g3 ~ dbeta(1,1)T(1e-4,0.9999)
g3mu <- logit(mean.g3)
g3tau <- 1/pow(sd.g3A, 2)
#sd.g3A <- sd.g3 + 0.1 # see Kruschke page 486 (don't want shrinkage to be too strong when data sparse)
sd.g3A <- sd.g3
sd.g3 ~ dt(0, 0.01, 1)T(0,)

} ## END MODEL
", fill = TRUE)
sink()
