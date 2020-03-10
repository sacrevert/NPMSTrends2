################################################################################
#### 5c. Create function to make model runs across multiple species easier.#####
          ################   CLUSTER VERSION ################   
################################################################################
# O.L. Pescott, olipes@ceh.ac.uk
# 10.03.2020
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
runModels_v1 <- function(i, file = 'scripts/JAGS/JAGS_mod4.3.txt') {
  x <- spForMods[[i]]
  #x <- (x %>% group_by(plot_id) %>% filter(any(!is.na(dominUnify))) %>% as.data.frame()) # added in v1: keep plots with at least one non NA only
  if ( is.factor(x$date.x) ) { x$date.x <- as.Date(as.character(x$date.x), format = "%d/%m/%Y") }
  if ( is.character(x$grazing) ) {x$grazing[is.na(x$grazing)] <- "Absent" }
  if ( is.factor(x$grazing) ) {x$grazing <- as.character(x$grazing); x$grazing[is.na(x$grazing)] <- "Absent"}
  x$grazing <- factor(x$grazing, ordered = T, levels = c("Absent", "Low", "Moderate", "High"))
  #levels(x$cutting) <- "cut/mow" # not needed yet!
  x$title <- factor(x$title, ordered = T, levels = c("Wildflower survey", "Indicator survey", "Inventory survey"))
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
  yOrig[yOrig==11] <- NA # remove these for the purposes of using yOrig as a covar in detection model
  yOrig[is.na(yOrig)] <- 0 # no cover = zero for covar purposes
  y.ind <- yOrig # also use Orig to create indicator to indicate whether there is a cover value or one needs creating
  y.ind[y.ind >= 1] <- 1 # yOrig cover value present
  y.ind[is.na(y.ind)] <- 0 # yOrig cover value present
  Data <- list(N = N,
               Y = Y,
               n.Plot.pos = n.Plot.pos,
               cposCens = cpos.Cens, # indicator (is censored?)
               cposLatent = cpos.Latent, # NA values for latent observations
               lims = lims,
               cuts = c(1e-16,1e-4,0.01,0.03,0.05,0.1,0.25,0.33,0.5,0.75,0.95,0.9999),
               year = year,
               V2 = V2,
               plotZ = plotZ,
               yearZ = yearZ,
               y = y,
               x2 = as.numeric(x2), # grazing
               x3 = as.numeric(x3), # survey level
               yOrig = yOrig, # make sure yOrig = 11 is deleted (used as covar on detection)
               y.ind = y.ind,
               g2Levs = length(unique(x2)), # sometimes some levels might be absent
               g3Levs = length(unique(x3))) # sometimes some levels might be absent
  zinit <- matrix(1, nrow = N, ncol = Y)
  inits.fn <- function() list(z = zinit,
                              sdC = runif(1,0,10),
                              sd0 = runif(1,0,10),
                              sd.m = runif(1,0,10),
                              sd.g2 = runif(1,0,10),
                              sd.g3 = runif(1,0,10),
                              mean.m = runif(1,0,1),
                              mean.p = runif(1,0,1),
                              mean.g2 = runif(1,0,1),
                              mean.g3 = runif(1,0,1),
                              phi = rep(runif(1,1,10), Y),
                              y.hat = rep(runif(1,0,1), V2),
                              g1 = runif(1,0,3),
                              # note extra initial value for dominUnify = 11 (represents a presence with unknown cover)
                              cposLatent = c(0.001,0.025,0.04,0.075,0.175,0.29,0.375,0.625,0.85,0.975,0.5)[spPos$dominUnify]
  )
  #for ref only
  cPos.Init <- c(0.001,0.025,0.04,0.075,0.175,0.29,0.375,0.625,0.85,0.975,0.5)[spPos$dominUnify]
  ### MAKE SURE YOU HAVE GIVEN THE RIGHT MODEL SCRIPT TO THE FUNCTION ###
  jagsModel <- rjags::jags.model(file = file, data = Data, inits = inits.fn, n.chains = 6, n.adapt = 100)
  ### MAKE SURE YOU HAVE THE RIGHT MODEL SCRIPT ###
  # Specify parameters for which posterior samples are saved
  #para.names <- c('psi')
  para.names <- c('mC', 'mPsi', 'mu', 'annOcc', 'avgOcc')
  # Continue the MCMC runs with 
  update(jagsModel, n.iter = 5000) # burn in!
  samples <- rjags::coda.samples(jagsModel, variable.names = para.names, n.iter = 5000, thin = 5)
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
sppModels <- lapply(seq_along(spForMods[i]), function(i, file = 'scripts/JAGS/JAGS_mod4.3.txt') runModels_v1(i))
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
sink('scripts/JAGS/JAGS_mod4.3.txt')
cat("
model{
for (j in 1:Y){ # years
    # positive cover
    alpha[j] <- mu[j] * phi[j]
    beta[j] <- (1 - mu[j]) * phi[j]
    cPos[j] ~ dbeta(alpha[j], beta[j]) T(1e-4,0.9999) # random draw of cPos for that year, based on data-informed beta distribution
  for (i in 1:N){ # plots
      # ZI cover
      C[i,j] <- z[i,j] * cPos[j] # combine random draw of cPos for year j with estimated plot-specific occupancy
      
      # occupancy
      z[i,j] ~ dbern(psi[i,j]) # true PA state of a plot within a year depends on occupancy probability psi
      logit(psi[i,j]) <- m[j] # year random intercept on occupancy prob
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
    cposLatent[k] ~ dbeta(alpha[year[k]], beta[year[k]])T(1e-4,0.9999) # recorded cover when present follows beta distribution
}

## Observation model for all plot visits ([within-year] detection within plots)
for (a in 1:V2){
    y[a] ~ dbern(py[a]) # detectability influences detection
    py[a] <- z[plotZ[a], yearZ[a]] * p.dec[a] # true state x detectability
    p.dec[a] <- min(max(1e-4, p.Dec[a]), 0.9999) # trick to stop numerical problems
    
    # detectability model: g1 is cover effect, g2 is grazing effect (with missing values), g3 is survey level (no missing values)
    logit(p.Dec[a]) <- g0 + g1*y.cov[a] + g3[x3[a]] + g2[x2[a]] 
    y.cov[a] <- yOrig[a]*y.ind[a] + y.int[a]*(1-y.ind[a]) # y.cov is a combination of original Domin scores, and predicted Domin scores where NA
    y.int[a] <- dinterval(y.hat[a], cuts)
    y.hat[a] ~ dbeta(alpha[yearZ[a]],beta[yearZ[a]])T(1e-4,0.9999)
}

### Priors ###
## Intercept for state model for occupancy
mean.m ~ dbeta(1,1)T(1e-4,0.9999)
m[1] ~ dnorm(logit(mean.m), tau0)
tau0 <- 1/pow(sd0, 2)
sd0 ~ dt(0, 1, 3)T(0,)
for(j in 2:Y){
  m[j] ~ dnorm(m[j-1], tau.m)T(logit(1e-4), logit(0.9999)) # fixing tau.m (e.g. at 4) results in much smoother trend
}  
tau.m <- 1/pow(sd.m, 2)
sd.m ~ dt(0, 1, 3)T(0,)

## Cover (with random walk prior)
phi[1] ~ dpar(1.5, 0.1) # alpha, c
mu[1] ~ dbeta(1, 1)T(1e-4,0.9999)
for (j in 2:Y){
  mu[j] <- ilogit(mInt[j]) # mean follows low-variation random walk on logit-normal scale (sd = 0.5, so tau = 4)
  mInt[j] ~ dnorm(logit(mu[j-1]), tauC)T(logit(1e-4), logit(0.9999)) # fixing tauC (e.g. at 4) results in much smoother trend
  phi[j] ~ dpar(1.5, 0.1) # alpha, c
}
tauC <- 1/pow(sdC, 2)
sdC ~ dt(0, 1, 3)T(0,) #-- weaker, density mostly -10 - 10
#sdC ~ dt(0, 0.001, 1)T(0,) # bit too strong (see 7_PriorPredChecks.R script) -- normally most density -5 to 5

## Detection model
mean.p ~ dbeta(1,1)T(1e-4,0.9999) # broad intercept on prob scale
g0 <- logit(mean.p) # transformed # note that this requires tweak to initial values
g1 ~ dnorm(0, 0.1)T(0,) # moderate belief that this parameter should be positive, increase in cover increases detectability, all else being equal

# Hierarchical prior for levels of grazing
for (j in 1:g2Levs) { g2[j] ~ dnorm(g2mu, g2tau) }
mean.g2 ~ dbeta(1,1)T(1e-3,0.999) # between -7 and 7
g2mu <- logit(mean.g2)
g2tau <- 1/pow(sd.g2, 2)
sd.g2 ~ dt(0, 1, 3)T(0,)

# Hierarchical prior for survey levels
for (j in 1:g3Levs) { g3[j] ~ dnorm(g3mu, g3tau) }
mean.g3 ~ dbeta(1,1)T(1e-3,0.999) # between 
g3mu <- logit(mean.g3)
g3tau <- 1/pow(sd.g3, 2)
sd.g3 ~ dt(0, 1, 3)T(0,) # dt(mu,tau,k)

} ## END MODEL
", fill = TRUE)
sink()
