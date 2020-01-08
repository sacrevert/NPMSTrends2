### 7_PriorPredChecks
# O.L. Pescott
# 03.01.2020
rm(list=ls())
######################################
library(R2jags)
library(dplyr)
######################################

load(file = "data/Achi_mille_grassSamples_20190909.Rdata") ## pre-processed v. 09 09 2019
#load("W:/PYWELL_SHARED/Pywell Projects/BRC/_BRC_projects/NPMS/Analyses/2018 08 - Per species trend analyses/r_proj/NPMStrends/data/Achi_mille_grassSamples_20180920.Rdata")
head(Achi_mill_PAN); tail(Achi_mill_PAN); unique(Achi_mill_PAN$dominUnify)
# reorder grazing
unique(Achi_mill_PAN$grazing) # "High"     "Low"      "Moderate"
Achi_mill_PAN$grazing <- factor(Achi_mill_PAN$grazing, levels(Achi_mill_PAN$grazing)[c(2,3,1)])
levels(Achi_mill_PAN$grazing)
Achi_mill_PAN$grazing <- as.numeric(Achi_mill_PAN$grazing)
# simplify cutting/mowing text
levels(Achi_mill_PAN$cutting) <- "cut/mow"
# reorder survey levels
levels(Achi_mill_PAN$title)# "Indicator survey"  "Inventory survey"  "Wildflower survey"
Achi_mill_PAN$title <- factor(Achi_mill_PAN$title, levels(Achi_mill_PAN$title)[c(3,1,2)])
levels(Achi_mill_PAN$title)
Achi_mill_PAN$title <- as.numeric(Achi_mill_PAN$title)

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
#lims = lims # the limits to the intervals used for censoring 
#(again, these should be of one type only, unless data collected under different schemes are combined) -- see Pescott et al. (2016)

## 1. Achillea millefolium in grassland samples
###############################################
# Change 25.03.2019 -- delete plots for which
# samples are entirely NA (i.e. no known states)
###############################################
Achi_mill_PAN <- (Achi_mill_PAN %>% 
                    group_by(plot_id) %>%
                    filter(any(!is.na(dominUnify))) # resolves to true if there is at least one non-NA value within a group of plot samples
                    )

# not sure it really matters, but sorting the data by year and then plot first may simplify downstream things
Achi_mill_PAN$year <- format(Achi_mill_PAN$date.x, "%Y") # create year column
Achi_mill_PAN <- Achi_mill_PAN[order(Achi_mill_PAN$year, Achi_mill_PAN$plot_id), ]
uniPlots <- unique(Achi_mill_PAN$plot_id) # unique plot IDs - 31/12/2018 = 933
# create unique plot index (useful for cross-referencing later on?)
plotIndex <- data.frame(plot = uniPlots, index = 1:length(uniPlots))
N <- length(uniPlots) # number of unique plots
# check that all years are in the range of 2015:System time
if ( max(unique(format(Achi_mill_PAN$date.x, "%Y"))) > format(Sys.time(), "%Y") ) {
    print ("error")
} else {
    print("OK")
}
# 
Y <- length(unique(format(Achi_mill_PAN$date.x, "%Y"))) # number of years
# All Domin data, including zeros and NAs
yOrig <- Achi_mill_PAN$dominUnify; hist(yOrig)
# All grazing data
x2 <- Achi_mill_PAN$grazing
# All survey level data
x3 <- Achi_mill_PAN$title
# number of samples with presences
n.Plot.pos <- length(Achi_mill_PAN$dominUnify[Achi_mill_PAN$dominUnify !='0' & !is.na(Achi_mill_PAN$dominUnify)]) # 372
## checks on n.Plot.pos
head(Achi_mill_PAN[!is.na(Achi_mill_PAN$dominUnify) & Achi_mill_PAN$dominUnify !='0',]) # data frame
head(Achi_mill_PAN$dominUnify[Achi_mill_PAN$dominUnify !='0' & !is.na(Achi_mill_PAN$dominUnify)]) # vector
## looks OK

cpos.Cens <- rep(1, n.Plot.pos) # all samples censored
cpos.Latent <- rep(NA, n.Plot.pos) # all require estimation of latent state
t <- c(1e-4, 0.01, # cut points (note that these need to correspond to inits.fn and cPos.init below for JAGS model)
       0.01, 0.03,
       0.03, 0.05,
       0.05, 0.1,
       0.1, 0.25,
       0.25, 0.33,
       0.33, 0.5,
       0.5, 0.75,
       0.75, 0.95, 
       0.95, 0.9999)
t = matrix(t, nrow = 10, ncol = 2, byrow = TRUE)
tdf <- as.data.frame(t)
colnames(tdf) <- c('L','U') # 'L'ower and 'U'pper bounds of categories
tdf$int <- c(1,2,3,4,5,6,7,8,9,10); tdf

# positive rows only
spPos <- Achi_mill_PAN[Achi_mill_PAN$dominUnify !='0' & !is.na(Achi_mill_PAN$dominUnify),]; head(spPos)
# add row number to allow resorting after merge (merge with sort = F does not actually do what we want, could also used plyer::join)
spPos$indexPos <- 1:nrow(spPos)
spPos <- merge(spPos, tdf, by.x = "dominUnify", by.y = "int", all.x = T, all.y = F)
spPos <- spPos[order(spPos$indexPos),]
lims <- spPos[,c("L","U")] # has the lower/upper cutpoints for all intervals

## Required data 2
#plot = plot # an indicator linking a percentage cover observation to its parent (spatially unique) plot
#year = year # an indicator linking a percentage cover observation to its parent year
#V2 = V2 # the total number of visits (samples), irrespective of whether there is a postive cover for a species or not
#plotZ = plotZ # an indicator linking a visit to its parent (spatially unique) plot
#yearZ = yearZ # an indicator linking a visit to its parent year
#x = x # visit-level detection history (binary)

#plot <- match(spPos$plot_id, Achi_mill_PAN$plot_id) # an indicator linking a percentage cover observation to its parent (spatially unique) plot
spPos <- merge(spPos, plotIndex, by.x = "plot_id", by.y = "plot", all.x = T, all.y = F)
spPos <- spPos[order(spPos$indexPos),]
plot <- spPos$index
#yDat <- data.frame(index = 1:length(unique(format(Achi_mill_PAN$date.x, "%Y"))), year = unique(format(Achi_mill_PAN$date.x, "%Y"))) # not really needed
spPos$year <- as.factor(spPos$year)
levels(spPos$year) <- 1:length(unique(format(Achi_mill_PAN$date.x, "%Y")))
year <- spPos$year # an indicator linking a percentage cover observation to its parent year
## All visit information
V2 <- nrow(Achi_mill_PAN) # the total number of visits (samples), irrespective of whether there is a positive cover for a species or not
#plotZ <- match(Achi_mill_PAN$plot_id, uniPlots)  # an indicator linking a visit to its parent (spatially unique) plot
Achi_mill_PAN <- merge(Achi_mill_PAN, plotIndex, by.x = "plot_id", by.y = "plot", all.x = T, all.y = F)
Achi_mill_PAN <- Achi_mill_PAN[order(Achi_mill_PAN$year, Achi_mill_PAN$plot_id), ]
plotZ <- Achi_mill_PAN$index
Achi_mill_PAN$year <- as.factor(Achi_mill_PAN$year)
levels(Achi_mill_PAN$year) <- 1:length(unique(format(Achi_mill_PAN$date.x, "%Y")))
yearZ <- Achi_mill_PAN$year # an indicator linking a visit to its parent year
x <- x1 <-  numeric()
for (i in 1:nrow(Achi_mill_PAN)) {
  if(is.na(Achi_mill_PAN$dominUnify[i])){ 
    x <- NA
    #print(c(i,x)) # this was just for catching a bug
  } else if (Achi_mill_PAN$dominUnify[i]==0) { 
    x <- 0
    #print(c(i,x))
  } else {
    x <- 1
    #print(c(i,x))
  }
  x1 <- c(x1,x)
}
x <- x1 # visit-level detection history (binary)

####################
## Send data to JAGS
####################
# Data list for passing to JAGS
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
             x = x,
             #x2 = x2, # grazing (with NAs)
             #x2 = c(rep(x = 1:3, times = 788) , 1),
             x3 = x3, # survey level
             yOrig = yOrig,
             g2Levs = 3, 
             g3Levs = 3,
             pi = c(0.33,0.33,0.33)) 
###################################### END OF DATA PREP

###########################################
## Initialisation of values for JAGS chains
###########################################
# Initial parameter values for JAGS
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
                            cposLatent = c(0.001,0.025,0.04,0.075,0.175,0.29,0.375,0.625,0.85,0.975)[spPos$dominUnify]
)
#for ref only
cPosInit <- c(0.001,0.025,0.04,0.075,0.175,0.29,0.375,0.625,0.85,0.975)[spPos$dominUnify]

######################################
## JAGS model
######################################
sink('scripts/JAGS/JAGS_mod4.0.txt')
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
    
    # detectability model: g1 is cover effect, g2 is grazing effect (with missing values), g3 is survey level (no missing values)
    ##############################################################
    ## Currently breks here when you add in categorical predictor g3
    ##############################################################
    logit(p.Dec[a]) <- g0 + g1*yOrig[a] + g3[x3[a]] + g2[x2[a]] 
    yOrig[a] ~ dunif(0,10) # this could be improved using the Wilson/Irvine approach (i.e. imputing new values where missing, 
                           # rather than relying on noise-inducing uniform draws)
    x2[a] ~ dcat(pi) # missing values for grazing info (or assume 0)
}

### Priors ###
## Intercept for state model for occupancy
mean.m ~ dbeta(1,1)
for(j in 1:Y){
  m[j] ~ dnorm(logit(mean.m), tau.m)
}  
tau.m <- 1/pow(sd.mA, 2)
sd.mA <- sd.m + 0.1 # see Kruschke page 486 (don't want shrinkage to be too strong when data sparse)
sd.m ~ dt(0, 0.1, 1)T(0,)

## Cover (with random walk prior)
phi[1] ~ dt(0, 0.01, 1)T(0,)
mu[1] ~ dbeta(1, 1)
for (j in 2:Y){
  mInt[j] ~ dnorm(logit(mu[j-1]), 4)
  mu[j] <- ilogit(mInt[j]) # mean follows low-variation random walk on logit-normal scale (sd = 0.5, so tau = 4)
  phi[j] ~ dt(0, 0.1, 1)T(0,) # phi's are independent half-Cauchy priors (make tau = 0.1 given likely values of a/b)
}

## Detection model
mean.p ~ dbeta(1,1) # broad intercept on prob scale
g0 <- logit(mean.p) # transformed # note that this requires tweak to initial values
g1 ~ dunif(-5, 5)

# Hierarchical prior for levels of grazing
for (j in 1:g2Levs) { g2[j] ~ dnorm(g2mu, g2tau) }
mean.g2 ~ dbeta(1,1)
g2mu <- logit(mean.g2)
g2tau <- 1/pow(sd.g2A, 2)
sd.g2A <- sd.g2 + 0.1 # see Kruschke page 486 (don't want shrinkage to be too strong when data sparse)
sd.g2 ~ dt(0, 0.01, 1)T(0,)

# Hierarchical prior for survey levels
for (j in 1:g3Levs) { g3[j] ~ dnorm(g3mu, g3tau) }
mean.g3 ~ dbeta(1,1)
g3mu <- logit(mean.g3)
g3tau <- 1/pow(sd.g3A, 2)
sd.g3A <- sd.g3 + 0.1 # see Kruschke page 486 (don't want shrinkage to be too strong when data sparse)
sd.g3 ~ dt(0, 0.01, 1)T(0,)

} ## END MODEL
", fill = TRUE)
sink()

jagsModel <- jags.model(file= 'scripts/JAGS/JAGS_mod3.0.txt', data = Data, inits = inits.fn, n.chains = 3, n.adapt= 500)

# Specify parameters for which posterior samples are saved
para.names <- c('mC', 'mPsi', 'mu', 'annOcc', 'avgOcc')
para.names <- c('mC')
# Continue the MCMC runs with sampling
samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 500)
## Inspect results
out <- summary(samples)
mu_se <- out$stat[,c(1,3)] #make columns for mean and se
#save(samples, file = "outputs/")
gelman.diag(samples)
plot(samples)

#############
#### dt stuff
#############
# below is for dt(0, 0.01, 1)T(0,)
x <- rnorm(n=1000, 0, sd = sqrt(1/0.1)) # tau = 0.01
s <- rgamma(n=1000, shape = 1/2, rate = 1/2) # k = 1; shape = rate = k/2 = 1/2
y <- x/sqrt(s)
hist(abs(y), breaks = 10000, xlim = c(0,250))

# cf. to used in Kruschke (much lower density around0-50 and much longer RHS tail)
x <- rnorm(n=1000, 0, sd = sqrt(1/0.001)) # tau = 0.01
s <- rgamma(n=1000, shape = 2/2, rate = 2/2) # k = 1; shape = rate = k/2 = 1/2
y <- x/sqrt(s)
hist(abs(y), breaks = 10000, xlim = c(0,250))
