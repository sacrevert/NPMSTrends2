### 5c_processDataForTrends_Example2_removeNAPlots.R
# now works for Achillea millefolium example
# O.L. Pescott
# 17.09.2018, fixed 31.12.2018
# 09.01.2019, annual derived values and Domin covar on detectability added
# 25.03.2019, remove plots for which only NA values have been recorded (could change in future to leave in ones within monads where species has known P/A?)
# 09.09.2019, updates to accommodate new data for detection model, plus new model for estimates of annual cover per species.
rm(list=ls())
######################################
library(R2jags)
library(dplyr)
######################################

#load(file = "data/Achi_mille_grassSamples_2020-01-16.Rdata") ## pre-processed
load(file = "data/Gymn_conop_grassSamples_2020-01-16.Rdata") ## pre-processed
#load(file = "data/Trif_camp_grassSamples_2020-01-24.Rdata") ## pre-processed
#dat <- Trif_camp_PAN
dat <- Gym_con_PAN
str(dat)
#dat$date.x <- as.Date(as.character(dat$date.x), format = "%d/%m/%Y")
#load("W:/PYWELL_SHARED/Pywell Projects/BRC/_BRC_projects/NPMS/Analyses/2018 08 - Per species trend analyses/r_proj/NPMStrends/data/Achi_mille_grassSamples_20180920.Rdata")
head(dat); tail(dat); unique(dat$dominUnify)
# reorder grazing
dat$grazing[is.na(dat$grazing)] <- "Absent"
dat$grazing <- as.factor(dat$grazing)
levels(dat$grazing) # "Absent" "High" "Low" "Moderate"
dat$grazing <- factor(dat$grazing, levels(dat$grazing)[c(1,3,4,2)]) # reorder
levels(dat$grazing)
dat$grazing <- as.numeric(dat$grazing)
# simplify cutting/mowing text
levels(dat$cutting) <- "cut/mow"
# reorder survey levels
dat$title <- as.factor(dat$title)
levels(dat$title)# "Indicator survey"  "Inventory survey"  "Wildflower survey"
dat$title <- factor(dat$title, levels(dat$title)[c(3,1,2)])
levels(dat$title)
dat$title <- as.numeric(dat$title)

# domin things
#domins <- read.csv(file = "data/dominScores.csv", header = T, stringsAsFactors = F)
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
#library(dplyr)
## This bit was only necessary when we were using the wrong prior for the beta distribution on cover
## Presence of totally NA samples should just increase uncertainty, because for those we are sampling from the priors (for cover)
## Detectability, however, influenced by regression parameter estimates
#dat <- (dat %>% 
#            group_by(plot_id) %>%
#            filter(any(!is.na(dominUnify))) # resolves to true if there is at least one non-NA value within a group of plot samples
#            )


# not sure it really matters, but sorting the data by year and then plot first may simplify downstream things
dat$year <- format(dat$date.x, "%Y") # create year column
dat <- dat[order(dat$year, dat$plot_id), ]
uniPlots <- unique(dat$plot_id) # unique plot IDs - 31/12/2018 = 933
# create unique plot index (useful for cross-referencing later on?)
plotIndex <- data.frame(plot = uniPlots, index = 1:length(uniPlots))
N <- length(uniPlots) # number of unique plots
# check that all years are in the range of 2015:System time
if ( max(unique(dat$year)) > format(Sys.time(), "%Y") ) {
    print ("error")
} else {
    print("OK")
}
# 
Y <- length(unique(format(as.Date(dat$date.x), "%Y"))) # number of years
# All Domin data, including zeros and NAs
#hist(yOrig, breaks = 10) # includes value of 999, so doesn't make sense now
yOrig <- dat$dominUnify
# All grazing data
x2 <- dat$grazing
# All survey level data
x3 <- dat$title
# number of samples with presences
n.Plot.pos <- length(dat$dominUnify[dat$dominUnify !='0' & !is.na(dat$dominUnify)]) # 372
## checks on n.Plot.pos
head(dat[!is.na(dat$dominUnify) & dat$dominUnify !='0',]) # data frame
head(dat$dominUnify[dat$dominUnify !='0' & !is.na(dat$dominUnify)]) # vector
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
       0.95, 0.9999,
       1e-4,0.9999) # row for presences without cover value
t = matrix(t, nrow = 11, ncol = 2, byrow = TRUE)
tdf <- as.data.frame(t)
colnames(tdf) <- c('L','U') # 'L'ower and 'U'pper bounds of categories
tdf$int <- c(1,2,3,4,5,6,7,8,9,10,11); tdf

# positive rows only
spPos <- dat[dat$dominUnify !='0' & !is.na(dat$dominUnify),]; head(spPos)
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

#plot <- match(spPos$plot_id, dat$plot_id) # an indicator linking a percentage cover observation to its parent (spatially unique) plot
spPos <- merge(spPos, plotIndex, by.x = "plot_id", by.y = "plot", all.x = T, all.y = F)
spPos <- spPos[order(spPos$indexPos),]
plot <- spPos$index
#yDat <- data.frame(index = 1:length(unique(format(dat$date.x, "%Y"))), year = unique(format(dat$date.x, "%Y"))) # not really needed
spPos$year <- as.factor(spPos$year)
levels(spPos$year) <- 1:length(unique(format(dat$date.x, "%Y")))
year <- spPos$year # an indicator linking a percentage cover observation to its parent year
## All visit information
V2 <- nrow(dat) # the total number of visits (samples), irrespective of whether there is a positive cover for a species or not
#plotZ <- match(dat$plot_id, uniPlots)  # an indicator linking a visit to its parent (spatially unique) plot
dat <- merge(dat, plotIndex, by.x = "plot_id", by.y = "plot", all.x = T, all.y = F)
dat <- dat[order(dat$year, dat$plot_id), ]
plotZ <- dat$index
dat$year <- as.factor(dat$year)
levels(dat$year) <- 1:length(unique(format(dat$date.x, "%Y")))
yearZ <- dat$year # an indicator linking a visit to its parent year
y <- y1 <-  numeric()
for (i in 1:nrow(dat)) {
  if(is.na(dat$dominUnify[i])){ 
    y<- NA
    #print(c(i,x)) # this was just for catching a bug
  } else if (dat$dominUnify[i]==0) { 
    y <- 0
    #print(c(i,x))
  } else {
    y <- 1
    #print(c(i,x))
  }
  y1 <- c(y1,y)
}
y <- y1 # visit-level detection history (binary)

####################
## Send data to JAGS
####################
# Data list for passing to JAGS
yOrig[yOrig==11] <- NA
yOrig[is.na(yOrig)] <- 0
y.ind <- yOrig
y.ind[y.ind >= 1] <- 1
y.ind[is.na(y.ind)] <- 0
Data <- list(N = N,
             Y = Y,
             n.Plot.pos = n.Plot.pos,
             cposCens = cpos.Cens, # indicator (is censored?)
             cposLatent = cpos.Latent, # NA values for latent observations
             lims = lims,
             cuts = c(1e-16,1e-4,0.01,0.03,0.05,0.1,0.25,0.33,0.5,0.75,0.95,0.9999),
             #plot = plot,
             year = year,
             V2 = V2,
             plotZ = plotZ,
             yearZ = yearZ,
             y = y,
             x2 = x2, # grazing
             #x2 = c(rep(x = 1:3, times = 788) , 1),
             x3 = x3, # survey level
             yOrig = yOrig, # make sure yOrig = 11 is deleted (used as covar on detection)
             y.ind = y.ind,
             g2Levs = 4, 
             g3Levs = 3) 
###################################### END OF DATA PREP

###########################################
## Initialisation of values for JAGS chains
###########################################
# Initial parameter values for JAGS
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
                            # note extra value for dominUnify = 11
                            cposLatent = c(0.001,0.025,0.04,0.075,0.175,0.29,0.375,0.625,0.85,0.975,0.5)[spPos$dominUnify]
)
#for ref only
cPosInit <- c(0.001,0.025,0.04,0.075,0.175,0.29,0.375,0.625,0.85,0.975,0.5)[spPos$dominUnify]

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

jagsModel <- jags.model(file= 'scripts/JAGS/JAGS_mod4.3.txt', data = Data, inits = inits.fn, n.chains = 6, n.adapt= 100)

# Specify parameters for which posterior samples are saved
para.names <- c('mC', 'mPsi', 'mu', 'annOcc', 'avgOcc', 'g0', 'g1', 'g2', 'g3')
#para.names <- c('mC')
# Continue the MCMC runs with sampling
samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 10000, thin = 10)
save(samples, file = "samples.Rdata")
## Inspect results
out <- summary(samples)
mu_se <- out$stat[,c(1,3)] #make columns for mean and se
#save(samples, file = "outputs/")
gelman.diag(samples)
plot(samples)

### Visualise ###
#library(tidybayes)
library(MCMCvis)
library(mcmcplots)
MCMCtrace(samples, 
          params = c('mu'),
          ISB = FALSE,
          iter = 500,
          ind = TRUE,
          pdf = FALSE)
mcmcplot(samples)
caterplot(mcmcout = samples, parms = 'cPos')
denplot(mcmcout = samples, parms = 'mu')
traplot(samples, parms = 'mu')
parcorplot(samples, parms = c('cPos', 'mu'))


#############
#### dt stuff -- now in script 7: 7_PriorPredChecks.
#############
