## 6. Simulation code
#
# O.L. Pescott
# 29.08.2019
rm(list=ls())

######################################
library(R2jags)
library(mcmcplots)
#library(BayesianTools)
######################################
set.seed(4323) # for reproducibility
######################################
## Define potentially useful functions
ilt<-function(x){ # inverse logit function
  exp(x)/(1+exp(x))
}
logit <- function(x){ # logit function
  log(x/(1 - x))
}
######################################

######################################
## Could put this in a simulation function (see Wright et al. 2017)
N <- 100 # number of spatially unique plots
J <- 2 # number of visits to a plot within a year (assume constant for the moment)
psi <- 1 # should be easier for model to retrieve true values of gamma0 and gamma1 if occupancy is 1 -- this appears to be true
Y <- 2 # total number of years of monitoring covered by the data
mu0 <- 0.50 # parameter for mean of cover beta distribution # 0.25
tau0 <-3 # parameter for 'precision' of cover distribution # 3
gamma0 <- -1.5 # intercept for detection logistic regression # -1.5
gamma1 <- 2 # slope for detection with %cover # 2
# e.g. plogis(-2 + 3*cover) makes for greater detectability range based on percent covers
# if you do this but only estimate gamma0 in the model, then biases in gamma0 and mu.C estimates increase
######################################

# array of plot covers per visit per year
y.array <- array(dim = c(N, J, Y)) # plots, visits, years
for(k in 1:Y){
      y.array[,,k] <- matrix(rbeta(N*J, mu0*tau0, (1-mu0)*tau0), 
                                 nrow=N, ncol=J)
    for (i in 1:N){
      y.array[i,,k] <- rep(rbinom(1,1,psi), J) * y.array[i,,k] # plot level occupancy within years
    }
}
# make a detection history matrix based on cover data
x.array <- array(dim = c(N, J, Y))
for(k in 1:Y){
  for(i in 1:N){
    for(j in 1:J){
      x.array[i, j, k] <- ifelse(y.array[i, j, k] > 0,
                                 rbinom(1, 1, 
                                       #plogis(gamma0 + gamma1*y.array[i, j, k])), # standard detection function with link between cover and detection
                                        #plogis(gamma0)), # intercept only
                                 1), # for perfect detection
                                 0)
    }
  }
}
y.arrayOrig <- y.array # keep original covers
y.array[which(x.array==0)] <- 0 # if the plant was not actually detected, then set the recorded cover to zero as well
y.array[which(y.array<1e-16 & y.array > 0)] <- 1e-16 # put hard limit on lower cover bound
###################################### END OF SIMS

##################################################
## Data/indicators required for running JAGS model
##################################################
# total number of visits with positive covers
cpos <- matrix(NA, nrow = length(as.vector(y.array)[which(as.vector(y.array) > 0)]), ncol = 3) # empty matrix for actual covers and cover classes
cpos[,1] <- as.vector(y.array)[which(as.vector(y.array) > 0)] # actual cover value for every positive cover
cpos[,2] <- 1 # indicator stating the observation censored
cpos[,3] <- NA # NA values for latent variable
# code borrowed from Pescott et al. 2016 Power paper:
t <- c(1e-16,0.05, # need to tweak lower categories if starting with low percentage cover like 0.01
       0.05,0.25,
       0.25,0.5,
       0.5,0.75,
       0.75,0.95,
       0.95,0.9999999999999999)
t = matrix(t, nrow = 6, ncol = 2, byrow = TRUE)
tdf <- as.data.frame(t)
colnames(tdf) <- c('L','U') # 'L'ower and 'U'pper bounds of categories
intervals <- c(1,2,3,4,5,6)
tdf$int <- intervals
# library(plyr) for SQL join function (like merge but keeps order) -- although actually could use merge with sort = F
tInt <- c(1e-16,0.05,0.25,0.5,0.75,0.95,0.9999999999999999) # need to tweak lower categories if starting with low percentage cover like 0.01
int <- findInterval(cpos[,1], tInt) # find corresponding interval for all data points
int <- as.data.frame(int)
m1 <- plyr::join(int, tdf, by = "int", type = "left", match = "all") # change to using merge at some point (although merge resorts, need to add row index to resort on)
lims <- m1[,2:3] # has the intervals for all points
check <- cbind(cpos, m1); head(check); tail(check) # just for quick visual check that all is well
cposCens <- rep(1, nrow(cpos))
cposLatent <- cpos[,3] # just NAs

n.Plot.pos <- nrow(cpos)
# indicators linking positive plot k to the correct plots and years; length = nrow(cpos)
out <- out1 <-  numeric()
for(k in 1:Y){
  for(j in 1:J){
    for(i in 1:N){
    out <- ifelse(y.array[i, j, k] > 0, i, NA)
    out1 <- c(out1, out)
    }
  }
}
plot <- out1[!is.na(out1)]

out <- out1 <-  numeric()
for(k in 1:Y){
  for(j in 1:J){
    for(i in 1:N){
    out <- ifelse(y.array[i, j, k] > 0, k, NA)
    out1 <- c(out1, out)
    }
  }
}
year <- out1[!is.na(out1)]

#############
### NEEDS CHANGING WHEN THERE ARE UNEVEN VISITS NUMBERs ETC. BETWEEN YEARS (something to consider for real data)
#############
V2 <- N*J*Y # total number of plot visits (# plots x # visits x # years)
V3 <- N*Y
# indicator linking every visit x to plot; length = V2 
plotZ <- rep(1:N, J*Y)
# indicator linking every visit x to year; length V2 
yearZ <- rep(1:Y, each = N*J)
x <- as.vector(x.array) # detection indicator for every visit (length V2)
# covers for all visits
#y <- as.vector(y.arrayOrig) # original covers used in detectability loop (inc. zeros) for every visit (length V2)
y <- as.vector(y.array) # covers after non-detections -- shouldn't use previous line, as this information is not available in reality
x3 <- sample(c(1:3), 400, replace = T)

# Data list for passing to JAGS
Data <- list(N = N,
            Y = Y,
            n.Plot.pos = n.Plot.pos,
            #cpos.Cens = cpos[,2], # indicator (is censored?)
            cposCens = cposCens, # indicator (is censored?)
            cposLatent = cposLatent, # NA values for latent observations
            lims = lims,
            #plot = plot,
            year = year,
            V2 = V2,
            #V3 = V3,
            plotZ = plotZ,
            yearZ = yearZ,
            x = x,
            yOrig = y,
            g2Levs = 3,
            g3Levs = 3,
            x3 = x3)
###################################### END OF DATA PREP

###########################################
## Initialisation of values for JAGS chains
###########################################
# Initial parameter values for JAGS
# To ensure conformity with the data all occupancies (ZI.*) are set to one
# Some other parameters are also fixed to avoid extemly small likelihoods
# but these can probably be relaxed a bit more.
zinit <- matrix(1, nrow = N, ncol = Y)
inits.fn <- function() list(z = zinit,
                            sd.m = runif(1,0,10),
                            sd.g2 = runif(1,0,10),
                            mean.g2 = runif(1,0,1),
                            sd.g3 = runif(1,0,10),
                            mean.g3 = runif(1,0,1),
                            mean.m = runif(1,0,1),
                            mean.p = runif(1,0,1),
                            phi = rep(runif(1,1,10), Y),
                            g1 = runif(1,-5,5),
                            cposLatent = c(0.025,0.15,0.375,0.625,0.85,0.975)[check$int]
                            )

######################################
## JAGS model (see end of script for model versions)
######################################
sink('scripts/JAGS/JAGS_simv8.txt')
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
    logit(p.Dec[a]) <- g0 + g1 * yOrig[a] + g3[x3[a]] #+ g2[x2[a]] 
    yOrig[a] ~ dunif(0,10) # this could be improved using the Wilson/Irvine approach (i.e. imputing new values where missing, 
                           # rather than relying on noise-inducing uniform draws)
    #x2[a] ~ dcat(pi) # (cover missing values for g2, or add in zeros?)
}

### Priors ###
## Intercept for state model for occupancy
mean.m ~ dbeta(1,1)
for(j in 1:Y){
  m[j] ~ dnorm(logit(mean.m), tau.m)
}  
tau.m <- 1/pow(sd.mA, 2)
sd.mA <- sd.m + .1 # see Kruschke page 486 (don't want shrinkage to be too strong when data sparse)
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
g1 ~ dnorm(0, 1) # weakly informative prior to help with quasi-complete separation in logistic regression

# Hierarchical prior for levels of grazing
#for (j in 1:g2Levs) { g2[j] ~ dnorm(0, g2tau) }
#g2tau <- 1/pow(sd.g2A, 2)
#sd.g2A <- sd.g2 + .1 # see Kruschke page 486 (don't want shrinkage to be too strong when data sparse)
#sd.g2 ~ dt(0, 0.01, 1)T(0,)

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

jagsModel <- jags.model(file= 'scripts/JAGS/JAGS_simv8.txt', data = Data, inits = inits.fn, n.chains = 3, n.adapt= 500)
para.names <- c('avgOcc', 'mu', 'm', 'C', 'cPos', 'z')
para.names <- c('mean.m', 'mPsi', 'mean.m0')
samples9 <- coda.samples(jagsModel, variable.names = para.names, n.iter = 500)
## Inspect results
out <- summary(samples9); out
mu_sd <- out$stat[,1:2] #make columns for mean and sd
q <- out$quantile[,c(3,1,5)] #make columns for median and CI
tableOut_m9 <- as.data.frame(cbind(mu_sd,q)) #make table
tableOut_m9
#save(tableOut_m9, samples9, file = "outputs/model9.Rdata")
#tableOut[grep(rownames(tableOut), pattern = 'gamma'),]
#tableOut[grep(rownames(tableOut), pattern = 'annOcc'),]
#write.csv(tableOut, file = "outputs/tests/test1_CS.csv")

plot(samples9)
gelman.diag(samples9)
## Using mcmcplots::mcmcplot
mcmcplots::mcmcplot(samples9, random = 10)
mcmcplots::caterplot(samples9, style = "plain")
samplesDF <- do.call(rbind.data.frame, samples9)
cor(samplesDF)
pairs(samplesDF)

#######################
# JAGS model versions
#######################
# JAGS_simv1.txt # with random walk prior/occupancy year effect (comments: variance for occupancy year effect in second year huge)
# JAGS_simv2.txt # with random walk prior, with occupancy year effect, but limiting variance
# JAGS_simv3.txt # with occupancy year effect, but not random walk prior
# JAGS_simv4.txt # with occupancy year effect, and detection model with intercept only
# JAGS_simv5.txt # withOUT occupancy year effect, and detection model with intercept only
# JAGS_simv6.txt # with occupancy year effect, but also hyperpriors providing shrinkage
# JAGS_simv8.txt # with occupancy year effect, but also hyperpriors providing shrinkage -- development/testing
# JAGS_simv9.txt # test of detection model only
