######################################
## JAGS model
######################################
sink('scripts/JAGS_mod4.3.txt')
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