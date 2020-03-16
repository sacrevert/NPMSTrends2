################################################################################
#### 5d. Create function to make model runs across multiple species easier.#####
          ################   CLUSTER VERSION ################   
################################################################################
# O.L. Pescott, olipes@ceh.ac.uk
# 10.03.2020
#rm(list=ls())
# v0 -- original version
# v1 -- amend function to remove all plots with no known states at time of processing
# v5c -- parrellelisation work -- this approach isn't going to work on cluster (at least not in a simple way)
# v5d -- use PBS or bash script to manage first bit of parrellelisation
######################################
#library(snowfall) # Use sfLapply( x, fun, ... )
library(R2jags) ## obviously requires that JAGS (preferably 4.3.0) is installed
#sfInit(parallel=TRUE, cpus=2)

######################################
### Required files for development ### # Change to readRDS assignments to objects for cluster
load(file = "outputs/sppPerHabL_P_09032020.Rdata") # Positive indicator species data per broad habitat
load(file = "outputs/sppPerHabL_N_09032020.Rdata") # Negative indicator species data per broad habitat
load(file = "data/focalSpp_P.Rdata") # Positive indicator list per hab
load(file = "data/focalSpp_N.Rdata") # Negative indicator list per hab
######################################

###### Future improvments (noted 10.01.2019)
# Need to combine some indicator species during processing (script 3)
# e.g. Thymus polytrichus and Ulex gallii -- currently I am creating separate trends for these
######

# Prep
domins <- read.csv(file = "data/dominScores.csv", header = T, stringsAsFactors = F)
### NEED TO ALTER JAGS SCRIPT LOCATION ONCE ON CLUSTER ###
source(file = "scripts/X2_sinkJAGSscript.R")
source(file = "scripts/X3_sinkClusterFunction.R")
file = 'scripts/JAGS_mod4.3.txt'

############################
## RUN MODELS IN PARALLEL ##
############################
### Apply JAGS model across broad habitats and **P** species (use sfLapply on cluster)
### Positive indicators
allPModels <- list()
ptm <- proc.time()
allPModels <- lapply(seq_along(sppPerHabL_P), # for each broad hab
                    function(x) lapply(focalSpp_P[[x]], # for each species within a broad hab
                                       # Note that "dat" object is also used within function to get species name
                                       # i = species name, dat = relevant list of habitat samples for sp., file = JAGS model
                                       function(i) runModels_v5c_CLUS(i = i, dat = sppPerHabL_P[[x]][i], file = file)
                                       )
                    )
proc.time() - ptm
# Add names to both list levels for convenience
names(allPModels) <- names(sppPerHabL_P)
for (i in 1:11){
  names(allPModels[[i]]) <- focalSpp_P[[i]]
}
save(allPModels, file = "outputs/allPModels.Rdata")
###

### Apply JAGS model across broad habitats and **N** species (use sfLapply on cluster)
### Negative indicators
allNModels <- list()
allNModels <- lapply(seq_along(sppPerHabL_N), # normally sequential lapply for each broad hab
                     function(x) lapply(focalSpp_N[[x]], # snowfall sfLapply for each species within a broad hab
                                        # Note that "dat" object is also used within function to get species name
                                        # i = species name, dat = relevant list of habitat samples for sp., file = JAGS model
                                        function(i) runModels_v5c_CLUS(i = i, dat = sppPerHabL_N[[x]][i], file = file)
                                        )
                    )
# Add names to both list levels for convenience
names(allNModels) <- names(sppPerHabL_N)
for (i in 1:11){
  names(allNModels[[i]]) <- focalSpp_N[[i]]
}
save(allNModels, file = "outputs/allNModels.Rdata")
###
