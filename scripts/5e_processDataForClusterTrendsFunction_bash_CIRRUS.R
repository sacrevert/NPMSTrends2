################################################################################
#### 5c. Create function to make model runs across multiple species easier.#####
          ################   CLUSTER VERSION ################   
################################################################################
# O.L. Pescott, olipes@ceh.ac.uk
# 10.03.2020
#rm(list=ls())
# v0 -- original version
# v1 -- amend function to remove all plots with no known states at time of processing
# v5c -- parrellelisation work -- this is probably inefficient for at least two reasons
# (1) first step is sequential (normal lapply statement)
# (2) exporting all data to all nodes, even though the nodes only need a subset of this
# test 2 -- use bash script to iterate across broad habitats
######################################
bHab <- commandArgs(trailingOnly = TRUE) # get broad habitat from command line
######################################
setwd("/home/pywell/olipes/npms2020_test2/")
library(snowfall) # Use sfLapply( x, fun, ... )
library(R2jags) ## requires that JAGS (preferably 4.3.0) is installed
# Initialise a the cluster across Cirrus nodes
# You do not need to edit these lines
hosts<-as.character(read.table(Sys.getenv('PBS_NODEFILE'),header=FALSE)[,1]) # read the nodes to use
sfSetMaxCPUs(length(hosts)) # ensure that snowfall can cope with this many hosts
sfInit(parallel=TRUE,type="MPI",cpus=length(hosts),useRscript=TRUE) # initialise the connection
sfLibrary(R2jags)

######################################
### Required files for development ###
sppPerHabL_P <- readRDS(file = "data/sppPerHabL_P_09032020.Rdata") # Positive indicator species data per broad habitat
sppPerHabL_N <- readRDS(file = "data/sppPerHabL_N_09032020.Rdata") # Negative indicator species data per broad habitat
focalSpp_P <- readRDS(file = "data/focalSpp_P.Rdata") # Positive indicator list per hab
focalSpp_N <- readRDS(file = "data/focalSpp_N.Rdata") # Negative indicator list per hab
######################################


###### Future improvments (noted 10.01.2019)
# Need to combine some indicator species during processing (script 3)
# e.g. Thymus polytrichus and Ulex gallii -- currently I am creating separate trends for these
######

# Prep
domins <- read.csv(file = "data/dominScores.csv", header = T, stringsAsFactors = F)
### NEED TO ALTER JAGS SCRIPT LOCATION ONCE ON CLUSTER ###
source(file = "scripts/X2_sinkJAGSscript.R") # write JAGS model to disk
source(file = "scripts/X3_sinkClusterFunction.R") # runModels_v5c_CLUS() function
file = 'scripts/JAGS_mod4.3.txt'
sfExportAll()

############################
## RUN MODELS IN PARALLEL ##
############################
### Apply JAGS model across broad habitats and **P** species (use sfLapply on cluster)
### Positive indicators
allPModels <- list()
ptm <- proc.time()
#allPModels <- lapply(seq_along(sppPerHabL_P[[1]]), # limit loop for testing
allPModels <- sfLapply(focalSpp_P[[bHab]],
                    function(i) runModels_v5c_CLUS(i = i, dat = sppPerHabL_P[[bHab]][i], file = file, bHab = bHab,
                      status = "P",
                      n.chains = 6, n.adapt = 100, n.iter = 5000, thin = 5)
                    )
Ptime <- (proc.time() - ptm); Ptime
# Add names to both list levels for convenience
names(allPModels) <- focalSpp_P[[bHab]]
save(allPModels, file = paste("outputs/allPModels_", bHab ,"_.Rdata", sep=""))
###

### Apply JAGS model across broad habitats and **N** species (use sfLapply on cluster)
### Negative indicators
allNModels <- list()
ptm <- proc.time()
allNModels <- sfLapply(focalSpp_N[[bHab]],
                       function(i) runModels_v5c_CLUS(i = i, dat = sppPerHabL_N[[bHab]][i], file = file, bHab = bHab,
                        status = "P",
                        n.chains = 6, n.adapt = 100, n.iter = 5000, thin = 5)
)
Ntime <- (proc.time() - ptm); Ntime
# Add names to both list levels for convenience
names(allNModels) <- focalSpp_N[[bHab]]
save(allNModels, file = paste("outputs/allNModels_", bHab, "_.Rdata", sep=""))
###
sfStop()
###