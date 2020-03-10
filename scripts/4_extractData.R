## 4. Run extraction and processing for broad habitat X and species Y
## prepare Domin data for modelling
# O.L. Pescott
# 21.08.2018, updated 10.01.2019
#rm(list=ls())

## Source file with datasets, getSamples() and spSamplePA() functions
## Might be slow, as script 3_process... itself sources script 2_getDataFrom..., which queries the Indicia database
##############################################################
################ Requires database password!!! ###############
############### or you can load existing data ################
############### from: file = "data/*"  #######################
######## See script 3_processDataFuns.R lines 10 and 11 ######
######## choose latest data files and run from there #########
##############################################################
source(file = "scripts/3_processDataFuns.R")

##############################################################
#### Single example of Achillea millefolium in grasslands ####
############################################################## #### collapse ####
#grasslands <- c("Neutral pastures and meadows", "Dry acid grassland", "Dry calcareous grassland", "Neutral damp grassland", "Lowland grassland")
#grassSamples <- getSamples(habsList = grasslands)
#Achi_mill_PAN <- spSamplePA(samples = grassSamples, species = "Achillea millefolium")
#head(Achi_mill_PAN); tail(Achi_mill_PAN)
#load(file = "data/Achi_mille_grassSamples_2020-01-03.Rdata")
##

## Following can be used if database extraction not run
#load(file = "data/npms1518_SamplesSpecies_2019-09-05.Rdata.Rdata")
#load(file = "data/npms1518_PlotsSamples_2019-09-05.Rdata.Rdata")


##############################################################
############ Generalised approach using functions ############
##############################################################
#### hab samples now in allHabSamples list, see script 3 ####
## General habitats list
#habs <- read.csv( file = "data/NPMShabitats_Surv_to_Broad.csv", h = T, strings = F)
# bsh = broad scale habitat
#getSamps2 <- function(bsh) { getSamples(habsList = unique(habs[habs$NPMS.broad.habitat==bsh,2]))}
bsh <- c("Lowland grassland", "Arable margins", "Bog and wet heath", "Broadleaved woodland, hedges and scrub", # same order as used for sample extraction
         "Coast", "Freshwater", "Heathland",  "Marsh and fen", "Upland grassland",
         "Native pinewood and juniper scrub", "Rock outcrops, cliffs and scree")
#habSamps <- getSamps2(bsh = "Lowland grassland") # includes relevant broad and fine scale habitat names

## List of positive indicator species per fine habitat set (species subset = spp)
inds <- read.csv(file = "data/npmsIndicatorsIndicia_Aug2018.csv", header = T, stringsAsFactors = F) # repeated from script 3
#
selectSp <- function(bsh, type = "positive") {tmp <- unique(inds[inds$broad.scale_habitat == bsh & inds$indicator_type == type,]$indiciaName)} # type can also be "negative"
# get species list for relevant broad habitat
focalSpp_P <- list()
for (i in 1:11){ 
  focalSpp_P[[i]] <- selectSp(bsh = bsh[i])
}
names(focalSpp_P) <- names(allHabSamples)

# for negative indicators
focalSpp_N <- list()
for (i in 1:11){ 
  focalSpp_N[[i]] <- selectSp(bsh = bsh[i], type = "negative")
}
names(focalSpp_N) <- names(allHabSamples)

sppPerHabL_P <- list() # lists of samples for POSITIVE species per hab
sppPerHabL_N <- list() # lists of samples for NEGATIVE species per hab

# error causing species will currently just be a character element to the list, rather than a nested df
# Function "spSamplePA_v1.1()" will print names of species with no data within a habitat to screen
sppPerHabL_P <- lapply(seq_along(allHabSamples), 
                        function(x) lapply(focalSpp_P[[x]], # Positive indicator species
                                            function(z) spSamplePA_v1.1(samples = allHabSamples[[x]], species = z)
                                           )
                      )
names(sppPerHabL_P) <- names(allHabSamples)
for (i in 1:11){
  names(sppPerHabL_P[[i]]) <- focalSpp_P[[i]]
}

sppPerHabL_N <- lapply(seq_along(allHabSamples), 
                       function(x) lapply(focalSpp_N[[x]], # Negative indicator species
                                          function(z) spSamplePA_v1.1(samples = allHabSamples[[x]], species = z)
                                          )
                      )
names(sppPerHabL_N) <- names(allHabSamples)
for (i in 1:11){
  names(sppPerHabL_N[[i]]) <- focalSpp_N[[i]]
}

#save(sppPerHabL_P, file = "outputs/sppPerHabL_P_09032020.Rdata")
#save(sppPerHabL_N, file = "outputs/sppPerHabL_N_09032020.Rdata")
#library(rlist)
#summary(sppPerHabL_P)
load(file = "outputs/sppPerHabL_P_09032020.Rdata")
load(file = "outputs/sppPerHabL_N_09032020.Rdata")

## Species with no data for reference
sppHabsNoData <- data.frame(noData =   
  unlist(
    lapply(seq_along(sppPerHabL_P),
         function(x) lapply(seq_along(focalSpp_P[[x]][]), 
          function(i) ifelse(which(!is.data.frame(sppPerHabL_P[[x]][[i]]))==1, paste(names(focalSpp_P[x]), sep ="_", focalSpp_P[[x]][[i]]), NULL) 
         )
    )
  )
)

### PREVIOUS CODE FOR ONE HABITAT ###
#sppDatList <- lapply(focalSpp, function(x) spSamplePA_v1.1(samples = habSamps, species = x)) # no data species are printed to screen
#names(sppDatList) <- focalSpp
#excludedSppHab_P <- unlist(lapply(seq_along(sppDatList), function(i) ifelse(which(!is.data.frame(sppDatList[[i]]))==1, names(sppDatList)[i], NULL)))
#excludedSppHab_N <- unlist(lapply(seq_along(sppDatList), function(i) ifelse(which(!is.data.frame(sppDatList[[i]]))==1, names(sppDatList)[i], NULL)))
##
########
# this now has cutting/mowing and grazing data, and is unified to indicator names
#save(sppDatList, excludedSpp, file = paste("outputs/grasslandsEg_", as.character(Sys.Date()), ".Rdata", sep = ""))
#load(file = "outputs/grasslandsEg_09 09 2019.Rdata")
#load(file = "outputs/grasslandsEg_26 03 2019.Rdata") # previous run (first indicator project) for comparison
#####
# The next step is to run the model over each species data frame
# (i.e. run script 5_processData..., but over multiple species, and collecting useful results)
#####
