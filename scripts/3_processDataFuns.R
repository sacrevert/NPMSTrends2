## 3. Process retrieved data (plot/sample level info and species per sample info)
## in order to create per species information for modelling
# O.L. Pescott
# 17.08.2018
#rm(list=ls())

#source(file = "scripts/2_getDataFromIndiciaExec.R")
## Load data (previously extracted using functions in script 1)
## Or will be in the global environment if sourcing script "2_getData..."
load(file = "data/npms_SamplesSpecies_2020-01-03.Rdata")
load(file = "data/npms_PlotsSamples_2020-01-03.Rdata")


## Read in the official list of indicators with indicia preferred names and TVKs
# file has one row for every species x fine-scale habitat association
# broad habitat lists can be derived from another column indicating fine::broad habitat associations
inds <- read.csv(file = "data/npmsIndicatorsIndicia_Aug2018.csv", header = T, stringsAsFactors = F)
domins <- read.csv(file = "data/dominScores.csv", header = T, stringsAsFactors = F)
## NOTE THAT THIS FILE MIGHT NEED TO BE OCCASIONALLY UPDATED SO THAT IT IS IN LINE WITH THE PREFERRED TAXON ACCORDING TO THE TAXA TABLE IN INDICIA! ##
# UNFORTUNATELY THIS IS CURRENTLY A MANUAL TASK, BECAUSE THE HABITAT:SPECIES INFORMATION FOR WILDFLOWER/INDICATOR SPECIES IS NOT CONTAINED IN THE DATABASE ##
# working is here: W:\PYWELL_SHARED\Pywell Projects\BRC\_BRC_projects\NPMS\Technical work\Species
# last file: NPMS_FINAL_CONCISE_withBRC_IndiciaUpdateAug2018.xslx
#

## Need to produce lookup table to match between taxa actually recorded by surveyors at Inventory level, and the indicator list
#uniIndsList <- unique(inds$guidanceSpecies); write.csv(uniIndsList, file = "data/uniIndsList.csv")
#uniSpRecorded <- unique(npms15_18spp$preferred_taxon); write.csv(uniSpRecorded, file = "data/uniSpRecorded.csv")
###########################################################################################################
### Note that this will need updating every year as new recorded names are added to the indicia export ####
###########################################################################################################
unifyNames <- read.csv(file = "data/unifyNames.csv", header = T, stringsAsFactors = F)
npms_spp <- merge(npms_spp, unifyNames, by.x = "preferred_taxon", by.y = "nameRec", all.x = T, all.y = F)

## Process fie/broad-scale indicator data for later function
head(inds)
indsF <- aggregate(guidanceSpecies ~ indiciaName + indiciaPrefTvk + fine.scale_habitat + indicator + indicator_type + wildflower,
                   data = inds, function(x) length(x))
indsB <- aggregate(guidanceSpecies ~ indiciaName + indiciaPrefTvk + broad.scale_habitat + indicator + indicator_type + wildflower,
                   data = inds, function(x) length(x))
indsF <- indsF[,-(7)]; indsB <- indsB[,-(7)]; # delete the unnecessary aggregated columns

## These steps are all about making a lookup table that can later tell us, for every possible indicator/wildflower v. broad/fine habitat combination,
# whether species X is on the list for that combination
indsF$combinedI <- ifelse(indsF$indicator == "y", paste(indsF$fine.scale_habitat,", Indicator survey", sep = ""), NA)
indsF$combinedW <- ifelse(indsF$wildflower == "y", paste(indsF$fine.scale_habitat,", Wildflower survey", sep = ""), NA)
indsB$combinedI <- ifelse(indsB$indicator == "y", paste(indsB$broad.scale_habitat,", Indicator survey", sep = ""), NA)
indsB$combinedW <- ifelse(indsB$wildflower == "y", paste(indsB$broad.scale_habitat,", Wildflower survey", sep = ""), NA)

# for fine-scale habs
indsF_tmp1 <- indsF[,c("indiciaName", "indiciaPrefTvk", "combinedI")]
indsF_tmp2 <- na.omit(indsF[,c("indiciaName", "indiciaPrefTvk", "combinedW")])
names(indsF_tmp1)[3] <- "combined"; names(indsF_tmp2)[3] <- "combined"
indsF_minimal <- rbind(indsF_tmp1, indsF_tmp2)

# for broad-scale habs
indsB_tmp1 <- indsB[,c("indiciaName", "indiciaPrefTvk", "combinedI")]
indsB_tmp2 <- na.omit(indsB[,c("indiciaName", "indiciaPrefTvk", "combinedW")])
names(indsB_tmp1)[3] <- "combined"; names(indsB_tmp2)[3] <- "combined"
indsB_minimal <- rbind(indsB_tmp1, indsB_tmp2)

indsLookup <- rbind(indsF_minimal, indsB_minimal) # combine

## function to select the relevant samples for any given set of habitats
getSamples <- function(habsList){ temp <- merge(npms_spp, npms_plots, by.x = 'sample_id', by.y = "sample", all.x = T, all.y = F)
                                  temp2 <- temp[temp$surv_habitat %in% habsList,]
                                  return(temp2)
}
## 2019 addition:
# Add unique() just to remove any database duplicates (one only for graze dat)
grazeDat <- unique(npms_plots[npms_plots$caption == "NPMS Grazing",])
names(grazeDat)[6] <- "grazing"
mowDat <- unique(npms_plots[npms_plots$surv_habitat == "Cutting / mowing",])
names(mowDat)[6] <- "cutting"

# END function
##### 03 09 2019: at some point I will need to harmonise the names that are taken directly from the database, e.g., at the moment infraspecific taxa etc.
##### within indicator taxa will be ignored. (03 01 2020: hang on, I did the unify names bit above, so this should be OK if that lookup list is up-to-date,
##### I think me from 03 09 2019 overlooked that step)

## Grasslands examples
grasslands <- c("Neutral pastures and meadows", "Dry acid grassland", "Dry calcareous grassland", "Neutral damp grassland", "Lowland grassland")
grassSamples <- getSamples(habsList = grasslands)
#

## Then process the filtered data so that species abundance information is as it should be (i.e. present/absent/NA)
# try Achillea millefolium first (outside of function, in order to develop principle)
#spSamplePA <- function(samples, species){
#                                         temp <- aggregate(preferred_taxon ~ sample_id, data = grassSamples, function(x) max(ifelse(x == "Achillea millefolium", 1, 0)))
#                                         presence <- temp[temp$preferred_taxon == 1,] # samples containing the taxon of interest
#                                         presSamps <- merge(presence, npms15_18plots[npms15_18plots$caption == "NPMS Habitat",], by.x = "sample_id", by.y = "sample", all.x = T, all.y = F) # useful for subsequent rbind
#                                         presSamps$combination <- paste(presSamps$surv_habitat,", ",presSamps$title, sep = "") # useful for subsequent rbind
#                                         absence <- temp[temp$preferred_taxon == 0,] # samples NOT containing the taxon of interest (but are these 0 or NA?)
#                                         absSamps <- merge(absence, npms15_18plots[npms15_18plots$caption == "NPMS Habitat",], by.x = "sample_id", by.y = "sample", all.x = T, all.y = F) # habitat/level info for every sample
#                                         absSamps$combination <- paste(absSamps$surv_habitat,", ",absSamps$title, sep = "") # create column for lookup to indsLookup
#                                         indsLookup_fil <- indsLookup[indsLookup$indiciaName == "Achillea millefolium",]
#                                         absSamps_AN <- merge(absSamps, indsLookup_fil, by.x = "combination", by.y = "combined", all.x = T, all.y = F)# absent or NA indicator
#                                         presSamps_AN <- merge(presSamps, indsLookup_fil, by.x = "combination", by.y = "combined", all.x = T, all.y = F)# again, for rbind
#                                         absSamps_AN$preferred_taxon <- ifelse(is.na(absSamps_AN$indiciaName), "NA", 0)
#                                         # but inventory samples should always be zero if absent (in theory at least)
#                                         absSamps_AN$preferred_taxon <- ifelse(absSamps_AN$title =="Inventory survey", 0, absSamps_AN$preferred_taxon)
#                                         samples_PAN <- rbind(presSamps_AN, absSamps_AN) # combine so we have a 'PAN' (pres, abs, NA) view on our focal species
                                        # return(samples_PAN)
                                        # }

###########################################
## Generalise to make the function function
###########################################
#spSamplePA <- function(samples, species){ tryCatch(
#                                          {temp <- aggregate(preferred_taxon ~ sample_id + date, data = samples, function(x) max(ifelse(x == species, 1, 0)))
#                                          presence <- temp[temp$preferred_taxon == 1,] # samples containing the taxon of interest
#                                          presSamps <- merge(presence, npms15_18plots, by.x = "sample_id", by.y = "sample", all.x = T, all.y = F) # useful for subsequent rbind
#                                          presSamps$combination <- paste(presSamps$surv_habitat,", ",presSamps$title, sep = "") # useful for subsequent rbind
#                                          absence <- temp[temp$preferred_taxon == 0,] # samples NOT containing the taxon of interest (but are these 0 or NA?)
#                                          absSamps <- merge(absence, npms15_18plots, by.x = "sample_id", by.y = "sample", all.x = T, all.y = F) # habitat/level info for every sample
#                                          absSamps$combination <- paste(absSamps$surv_habitat,", ",absSamps$title, sep = "") # create column for lookup to indsLookup
#                                          indsLookup_fil <- indsLookup[indsLookup$indiciaName == species,]
#                                          absSamps_AN <- merge(absSamps, indsLookup_fil, by.x = "combination", by.y = "combined", all.x = T, all.y = F)# absent or NA indicator
#                                          presSamps_AN <- merge(presSamps, indsLookup_fil, by.x = "combination", by.y = "combined", all.x = T, all.y = F)# again, for rbind
#                                          absSamps_AN$preferred_taxon <- ifelse(is.na(absSamps_AN$indiciaName), "NA", 0)
#                                          # but inventory samples should always be zero if absent
#                                          absSamps_AN$preferred_taxon <- ifelse(absSamps_AN$title =="Inventory survey", 0, absSamps_AN$preferred_taxon)
#                                          samples_PAN <- rbind(presSamps_AN, absSamps_AN) # combine so we have a PAN (pres, abs, NA) view on our focal species
#                                          names(samples_PAN)[4] <- "PAN" # rename to avoid confusion
#                                          samples_PAN$PAN <- as.numeric(samples_PAN$PAN) 
#                                          # add in domins and dates etc. from original sample/species data
#                                          samples_PAN <- merge(samples_PAN, npms15_18spp[npms15_18spp$preferred_taxon == species,], by.x = "sample_id", by.y = "sample_id", all.x = T, all.y = F)
#                                          samples_PAN$domin <- ifelse(samples_PAN$PAN == 0, 0, samples_PAN$domin)
#                                          samples_PAN <- merge(samples_PAN, domins, by.x = "domin", by.y = "dominOrig", all.x = T, all.y = F)
#                                          return(samples_PAN)}, 
#                                          error = function(err) print(species) )
#}

# update so using "UnifiedPlusInd" column in npms15_18spp, and incorporates info on grazing or cutting
spSamplePA_v1.1 <- function(samples, species){ tryCatch(
    {temp <- aggregate(UnifiedPlusInd ~ sample_id + date, data = samples, function(x) max(ifelse(x == species, 1, 0)))
    presence <- temp[temp$UnifiedPlusInd == 1,] # samples containing the taxon of interest
    presSamps <- merge(presence, npms_plots[npms_plots$caption == "NPMS Habitat",], by.x = "sample_id", by.y = "sample", all.x = T, all.y = F) # useful for subsequent rbind
    presSamps$combination <- paste(presSamps$surv_habitat,", ",presSamps$title, sep = "") # useful for subsequent rbind
    absence <- temp[temp$UnifiedPlusInd == 0,] # samples NOT containing the taxon of interest (but are these 0 or NA?)
    absSamps <- merge(absence, npms_plots[npms_plots$caption == "NPMS Habitat",], by.x = "sample_id", by.y = "sample", all.x = T, all.y = F) # habitat/level info for every sample
    absSamps$combination <- paste(absSamps$surv_habitat,", ",absSamps$title, sep = "") # create column for lookup to indsLookup
    indsLookup_fil <- indsLookup[indsLookup$indiciaName == species,]
    absSamps_AN <- merge(absSamps, indsLookup_fil, by.x = "combination", by.y = "combined", all.x = T, all.y = F)# absent or NA indicator
    presSamps_AN <- merge(presSamps, indsLookup_fil, by.x = "combination", by.y = "combined", all.x = T, all.y = F)# again, for rbind
    absSamps_AN$UnifiedPlusInd <- ifelse(is.na(absSamps_AN$indiciaName), "NA", 0)
    # but inventory samples should always be zero if absent
    absSamps_AN$UnifiedPlusInd <- ifelse(absSamps_AN$title =="Inventory survey", 0, absSamps_AN$UnifiedPlusInd)
    samples_PAN <- rbind(presSamps_AN, absSamps_AN) # combine so we have a PAN (pres, abs, NA) view on our focal species
    names(samples_PAN)[4] <- "PAN" # rename to avoid confusion
    samples_PAN$PAN <- as.numeric(samples_PAN$PAN) 
    # add in domins and dates etc. from original sample/species data
    ## change here from v1.0 (use 'UnifiedPlusInd' rather than original "preferred_taxon")
    samples_PAN <- merge(samples_PAN, npms_spp[npms_spp$UnifiedPlusInd == species,], by.x = "sample_id", by.y = "sample_id", all.x = T, all.y = F)
    samples_PAN$domin <- ifelse(samples_PAN$PAN == 0, 0, samples_PAN$domin)
    samples_PAN <- merge(samples_PAN, domins, by.x = "domin", by.y = "dominOrig", all.x = T, all.y = F)
    samples_PAN <- merge(samples_PAN, grazeDat[,c("sample","grazing")], by.x = "sample_id", by.y = "sample", all.x = T, all.y = F) # add grazing data
    samples_PAN$grazing <- droplevels(samples_PAN$grazing)
    samples_PAN <- merge(samples_PAN, mowDat[,c("sample","cutting")], by.x = "sample_id", by.y = "sample", all.x = T, all.y = F) # add mowing/cutting data
    samples_PAN$cutting <- droplevels(samples_PAN$cutting)
    return(samples_PAN)}, 
    error = function(err) print(species) )
}

## Let's see if it works! (looks ok - 05 09 2019)
grasslands <- c("Neutral pastures and meadows", "Dry acid grassland", "Dry calcareous grassland", "Neutral damp grassland", "Lowland grassland")
grassSamples <- getSamples(habsList = grasslands)
Achi_mill_PAN <- spSamplePA_v1.1(samples = grassSamples, species = "Achillea millefolium") # Seems good (now 3979 rows)
save(Achi_mill_PAN, file = paste("data/Achi_mille_grassSamples_", as.character(Sys.Date()), ".Rdata", sep = ""))
#head(Achi_mill_PAN)
#load(file = "data/Achi_mille_grassSamples_20180920.Rdata") # old data (trends proj 1) for comparison
