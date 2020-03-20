# Produce species plots from the bash cluster approach (script v. 5e)
# 16.03.2020
# O.L. Pescott, olipes@ceh.ac.uk
#rm(list=ls())
library(ggplot2)
library(coda)
## recall that files saved as e.g.: save(allNModels, file = paste("outputs/allNModels_", bHab, "_.Rdata", sep=""))
## Remember list of broad hab names
#bHab <- c("lowGrass", "arable", "bogWHeath", "woodsEtc", "coast", "freshwater", "heaths", "marshFen",
#          "highGrass", "pineWoods", "rocks")
# Load per hab species lists
focalSpp_P <- readRDS(file = "data/focalSpp_P.Rdata") # Positive indicator list per hab
focalSpp_N <- readRDS(file = "data/focalSpp_N.Rdata") # Negative indicator list per hab

## Single habitat test
#load(file = "outputs/test2_bash/allPModels_lowGrass_.Rdata")
#
#allPOuts <- allPModels
##
allPOuts <- list()
allNOuts <- list()
pFiles <- list.files("outputs/test2_bash/P")
nFiles <- list.files("outputs/test2_bash/N")
for (i in 1:11){
  allPOuts[i] <- mget(load(file = paste("outputs/test2_bash/P/", pFiles[i], sep="")))
  allNOuts[i] <- mget(load(file = paste("outputs/test2_bash/N/", nFiles[i], sep="")))
}

## Within allPOuts list, we have the following levels (L)
# L1 = broad H
# L2 = species
# L3 = things returned from JAGS

### Function for plotting for an individual species (i.e. working at L3 of the above) ###
# type = P or N species
mkPlot <- function(spp = i, bHab, type, dat, par = "mC", yrRange = c(2015:2019)) {
  dat2 <- dat[[3]] # get samples
  dat2S <- summary(dat2)
  allMuSe <- data.frame(dat2S$stat[,c(1,4)]) # mu and se from summary
  mC2Plot <- allMuSe[grep(par, rownames(allMuSe), value = T), ] # select correct row output
  mC2Plot$Year <- yrRange
  mC2Plot$low <- mC2Plot$Mean - mC2Plot$Time.series.SE
  mC2Plot$high <- mC2Plot$Mean + mC2Plot$Time.series.SE
  mC2Plot$species <- spp
  ggplot(mC2Plot, aes_string(x = "Year", y = "Mean")) +
    theme_bw() +
    geom_ribbon(data = mC2Plot,
                aes_string(group = 1, ymin = "low", ymax = "high"),
                alpha = 0.2) +
    geom_line(size = 1, col = "black") +
    ylab("Mean ZI cover") +
    xlab("Time period") +
    scale_y_continuous(limits = if( min(mC2Plot$low - 0.1) > 0) { 
      c( min(mC2Plot$low - 0.1), max(mC2Plot$high + 0.1)) } else {
        c(0,0.15)} ) +
    ggtitle(paste(spp, "--", bHab, "-- ", type)) + 
    theme(plot.title = element_text(lineheight = .8, face = "bold"),
          legend.position = 'bottom')
}

##########################################################################################

##############################
bHabNames <- c("arable", "bogWHeath","coast", "freshwater", "heaths", "highGrass", "lowGrass", "marshFen",
               "pineWoods", "rocks", "woodsEtc")
## Reorder lists of species alphabetically by bHab
alphFocalSpp_P <- focalSpp_P[sort(names(focalSpp_P))]
alphFocalSpp_N <- focalSpp_N[sort(names(focalSpp_N))]
###
### Create plots for P species ###
pPlots <- lapply(seq_along(allPOuts), # for each numbered broad hab
       function(x) { lapply(alphFocalSpp_P[[x]], 
       #function(x) { lapply(alphFocalSpp_P[[bHab]][1:2], # Test
                            function(i) tryCatch( { mkPlot(spp = i, bHab = bHabNames[x], type = "P", par = "mC",
                                                      dat = allPOuts[[x]][[i]], yrRange = c(2015:2019)) },
                                                  ## insert NULL if, for some reason, the JAGS model failed
                                                  error = function(err) NULL
                                                  )
                            )
         }
       )
### NAMING ###
## These are alphabetical as they were read back in from file
names(pPlots) <- c("arable", "bogWHeath","coast", "freshwater", "heaths", "highGrass", "lowGrass", "marshFen",
                      "pineWoods", "rocks", "woodsEtc")
for (i in 1:11){
  names(pPlots[[i]]) <- alphFocalSpp_P[[i]]
}
#save(pPlots, file = "outputs/plots/test2_bash/pPlots.Rdata")
#load(file = "outputs/plots/test2_bash/pPlots.Rdata")
##############################
### Create plots for N species ###
nPlots <- lapply(seq_along(allNOuts), # for each numbered broad hab
       function(x) { lapply(alphFocalSpp_N[[x]], 
                            function(i) tryCatch( { mkPlot(spp = i, bHab = bHabNames[x], type = "N", par = "mC",
                                                    dat = allNOuts[[x]][[i]], yrRange = c(2015:2019)) },
                                                  ## insert NULL if, for some reason, the JAGS model failed
                                                  error = function(err) NULL
                                                  )
                            )
       }
     )
### NAMING ###
## These are alphabetical as they were read back in from file
names(nPlots) <- c("arable", "bogWHeath","coast", "freshwater", "heaths", "highGrass", "lowGrass", "marshFen",
                     "pineWoods", "rocks", "woodsEtc")
for (i in 1:11){
  names(nPlots[[i]]) <- alphFocalSpp_N[[i]]
}
#save(nPlots, file = "outputs/plots/test2_bash/nPlots.Rdata")
#load(file = "outputs/plots/test2_bash/nPlots.Rdata")

#################### SAVING PLOTS TO FILES BELOW ####################
##########################
bHabNames <- c("arable", "bogWHeath","coast", "freshwater", "heaths", "highGrass", "lowGrass", "marshFen",
               "pineWoods", "rocks", "woodsEtc")
##########################
### Save plot function ###
savePlots <- function(x = x, spp = i, bHab = bHab, type, subfolder = "test2_bash") { 
  file_name <- paste("outputs/plots/", subfolder, "/", spp , "_", bHab, "_", type, ".png", sep="")
  png(file_name)
  if (type == "P") {
  print(pPlots[[x]][[spp]])
  } else {
  print(nPlots[[x]][[spp]])
  }
  dev.off()
}
################################################
#################### SLOW!! ####################
### Save P plots to png - separate file for each plot ###
lapply(seq_along(allPOuts), # for each numbered broad hab
       function(x) { lapply(alphFocalSpp_P[[x]], 
                            function(i) savePlots(x = x, spp = i, bHab = bHabNames[x], type = "P")
        )
      }
)

### Save N plots to png - separate file for each plot ###
lapply(seq_along(allNOuts), # for each numbered broad hab
       function(x) { lapply(alphFocalSpp_N[[x]], 
                            function(i) savePlot(spp = i, bHab = bHabNames[x], type = "P")
        )
      }
)
################################################
