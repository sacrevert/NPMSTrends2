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

## Test!
load(file = "outputs/allPModels_arable_.Rdata")
allPOuts <- list()
allPOuts[[1]] <- allPModels
##

#allPOuts <- lapply(Sys.glob("outputs/allPModels_*.Rdata"), load)
#allNOuts <- lapply(Sys.glob("outputs/allNModels_*.Rdata"), load)

## Within allPOuts list, we have the following levels (L)
# L1 = broad H
# L2 = species
# L3 = things returned from JAGS

# Function for plotting for an individual species (i.e. working at L3 of the above)
# type = P or N species
mkPlot <- function(spp = i, bHab, type, dat, yrRange, par = "mC", parlims = c(0,1)) {
  dat2 <- dat[[3]] # get samples
  dat2S <- summary(dat2S)
  allMuSe <- data.frame(dat2S$stat[,c(1,3)]) # mu and se from summary
  mC2Plot <- allMuSe[rownames(allMuSe) %like% par, ] # select correct row output
  mC2Plot$Year <- yrRange
  mC2Plot$low <- mC2Plot$Mean - mC2Plot$Naive.SE
  mC2Plot$high <- mC2Plot$Mean + mC2Plot$Naive.SE
  mC2Plot$species <- spp
  pdf(file = paste("outputs/plots/", bHab, "_", i, "_", type, ".pdf"), width = 6, height = 4, onefile = T)
  thePlot <- ggplot(mC2Plot, aes_string(x = "Year", y = "Mean")) +
    theme_bw() +
    geom_ribbon(data = mC2Plot,
                aes_string(group = 1, ymin = "low", ymax = "high"),
                alpha = 0.2) +
    geom_line(size = 1, col = "black") +
    ylab("Mean ZI cover") +
    xlab("Time period") +
    scale_y_continuous(limits = parlims) +
    ggtitle(mC2Plot$species) + 
    theme(plot.title = element_text(lineheight = .8, face = "bold"),
          legend.position = 'bottom')
  dev.off()
}

## Output graphs for P species
lapply(seq_along(allPOuts), # for each numbered broad hab
       function(x) { lapply(focalSpp_P[[x]], 
                            function(i) mkPlot(spp = i, bHab = names(focalSpp_P[x]), type = "P",
                                               dat = allPOuts[[x]][i]), yrRange = c(2015:2019))
         })

## Output graphs for N species
lapply(seq_along(allNOuts), # for each numbered broad hab
       function(x) { lapply(focalSpp_N[[x]], 
                            function(i) mkPlot(spp = i, bHab = names(focalSpp_N[x]), type = "N",
                                               dat = allNOuts[[x]][i]), yrRange = c(2015:2019))
       })





