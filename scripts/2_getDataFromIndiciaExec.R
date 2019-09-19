## 2. Execute function to retrieve plot/sample/habitat/level data from Indicia/NPMS database (to avoid hardcoding passwords in script)
# O.L. Pescott
# 17.08.2018
# 09.01.2019, changes include automatic dating of saved sample files, and database password request (in RStudio only)

# SQL query is also in W:\PYWELL_SHARED\Pywell Projects\BRC\_BRC_projects\NPMS\Analyses\2018 08 - Per species trend analyses\SQL\extractData_v0.0.sql

####
#### Warning! This is slow (11 mins if retrieving data up between 2015-mid-2018) and will get increasingly slower
#### Need to investigate possible speed-ups/new cache table that runs this automatically on some periodic basis
####
#rm(list=ls())

## Check for rstudioapi and install if not available
list.of.packages <- c("rstudioapi")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(rstudioapi)
##
#source(file = "scripts/1_getDataFromIndiciaFuns.R") # source SQL functions
##
# Get db password from user (RStudio only)
password <- rstudioapi::askForPassword("Database password")
##

## Get plot and sample data!
# columns: plot_id, monad, sample, title (survey), surv_habitat
npms15_18plots <- getNpmsData_PlotsSamples_v1.1(password = password)
save(npms15_18plots, file = paste("data/npms1518_PlotsSamples_", as.character(Sys.Date()), ".Rdata", sep = ""))

## Get taxon data across samples
npms15_18spp <- getNpmsData_SamplesSpecies(password = password)
# change this so that file name automatically updates with system date - 09.01.2019
#save(npms15_18spp, file = "data/npms1518_SamplesSpecies_17Aug2018.Rdata") # 101,508 -- old version - 17.08.2018
save(npms15_18spp, file = paste("data/npms1518_SamplesSpecies_", as.character(Sys.Date()), ".Rdata", sep =""))
#load(file = "data/npms1518_SamplesSpecies_2019-09-04.Rdata")
rm(password)
# regarding the row count discrepancies, see the two versions of the SQL in 
# W:\PYWELL_SHARED\Pywell Projects\BRC\_BRC_projects\NPMS\Analyses\2018 08 - Per species trend analyses\SQL\extractData_v0.0.sql

## END