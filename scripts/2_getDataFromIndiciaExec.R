## 2. Execute function to retrieve plot/sample/habitat/level data from Indicia/NPMS database (to avoid hardcoding passwords in script)
# O.L. Pescott
# 17.08.2018
# 09.01.2019, changes include automatic dating of saved sample files, and database password request (in RStudio only)
# SQL query is also in W:\PYWELL_SHARED\Pywell Projects\BRC\_BRC_projects\NPMS\Analyses\2018 08 - Per species trend analyses\SQL\extractData_v0.0.sql

#rm(list=ls())
##
source(file = "scripts/1_getDataFromIndiciaFuns.R") # source SQL functions
##
##

## Get plot and sample data!
# columns: plot_id, monad, sample, title (survey), surv_habitat
## Slighly nervous about the lack of any increase in plot numbers from 2018 dataset (30,519), but I suppose it has to plateau somewhere, and the dates in the functions have definitely changed!
# e.g.
getNpmsData_PlotsSamples_v1.1 # look at stored function
#npms_plots <- getNpmsData_PlotsSamples_v1.1(password = password) #6922 for 2019 in theory
npms_plots <- getNpmsData_PlotsSamples_v1.1(Connection = Connection)
# manually load as there is an issue with the SQL script not returning the latest (2019) data through RODBC
#npms_plots <- read.csv(file = "data/npms_PlotsSamples_2020-01-03.csv", header = T, stringsAsFactors = T) #strings as factors T to match db extraction
#save(npms_plots, file = "data/npms_PlotsSamples_2020-01-03.Rdata")
#load(file = "data/npms_PlotsSamples_2020-03-04.Rdata") # inc. 2019 data
save(npms_plots, file = paste("data/npms_PlotsSamples_", as.character(Sys.Date()), ".Rdata", sep = ""))

## Get taxon data across samples
npms_spp <- getNpmsData_SamplesSpecies(Connection = Connection)
#npms_spp <- read.csv(file = "data/npms_SamplesSpecies_2020-01-03.csv", header = T, stringsAsFactors = T)
#save(npms_spp, file = "data/npms_SamplesSpecies_2020-01-03.Rdata")
#load(file = "data/npms_SamplesSpecies_2020-03-04.Rdata")  # inc. 2019 data
# change this so that file name automatically updates with system date - 09.01.2019
#save(npms_spp, file = "data/npms_SamplesSpecies_17Aug2018.Rdata") # 101,508 -- old version - 17.08.2018
save(npms_spp, file = paste("data/npms_SamplesSpecies_", as.character(Sys.Date()), ".Rdata", sep =""))
#load(file = "data/npms_SamplesSpecies_2020-03-04.Rdata")
rm(pwd)
# regarding the row count discrepancies, see the two versions of the SQL in 
# W:\PYWELL_SHARED\Pywell Projects\BRC\_BRC_projects\NPMS\Analyses\2018 08 - Per species trend analyses\SQL\extractData_v0.0.sql

## END