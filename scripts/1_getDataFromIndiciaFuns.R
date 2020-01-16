## 1. Functions to retrieve plot/sample/habitat/level data from Indicia/NPMS database (to avoid hardcoding passwords in script)
# O.L. Pescott
# 29.08.2019

# SQL query is also in W:\PYWELL_SHARED\Pywell Projects\BRC\_BRC_projects\NPMS\Analyses\2018 08 - Per species trend analyses\SQL\extractData_v0.0.sql

####
#### Warning! This is slow (11 mins if retrieving data up between 2015-mid-2018) and will get increasingly slower
#### Need to investigate possible speed-ups/new cache table that runs this automatically on some periodic basis
####
list.of.packages <- c("RPostgreSQL", "rstudioapi")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(RPostgreSQL)
library(rstudioapi)

drv <- dbDriver("PostgreSQL")
pwd <- rstudioapi::askForPassword("Database password")
# Create a connection to the postgres database
Connection <- dbConnect(drv,
                        dbname = "warehouse1",
                        host = "192.171.199.208", port = 5432,
                        user = "brc_read_only", password = pwd
)
############# Function to retrieve plot and sample level attributes
## Change dates as appropriate! (Lines 23 and 36) ##
# dsn = "PostgreSQL30"; user = "brc_read_only"
## v1.1 incorporates grazing and management data
getNpmsData_PlotsSamples_v1.1 <- function(Connection = Connection) {
                                                                    select_query <- paste(
                                                                    "select sq.plot_id, sq.monad, occs.sample, occs.title, habs.caption, habs.term as surv_habitat",
                                                                    "from",
	                                                                    "(select s.id as sample, su.title",
                                                                      "from indicia.samples s",
                                                                      "left join indicia.cache_occurrences co",
                                                                      "on co.sample_id = s.id",
                                                                      "left join indicia.surveys su",
                                                                      "on su.id = s.survey_id",
                                                                      "where s.date_start between '2015-01-01' and '2019-12-31'",
                                                                      "and s.survey_id in (87,154,155)",
                                                                      "group by s.id, su.title",
                                                                      ") occs",
                                                                    "join",
                                                                      "(select s.id sample_id, sa.caption, t.term",
                                                                      "from indicia.samples s, indicia.sample_attributes sa, indicia.sample_attribute_values sav, indicia.termlists_terms tt, indicia.terms t",
                                                                      "where s.id = sav.sample_id",
                                                                      "and sav.sample_attribute_id = sa.id",
                                                                      "and tt.termlist_id = sa.termlist_id",
                                                                      "and tt.id = sav.int_value",
                                                                      "and t.id = tt.term_id",
                                                                      "and s.survey_id in (87,154,155)",
                                                                      "and s.date_start between '2015-01-01' and '2019-12-31'",
                                                                      "and int_value is not null",
                                                                      "and sa.termlist_id is not null",
                                                                      "and sa.caption in ('NPMS Habitat', 'NPMS Grazing', 'NPMS Management')",
                                                                      "order by sample_attribute_id, sample_id",
                                                                      ") habs",
                                                                    "on occs.sample = habs.sample_id",
                                                                    "join",
                                                                      "(select distinct square.centroid_sref as monad, plot.id as plot_id, visit.id as sample_id",
                                                                      "from indicia.samples visit",
                                                                      "join indicia.locations plot on plot.id = visit.location_id AND plot.deleted=false",
                                                                      "join indicia.locations square on square.id = plot.parent_id AND square.deleted=false",
                                                                      "where visit.survey_id in (87,154,155) AND visit.deleted=false",
                                                                      "group by square.centroid_sref, plot.id, visit.id",
                                                                      ") sq",
                                                                    "on occs.sample = sq.sample_id"
                                                                    )
                                                                    temp <- dbGetQuery(conn = Connection, statement = select_query)
                                                                    return(temp) }

############# Function to retrieve taxon data for samples
## Change dates as appropriate! (Lines 67 and 76)
## 03.09.2019 -- update to exclude rejected records
getNpmsData_SamplesSpecies <- function(Connection = Connection) {
                                                                query1 <- paste( ## Wildflower data
                                                                  "select co.id, co.sample_id, co.date_start as date, co.preferred_taxon, co.taxa_taxon_list_external_key as tvk,  oav.int_value as domin, co.record_status as record_status, co.sensitivity_precision",
                                                                  "from indicia.cache_occurrences co",
                                                                  "LEFT join indicia.occurrence_attribute_values oav on oav.occurrence_id = co.id and oav.deleted=false",
                                                                  "LEFT join indicia.occurrence_attributes oa on oa.id=oav.occurrence_attribute_id",
                                                                  "AND oa.id = 104 and oav.deleted=false Where co.survey_id in (87) and co.training = 'f'",
                                                                  "and co.date_start between '2015-01-01' and '2019-12-31'",
                                                                  "and co.record_status in ('V','C')", ## exclude rejected data
                                                                  "group by co.id, co.sample_id, co.date_start, co.preferred_taxon, co.taxa_taxon_list_external_key, oav.int_value, co.record_status, co.sensitivity_precision;"
                                                                  )
                                                                temp1 <- dbGetQuery(conn = Connection, statement = query1)
                                                                query2 <- paste( ## Indicator and Inventory data
                                                                  "select co.id, co.sample_id, co.date_start as date, co.preferred_taxon, co.taxa_taxon_list_external_key as tvk,  terms.term as domin, co.record_status as record_status, co.sensitivity_precision",
                                                                  "from indicia.cache_occurrences co",
                                                                  "LEFT join indicia.occurrence_attribute_values oav on oav.occurrence_id = co.id  and oav.deleted=false",
                                                                  "LEFT join indicia.occurrence_attributes oa on oa.id=oav.occurrence_attribute_id AND oa.id = 214 and oav.deleted=false",
                                                                  "LEFT join indicia.termlists_terms tt on tt.id=oav.int_value and tt.deleted=false",
                                                                  "LEFT join indicia.terms terms on terms.id = tt.term_id AND terms.deleted=false Where co.survey_id in (154,155) and co.training = 'f'",
                                                                  "and co.date_start between '2015-01-01' and '2019-12-31'",
                                                                  "and co.record_status in ('V','C')", ## exclude rejected data
                                                                  "group by co.id, co.sample_id, co.date_start, co.preferred_taxon, co.taxa_taxon_list_external_key, terms.term, co.record_status, co.sensitivity_precision;"
                                                                  )
                                                                temp2 <- dbGetQuery(conn = Connection, statement = query2)
                                                                temp3 <- rbind(temp1, temp2[-1,]) ## rbind but remove header from second file
                                                                return(temp3) }
