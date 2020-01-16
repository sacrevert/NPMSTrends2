list.of.packages <- c("RODBC")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(RODBC)
############# Function to retrieve plot and sample level attributes
## Change dates as appropriate! (Lines 23 and 36) ##
# dsn = "PostgreSQL30"; user = "brc_read_only"
## v1.1 incorporates grazing and management data
getNpmsData_PlotsSamples_v1.1 <- function(dsn = "PostgreSQL30", user = "brc_read_only", password) {   
  channel <- RODBC::odbcConnect(dsn = dsn, uid = user, pwd = password)
  temp <- RODBC::sqlQuery(channel, paste(
    "select sq.plot_id, sq.monad, occs.sample, occs.title, habs.caption, habs.term as surv_habitat",
    "from",
    "(select s.id as sample, su.title",
    "from indicia.samples s",
    "left join indicia.cache_occurrences co",
    "on co.sample_id = s.id",
    "left join indicia.surveys su",
    "on su.id = s.survey_id",
    "where s.date_start between '2019-01-01' and '2019-12-31'",
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
    "and s.date_start between '2019-01-01' and '2019-12-31'",
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
    "on occs.sample = sq.sample_id"))
  RODBC::odbcClose(channel)
  rm(channel)
  return(temp) }

npms_plots <- getNpmsData_PlotsSamples_v1.1(password = "ozpc0829ULBu#E4eVW20")
nrow(npms_plots)