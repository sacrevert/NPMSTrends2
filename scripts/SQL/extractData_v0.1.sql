-- ### Extract species data with appropriate 1/0/NA values, as befitting survey and species ### --
-- v.0.1
-- 05 09 2019
-- O.L. Pescott, olipes@ceh.ac.uk

--  ### Table of samples with habitat and level information ### --
-- can this be a cache-built table?
-----------------------------------
-- ## file/function for retrieval: W:\PYWELL_SHARED\Pywell Projects\BRC\_BRC_projects\NPMS\Analyses\2019 03 -- Trend Analyses 2
-----------------------------------
select sq.plot_id, sq.monad, occs.sample, occs.title, habs.caption, habs.term as surv_habitat from
	                                                                  (select s.id as sample, su.title
                                                                      from indicia.samples s
                                                                      left join indicia.cache_occurrences co
                                                                      on co.sample_id = s.id
                                                                      left join indicia.surveys su
                                                                      on su.id = s.survey_id
                                                                      where s.date_start between '2018-01-01' and '2018-12-31'
                                                                      and s.survey_id in (87,154,155)
                                                                      group by s.id, su.title
                                                                      ) occs
                                                                    join
                                                                      (select s.id sample_id, sa.caption, t.term
                                                                      from indicia.samples s, indicia.sample_attributes sa, indicia.sample_attribute_values sav, indicia.termlists_terms tt, indicia.terms t
                                                                      where s.id = sav.sample_id
                                                                      and sav.sample_attribute_id = sa.id
                                                                      and tt.termlist_id = sa.termlist_id
                                                                      and tt.id = sav.int_value
                                                                      and t.id = tt.term_id
                                                                      and s.survey_id in (87,154,155)
                                                                      and s.date_start between '2018-01-01' and '2018-12-31'
                                                                      and int_value is not null
                                                                      and sa.termlist_id is not null
                                                                      and sa.caption in ('NPMS Habitat', 'NPMS Grazing', 'NPMS Management')
                                                                      order by sample_attribute_id, sample_id
                                                                      ) habs
                                                                    on occs.sample = habs.sample_id
                                                                    join
                                                                      (select distinct square.centroid_sref as monad, plot.id as plot_id, visit.id as sample_id
                                                                      from indicia.samples visit
                                                                      join indicia.locations plot on plot.id = visit.location_id AND plot.deleted=false
                                                                      join indicia.locations square on square.id = plot.parent_id AND square.deleted=false
                                                                      where visit.survey_id in (87,154,155) AND visit.deleted=false
                                                                      group by square.centroid_sref, plot.id, visit.id
                                                                      ) sq
                                                                    on occs.sample = sq.sample_id;
																	
--## Species data ##--
-- Wildflower
select co.id, co.sample_id, co.date_start as date, co.preferred_taxon, co.taxa_taxon_list_external_key as tvk,  oav.int_value as domin, co.record_status as record_status, co.sensitivity_precision
                                                                                    from indicia.cache_occurrences co
                                                                                    left join indicia.occurrence_attribute_values oav on oav.occurrence_id = co.id and oav.deleted=false
                                                                                    left join indicia.occurrence_attributes oa on oa.id=oav.occurrence_attribute_id
																					left join indicia.
                                                                                    and oa.id = 104 and oav.deleted=false Where co.survey_id in (87) and co.training = 'f'
                                                                                    and co.date_start between '2018-01-01' and '2018-12-31'
                                                                                    and co.record_status in ('V','C')
                                                                                    group by co.id, co.sample_id, co.date_start, co.preferred_taxon, co.taxa_taxon_list_external_key, oav.int_value, co.record_status, co.sensitivity_precision
																					limit 10;
-- Indicator and Inventory
select co.id, co.sample_id, co.date_start as date, co.preferred_taxon, co.taxa_taxon_list_external_key as tvk,  terms.term as domin, co.record_status as record_status, co.sensitivity_precision
                                                                                    from indicia.cache_occurrences co
                                                                                    LEFT join indicia.occurrence_attribute_values oav on oav.occurrence_id = co.id  and oav.deleted=false
                                                                                    LEFT join indicia.occurrence_attributes oa on oa.id=oav.occurrence_attribute_id AND oa.id = 214 and oav.deleted=false
                                                                                    LEFT join indicia.termlists_terms tt on tt.id=oav.int_value and tt.deleted=false
                                                                                    LEFT join indicia.terms terms on terms.id = tt.term_id AND terms.deleted=false Where co.survey_id in (154,155) and co.training = 'f'
                                                                                    and co.date_start between '2015-01-01' and '2018-12-31'
                                                                                    and co.record_status in ('V','C')
                                                                                    group by co.id, co.sample_id, co.date_start, co.preferred_taxon, co.taxa_taxon_list_external_key, terms.term, co.record_status, co.sensitivity_precision;