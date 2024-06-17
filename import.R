# Import of clinical, transcriptome, mutation and CNA data.
# ComBat adjustment of log2 expression values - those data sets will be used for
# clustering.
# Computation of infiltration estimates (QunaTIseq, xCell, and MCP Counter),
# ssGSEA scores for Reactome pathways gene signatures, and Recon metabolic
# subsystem gene signatures.

# tools ------

  library(tidyverse)
  library(trafo)
  library(rlang)
  library(stringi)
  library(microViz)

  library(gseaTools)
  library(htGLMNET)
  library(biggrExtra)
  library(immunedeconv)

  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GEOquery)

  library(furrr)

  library(soucer)

  insert_head()

  select <- dplyr::select
  reduce <- purrr::reduce

  c('./tools/globals.R',
    './tools/functions.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Import scripts -------

  insert_msg('Import scripts')

  list(cache_path = c('./data/tcga.RData',
                      './data/gse3218.RData',
                      './data/gse99420.RData'),
       script_path = c('./import scripts/tcga.R',
                       './import scripts/gse3213.R',
                       './import scripts/gse99420.R'),
       message = paste('Loading cached',
                       c('TCGA data sets',
                         'GSE3218 data sets',
                         'GSE99420 data sets'))) %>%
    pwalk(access_cache)

# Immunodeconvolution and ssGSEA scores -------

  insert_msg('Immunodeconvolution and ssGSEA scores')

  ## infiltration

  list(cache_path = c('./data/quantiseq.RData',
                      './data/xcell.RData',
                      './data/mcp.RData'),
       message = paste('Loading cached',
                       c('QuanTIseq results',
                         'xCell results',
                         'MCP results'))) %>%
    pmap(access_cache,
         script_path = './import scripts/infiltration.R')

  ## ssGSEA scores

  list(cache_path = c('./data/recon.RData',
                      './data/reactome.RData'),
       script_path = c('./import scripts/recon.R',
                       './import scripts/reactome.R'),
       message = paste('Loading cached ssGSEA scores for',
                       c('Recon metabolic subsystems',
                         'Reatome pathways'))) %>%
    pwalk(access_cache)

# ComBat ------

  insert_msg('Batch adjustment')

  access_cache(cache_path = './data/combat.RData',
               script_path = './import scripts/combat.R',
               message = 'Loading cached ComBat-adjusted expression data')

# Prediction of drug sensitivity ------

  insert_msg('Prediction of drug sensitivity')

  access_cache(cache_path = './data/drugs.RData',
               script_path = './import scripts/drugs.R',
               message = 'Loading cached drug sensitivity predictions')

  ## merging the lexicons and predictions or the GDSC1 and GDSC2 experiments:
  ## they will be analyzed jointly, manual polishing of the molecular targets

  drugs$lexicons$gdsc <- drugs$lexicons[c("gdsc1", "gdsc2")] %>%
    reduce(rbind)

  drugs$lexicons <- drugs$lexicons[c("ctrp2", "gdsc")] %>%
    map(mutate,
        targets = map(targets,
                      stri_replace,
                      regex = '^\\s+',
                      replacement = ''),
        targets = map(targets,
                      stri_replace,
                      regex = '\\s+$',
                      replacement = ''),
        targets = map(targets,
                      stri_replace,
                      regex = '\\s+.*$',
                      replacement = ''))

  ## CTRP2: padding th missing drug names (combinations!)

  drugs$lexicons$ctrp2 <- drugs$lexicons$ctrp2 %>%
    mutate(drug_name = ifelse(!is.na(drug_name),
                              drug_name,
                              stri_replace(variable,
                                           regex = '\\s{1}\\(.*$',
                                           replacement = '')),
           drug_name = stri_capitalize(drug_name))

  drugs$predictions$gdsc <-
    map2(drugs$predictions$gdsc1,
         drugs$predictions$gdsc2,
         left_join, by = 'sample_id')

  drugs$predictions <- drugs$predictions[c("ctrp2", "gdsc")]

# END ------

  insert_tail()
