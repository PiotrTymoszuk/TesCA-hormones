# Definition and validation of hormonal clusters of testicular carcinoma

# tools --------

  library(tidyverse)
  library(rlang)
  library(trafo)
  library(stringi)

  library(exda)
  library(microViz)
  library(clustTools)

  library(ggrepel)
  library(figur)
  library(ggtext)

  library(furrr)

  library(soucer)

  insert_head()

  select <- dplyr::select
  reduce <- purrr::reduce
  explore <- exda::explore
  set_rownames <- trafo::set_rownames
  map <- purrr::map
  rename <- dplyr::rename

  c('./tools/globals.R',
    './tools/functions.R') %>%
    source_all(message = TRUE, crash = TRUE)

# analysis globals --------

  insert_msg('Analysis globals')

  c('./clustering scripts/globals.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Development of the clusters, predictions and evaluation -------

  insert_msg('Development of the clusters, predictions and evaluation')

  ## tuning of the clustering solution in the TCGA cohort

  access_cache(cache_path = './cache/clust_dev.RData',
               script_path = './clustering scripts/development.R',
               message = 'Cached cluster tuning results')

  c('./clustering scripts/predictions.R',
    './clustering scripts/evaluation.R',
    './/clustering scripts/features.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END -------

  insert_tail()
