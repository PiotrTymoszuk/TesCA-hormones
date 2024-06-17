# Exploratory data analysis:
#
# 1) Characteristic of the study cohorts
#
# 2) Comparison of expression of the hormone-related genes between
# cancer histologies
#
# 3) Analysis of gene co-expression by pairwise correlations and PCA.
# Assessment of spontaneous clustering tendency by Hopkins stat.
#
# 4) Identification of hubs of cancer sex hormone metabolism by a network
# analysis of co-expression graphs.


# tools -------

  library(tidyverse)
  library(rlang)
  library(trafo)
  library(stringi)

  library(exda)
  library(microViz)
  library(clustTools)

  library(igraph)
  library(ggnetwork)

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

  c('./tools/globals.R',
    './tools/functions.R') %>%
    source_all(message = TRUE, crash = TRUE)

# analysis scripts -------

  insert_msg('Analysis scripts')

  ## characteristic of the study cohorts,
  ## analysis of distribution of expression of the hormone-related genes
  ## differences in clinical variables and expression of the hormone-related
  ## genes between the cohorts

  c('./exploration scripts/cohorts.R',
    './exploration scripts/distribution.R',
    './exploration scripts/expression.R') %>%
    source_all(message = TRUE, crash = TRUE)

  ## analysis of differential and co-expression of the hormone-related genes
  ## variables with sufficient variability and expression levels identified by
  ## the './exploration scripts/distribution.R' script are used in the analyses

  c('./exploration scripts/histology.R',
    './exploration scripts/histology_icd.R',
    './exploration scripts/co_expression.R',
    './exploration scripts/pca.R',
    './exploration scripts/network.R') %>%
    source_all(message = TRUE, crash = TRUE)

  ## general frequency of somatic mutations, gene deletions and amplifications
  ## distribution tests for drug sensitivity estimates

  c('./exploration scripts/genetics.R',
    './exploration scripts/drugs.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END -----

  insert_tail()
