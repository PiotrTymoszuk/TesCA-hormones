# Characteristic of the hormonal clusters

# tools --------

  library(tidyverse)
  library(rlang)
  library(trafo)
  library(stringi)
  library(rvest)

  library(exda)
  library(microViz)
  library(clustTools)
  library(GOSemSim)
  library(decoupleR)
  library(SPIA)

  library(BiGGR)
  library(biggrExtra)
  library(MatrisomeAnalyzeR)

  library(org.Hs.eg.db)
  library(AnnotationDbi)

  library(survival)
  library(survminer)
  library(coxExtensions)
  library(glmnet)
  library(caret)

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
  components <- generics::components
  transpose <- purrr::transpose

  c('./tools/globals.R',
    './tools/functions.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Analysis globals --------

  insert_msg('Analysis globals')

  c('./characteristic scripts/globals.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Clinical characteristic, survival and relapse rate --------

  insert_msg('Clinical characteristic')

  c('./characteristic scripts/clinic.R',
    './characteristic scripts/survival.R',
    './characteristic scripts/relapse.R',
    './characteristic scripts/survival_modeling.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Infiltration, GSVA and differential gene expression: transcriptome analyses ---------

  insert_msg('Infiltration, GSVA, and differential gene expression')

  c('./characteristic scripts/infiltration.R',
    './characteristic scripts/reactome.R',
    './characteristic scripts/recon.R') %>%
    source_all(message = TRUE, crash = TRUE)

  ## differential gene expression, GO term enrichment,
  ## and activity of metabolic reactions:
  ## resorting to cached results, because computation takes a while

  list(cache_path = c('./cache/bcg_dge.RData',
                      './cache/bcg_go.RData',
                      './cache/bcg_goclust.RData',
                      './cache/bcg_meta.RData',
                      './cache/bcg_spia.RData'),
       script_path = c('./characteristic scripts/dge.R',
                       './characteristic scripts/go.R',
                       './characteristic scripts/go_clustering.R',
                       './characteristic scripts/metabolism.R',
                       './characteristic scripts/spia.R'),
       message = paste('Loading cached',
                       c('differential expression analysis results',
                         'GO enrichment analysis results',
                         'GO clustering results',
                         'activity of metabolic reactions',
                         'signaling by SPIA'))) %>%
    pwalk(access_cache)

  ## distances between the hormonal clusters in respect to the
  ## differential gene expression

  c('./characteristic scripts/dge_distance.R') %>%
    source_all(message = TRUE, crash = TRUE)

  ## plots for the differential gene expression analysis, GO enrichment and
  ## clustering, metabolism, and SPIA signaling pathways

  c('./characteristic scripts/dge_plots.R',
    './characteristic scripts/go_plots.R',
    './characteristic scripts/go_cluster_plots.R',
    './characteristic scripts/subsystem_plots.R',
    './characteristic scripts/metabolism_plots.R',
    './characteristic scripts/spia_plots.R') %>%
    source_all(message = TRUE, crash = TRUE)

  ## analysis of collecTRI regulons and PROGENy signaling pathways
  ## with the decopuleR tools

  c('./characteristic scripts/collectri.R',
    './characteristic scripts/progeny.R') %>%
    source_all(message = TRUE, crash = TRUE)

  ## detailed analysis of selected parts of transcriptome:
  ## genes associated  with immune checkpoint and T cell exhaustion
  ## genes of testis and cancer antigens (CTA)
  ## matrisome genes
  ## steroid hormone receptor genes
  ## genes of ERBB and FGFR families and theis ligands

  c('./characteristic scripts/immune_checkpoint.R',
    './characteristic scripts/antigens.R',
    './characteristic scripts/matrisome.R',
    './characteristic scripts/hormone_receptors.R',
    './characteristic scripts/errb_fgfr.R',
    './characteristic scripts/estrogen_responsive_genes.R',
    './characteristic scripts/androgen_responsive_genes.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Genetic background of the clusters and differential protein expression -------

  insert_msg('Genetics and protein expression')

  c('./characteristic scripts/protein.R',
    './characteristic scripts/protein_networks.R',
    './characteristic scripts/genetics.R',
    './characteristic scripts/genetic_numbers.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Differences in drug sensitivity and treatment between the clusters -------

  insert_msg('Drug sensitivity and treatment in the clusters')

  c('./characteristic scripts/drugs.R',
    './characteristic scripts/drug_targets.R',
    './characteristic scripts/drug_plots.R',
    './characteristic scripts/approved_drugs.R',
    './characteristic scripts/selected drugs.R',
    './characteristic scripts/therapy_type.R',
    './characteristic scripts/therapy_response.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END -------

  insert_tail()
