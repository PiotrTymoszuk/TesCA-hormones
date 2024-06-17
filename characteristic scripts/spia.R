# Analysis of modulation of signaling pathways (KEGG) by SPIA algorithm.
# The procedure is fed with log2-fold regulation estimates of expression
# as compared with the cohort mean obtained for the differentially regulated
# genes

  insert_head()

# container -------

  bcg_spia <- list()

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis globals ------

  insert_msg('Analysis globals')

  ## vectors with all investigated Entrez IDs

  bcg_spia$all_vectors <- bcg_dge$dev_test %>%
    map(filter,
        !is.na(entrez_id),
        entrez_id != '',
        !duplicated(entrez_id)) %>%
    map(~.x$entrez_id) %>%
    map(unname)

  ## vectors with log2 fold-regulation estimates
  ## for the differentially regulated genes named with Entrez IDs

  bcg_spia$de_vectors <-  bcg_dge$dev_test %>%
    map(blast, clust_id) %>%
    map(map,
        filter,
        !is.na(entrez_id),
        entrez_id != '',
        !duplicated(entrez_id)) %>%
    map(map,
        ~set_names(.x$deviation_center, .x$entrez_id))

# Perturbation and enrichment analysis --------

  insert_msg('Perturbation and enrichment analysis')

  for(i in names(bcg_spia$de_vectors)) {

    bcg_spia$test[[i]] <- bcg_spia$de_vectors[[i]] %>%
      future_map(spia,
                 all = bcg_spia$all_vectors[[i]],
                 verbose = FALSE,
                 .options = furrr_options(seed = TRUE)) %>%
      map(as_tibble)

  }

# Formatting the analysis results -------

  insert_msg('Formatting')

  bcg_spia$test <- bcg_spia$test %>%
    transpose %>%
    map(map,
        mutate,
        regulation = ifelse(pGFdr >= 0.05, 'ns',
                            ifelse(tA > 0, 'activated',
                                   ifelse(tA < 0,
                                          'inhibited', 'ns'))),
        regulation = factor(regulation, c('activated', 'inhibited', 'ns')))

# Significant effects ------

  insert_msg('Significant effects')

  ## in single cohorts

  bcg_spia$significant <- bcg_spia$test %>%
    map(map,
        filter,
        regulation %in% c('activated', 'inhibited')) %>%
    map(map, blast, regulation) %>%
    map(compact) %>%
    map(transpose) %>%
    map(map, map, ~.x$Name)

  ## common regulated pathways: TCGA and GSE99420

  bcg_spia$common_significant <- bcg_spia$significant %>%
    map(map, ~.x[c('tcga', 'gse99420')]) %>%
    map(map, reduce, intersect)

# Caching the results -------

  insert_msg('Caching')

  bcg_spia <- bcg_spia[c('test', 'significant', 'common_significant')]

  save(bcg_spia, file = './cache/bcg_spia.RData')

# END ------

  plan('sequential')

  insert_tail()
