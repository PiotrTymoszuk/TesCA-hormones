# Differential regulation of activity of metabolic reactions.
# Enrichment of significantly activated and inhibited metabolic reactions
# in the Recon subsystems.
#
# The analysis is done with tools provided by biggR and biggrExtra.
# In Monte Carlo modeling, regulation estimates (cluster versus cohort mean)
# with standard error for all available genes are used.

  insert_head()

# container -------

  bcg_meta <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  ## log2 fold-regulation estimates of differential expression
  ## with SEM

  bcg_meta$dge <- bcg_dge$dev_test %>%
    map(select,
        clust_id,
        gene_symbol,
        entrez_id,
        deviation_center,
        deviation_sem) %>%
    map(mutate,
        entrez_id = unname(entrez_id)) %>%
    map(blast, clust_id)

# Modeling --------

  insert_msg('Monte Carlo modeling')

  for(i in names(bcg_meta$dge)) {

    bcg_meta$models[[i]] <- bcg_meta$dge[[i]] %>%
      map(~build_geneSBML(x = set_names(.x$deviation_center,
                                        .x$entrez_id),
                          err = set_names(.x$deviation_sem,
                                          .x$entrez_id),
                          database = Recon2D,
                          scale = 'log2',
                          or_fun = 'mean',
                          and_fun = 'min',
                          x_default = 1,
                          err_method = 'mc',
                          n_iter = 3000,
                          ci_method = 'bca',
                          save_memory = TRUE,
                          .parallel = TRUE))

  }

  bcg_meta$models <- transpose(bcg_meta$models)

# Reaction regulation estimates -------

  insert_msg('Reaction regulation estimates, significant reactions')

  bcg_meta$estimates <- bcg_meta$models %>%
    map(map, components, 'regulation') %>%
    map(map,
        mutate,
        regulation = ifelse(p_adjusted >= 0.05, 'ns',
                            ifelse(fold_reg < 1, 'inhibited',
                                   ifelse(fold_reg > 0,
                                          'activated', 'ns'))),
        regulation = factor(regulation, c('activated', 'inhibited', 'ns')))

  ## significantly regulated reactions
  ## common significant: reactions shared by the TCGA and GSE99420 cohorts

  bcg_meta$significant <- bcg_meta$estimates %>%
    map(map, filter, regulation %in% c('activated', 'inhibited')) %>%
    map(map, blast, regulation) %>%
    map(transpose) %>%
    map(map, map, ~.x$react_id)

  bcg_meta$common_significant <- bcg_meta$significant %>%
    map(map, ~.x[c('tcga', 'gse99420')]) %>%
    map(map, reduce, intersect)

# Numbers of differentially regulated reactions ------

  insert_msg('Numbers of differenctially regulated reactcions')

  bcg_meta$react_numbers <- bcg_meta$models %>%
    map(map, count) %>%
    map(map,
        filter,
        subsystem == 'All reactions',
        status %in% c('activated', 'inhibited')) %>%
    map(map, select, -subsystem) %>%
    map(compress, names_to = 'cohort') %>%
    compress(names_to = 'clust_id')

# Subsystem enrichment analysis ------

  insert_msg('Subsystem enrichment analysis')

  bcg_meta$enrichment <- bcg_meta$models %>%
    map(map,
        suba,
        signif_type = 'fdr',
        method = 'simulation',
        n_iter = 100000,
        .parallel = TRUE)

# Significant metabolic subsystems ------

  insert_msg('Significantly enriched subsystems')

  ## significant subsystem enrichment in single cohort,
  ## I'm working with unadjusted p values, because pFFDR depends on
  ## iteration numbers

  bcg_meta$enrichment_significant <- bcg_meta$enrichment %>%
    map(map,
        filter,
        p_value < 0.05,
        status %in% c('activated', 'inhibited')) %>%
    map(map, blast, status) %>%
    map(transpose) %>%
    map(map, map, ~.x$subsystem) %>%
    map(map, map, as.character)

  ## common significantly enriched effects shared by
  ## the TCGA and GSE9940 cohort

  bcg_meta$enrichment_common_significant <-
    bcg_meta$enrichment_significant %>%
    map(map, ~.x[c('tcga', 'gse99420')]) %>%
    map(map, reduce, intersect)

# Caching the results ------

  insert_msg('Caching')

  bcg_meta$dge <- NULL

  bcg_meta <- compact(bcg_meta)

  save(bcg_meta, file = './cache/bcg_meta.RData')

# END -------

  rm(i)

  insert_tail()
