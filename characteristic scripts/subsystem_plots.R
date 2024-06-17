# Plots for the results of metabolic subsystem enrichment analysis.
# Metabolic subsystems found to be significantly enriched with activated or
# or inhibited reactions in the TCGA and GSE99420 cohorts are presented

  insert_head()

# container -------

  bcg_subplots <- list()

# analysis data ------

  insert_msg('Analysis data')

  bcg_subplots[c("data",
                 "common_significant")] <-
    bcg_meta[c("enrichment",
               "enrichment_common_significant")]

  ## enrichment data: re-formatting for plotting

  bcg_subplots$variables <- bcg_subplots$common_significant %>%
    unlist %>%
    unname %>%
    unique

  bcg_subplots$data <- bcg_subplots$data %>%
    transpose %>%
    map(compress,
        names_to = 'clust_id') %>%
    map(filter,
        subsystem %in% bcg_subplots$variables,
        status %in% c('activated', 'inhibited')) %>%
    map(mutate,
        clust_id = factor(clust_id,
                          levels(bcg_globals$assignment[[1]]$clust_id)),
        regulation = ifelse(p_value >= 0.05, 'ns', status),
        regulation = factor(regulation, c('activated', 'inhibited', 'ns')))

  ## metabolic subsystem classification per hand

  bcg_subplots$classification <-
    c("Arachidonic acid metabolism" = 'lipids',
      "Heparan sulfate degradation" = 'ECM/glycoproteins',
      "Inositol phosphate metabolism" = 'energy',
      "Nucleotide interconversion" = 'nucleotides',
      'Vitamin D metabolism' = 'vitamins',
      'Xenobiotics metabolism' = 'detoxification',
      'Purine synthesis' = 'nucleotides',
      'Transport, extracellular' = 'transport',
      'Folate metabolism' = 'vitamins',
      'Transport, golgi apparatus' = 'transport',
      'Vitamin E metabolism' = 'vitamins',
      'Blood group synthesis' = 'ECM/glycoproteins',
      'Chondroitin synthesis' = 'ECM/glycoproteins',
      'Keratan sulfate degradation' = 'ECM/glycoproteins',
      'Transport, peroxisomal' = 'transport',
      'Androgen and estrogen synthesis and metabolism' = 'steroids',
      'Fatty acid synthesis' = 'lipids',
      'Methionine and cysteine metabolism' = 'aminio acids',
      'Aminosugar metabolism' = 'ECM/glycoproteins',
      'Citric acid cycle' = 'energy',
      'Oxidative phosphorylation' = 'energy',
      'ROS detoxification' = 'detoxification',
      'Sphingolipid metabolism' = 'lipids',
      'Starch and sucrose metabolism' = 'energy',
      'Transport, lysosomal' = 'transport',
      'Vitamin B6 metabolism' = 'vitamins',
      'Eicosanoid metabolism' = 'lipids',
      'Nucleotide salvage pathway' = 'nucleotides',
      'Pyruvate metabolism' = 'energy',
      'Steroid metabolism' = 'steroids',
      'Cytochrome metabolism' = 'energy',
      'Fatty acid oxidation' = 'energy') %>%
    compress(names_to = 'subsystem',
             values_to = 'class')

  bcg_subplots$data <- bcg_subplots$data %>%
    map(left_join,
        bcg_subplots$classification,
        by = 'subsystem')

# Bubble plots for particular cohorts --------

  insert_msg('Bubble plots for single cohorts')

  bcg_subplots$bubble_plots <-
    list(x = bcg_subplots$data,
         y = globals$cohort_labs[names(bcg_subplots$data)]) %>%
    pmap(function(x, y) x %>%
           ggplot(aes(x = clust_id,
                      y = subsystem,
                      fill = regulation,
                      size = OR)) +
           geom_point(shape = 21) +
           facet_grid(class ~ status,
                      scales = 'free',
                      space = 'free',
                      labeller = labeller(.cols = c(activated = 'enriched\nwith activated\nreactions',
                                                    inhibited = 'enriched\nwith inhibited\nreactions'))) +
           scale_y_discrete(labels = subsystem_labeller) +
           scale_fill_manual(values = c(activated = 'firebrick',
                                        inhibited = 'steelblue',
                                        ns = 'gray70'),
                             labels = c(activated = 'activated\nreactions',
                                        inhibited = 'inhibited\nreactions',
                                        ns = 'ns'),
                             name = 'enrichment') +
           globals$common_theme +
           theme(axis.title.y = element_blank()) +
           labs(title = y,
                subtitle = 'Shared significant Recon metabolic subsystems',
                x = 'cluster'))

# Bar plots with the enrichment ORs -------

  insert_msg('Bar plots with the enrichment ORs')

  ## done for common significantly enriched metabolic subsystems
  ## in the TCGA and GSE99420 cohorts

  bcg_subplots$common_significant <-
    bcg_subplots$common_significant %>%
    map(reduce, union)

  bcg_subplots$common_data <- bcg_subplots$data[c("tcga", "gse99420")] %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = factor(cohort, c("tcga", "gse99420"))) %>%
    blast(clust_id) %>%
    map2(., bcg_subplots$common_significant,
         ~filter(.x,
                 subsystem %in% .y,
                 regulation %in% c('activated', 'inhibited'))) %>%
    map(select,
        cohort,
        subsystem,
        status,
        class,
        OR)

  ## bar plots

  bcg_subplots$bar_plots <-
    list(x = bcg_subplots$common_data,
         y = paste('RECON metabolic subsystems, cluster',
                   names(bcg_subplots$common_data))) %>%
    pmap(function(x, y) x %>%
           ggplot(aes(x = OR,
                      y = reorder(subsystem, OR),
                      fill = status)) +
           facet_grid(status ~ cohort,
                      scales = 'free',
                      space = 'free_y',
                      labeller = labeller(.cols = globals$cohort_labs)) +
           geom_vline(xintercept = 0,
                      linetype = 'dashed') +
           geom_bar(stat = 'identity',
                    color = 'black') +
           scale_fill_manual(values = c(activated = 'firebrick',
                                        inhibited = 'steelblue',
                                        ns = 'gray60'),
                             labels = c(activated = 'activated reactions',
                                        inhibited = 'inhibited reactions'),
                             name = 'Endichment with') +
           scale_y_discrete(labels = subsystem_labeller) +
           globals$common_theme +
           theme(axis.title.y = element_blank()) +
           labs(title = y,
                subtitle = 'Common significantly enriched subsystems',
                x = 'enrichment OR'))

# END -----

  bcg_subplots$data <- NULL
  bcg_subplots$common_significant <- NULL
  bcg_subplots$common_data <- NULL

  bcg_subplots <- compact(bcg_subplots)

  insert_tail()
