# Differences in ssGSEA scores of the Reactome pathway gene signatures
# between the clusters.
#
# Statistica significance is determined by one-way ANOVA with eta-square effect
# size statistic.
#
# Significant effects are defined by pFDR < 0.05 and eta-square >= 0.14.

  insert_head()

# container ------

  bcg_reactome <- list()

# Parallel backend ------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis data and variables ------

  insert_msg('Analysis data and variables')

  bcg_reactome$data <-
    map2(bcg_globals$assignment[names(reactome$scores)],
         reactome$scores,
         inner_join, by = 'sample_id')

  bcg_reactome$lexicon <- reactome$lexicon

# N numbers -------

  insert_msg('N numbers')

  ## ready-to-use plot captions

  bcg_reactome$n_captions <- bcg_reactome$data %>%
    map(count, clust_id) %>%
    map(~map2_chr(.x[[1]], .x[[2]],
                  paste, sep = ': n = ')) %>%
    map(paste, collapse = ', ')

# descriptive stats --------

  insert_msg('Descriptive stats')

  bcg_reactome$stats <- bcg_reactome$data %>%
    map(fast_num_stats,
        split_factor = 'clust_id')

# Differences between the clusters -------

  insert_msg('Testing')

  bcg_reactome$test <- bcg_reactome$data %>%
    future_map(test_anova,
               variables = bcg_reactome$lexicon$variable,
               split_fct = 'clust_id',
               .parallel = FALSE,
               .options = furrr_options(seed = TRUE)) %>%
    map(~.x$anova)

# Formatting the results and significant effects -------

  insert_msg('Formatting the results and significant effects')

  bcg_reactome$test <- bcg_reactome$test %>%
    ## to get nicely formatted significance
    map(re_adjust, method = 'none') %>%
    map(mutate,
        variable = response,
        eff_size = paste('\u03B7\u00B2 =', signif(effect_size, 2)),
        ## axis labels to be used (optionally) at heat map Y axes
        ## significant effects are highlighted in bold
        axis_lab = exchange(variable,
                            bcg_reactome$lexicon),
        axis_lab = ifelse(p_adjusted < 0.05 & effect_size >= 0.14,
                          html_bold(axis_lab), axis_lab))

  ## significant in single cohorts
  ## and common significant effect shared by the TCGA and GSE99420 cohorts

  bcg_reactome$significant <- bcg_reactome$test %>%
    map(filter,
        p_adjusted < 0.05,
        effect_size >= 0.14) %>%
    map(~.x$variable)

  bcg_reactome$common_significant <-
    bcg_reactome$significant[c("tcga", "gse99420")] %>%
    reduce(intersect)

# Classification of the common regulated signatures, TCGA --------

  insert_msg('Classification of the common regulated signatures')

  ## done for the TCGA cohort and applied to the remaining collectives

  bcg_reactome$classification_tcga <- bcg_reactome$data$tcga %>%
    classify(variables = bcg_reactome$common_significant,
             split_fct = 'clust_id') %>%
    .$classification

  ## description of the signatures specific for particular clusters

  bcg_reactome$classification_descr <-
    list('#1' = c('FGFR/ERBB/PI3K',
                  'EPH/NOTCH',
                  'GPCR/Ca2+',
                  'ECM/keratin',
                  'fatty acids',
                  'ketone bodies'),
         '#2' = c('Ag-presentation',
                  'apoptosis',
                  'DNA repair',
                  'inflammation'),
         '#3' = c('lipid metabolism',
                  'PPARA/energy',
                  'FGFR/ERBB/PI3K/AKT',
                  'steroid hormones',
                  'estrogens'),
         '#4' = c('androgens',
                  'fatty acids',
                  'nucleotides',
                  'gluconeogenesis, TCA',
                  'cell cycle, mitosis',
                  'RNA turnover',
                  #'mitochondria',
                  'autophagy',
                  'apoptosis')) %>%
    map_chr(paste, collapse = '\n') %>%
    compress(names_to = 'clust_id',
             values_to = 'description')

  bcg_reactome$classification_tcga <-
    left_join(bcg_reactome$classification_tcga,
              bcg_reactome$classification_descr,
              by = 'clust_id')

# Heat maps of ssGSEA scores of the common regulated signatures -----

  insert_msg('Heat maps of ssGSEA scores')

  bcg_reactome$hm_plots <-
    list(data = bcg_reactome$data,
         plot_title = paste('Reactome pathways,',
                            globals$cohort_labs[names(bcg_reactome$data)])) %>%
    pmap(heat_map,
         variables = bcg_reactome$common_significant,
         split_fct = 'clust_id',
         normalize = FALSE,
         cust_theme = globals$common_theme +
           theme(axis.text.x = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks = element_blank(),
                 strip.text.y = element_text(angle = 0,
                                             size = 6.5,
                                             hjust = 0)),
         variable_classification = bcg_reactome$classification_tcga[, c("variable", "description")],
         limits = c(-1, 1),
         midpoint = 0,
         oob = scales::squish,
         x_lab = 'cancer sample',
         y_lab = 'Reactome pathway gene signature',
         fill_lab = 'ssGSEA\nscore')

# Distance between the clusters by the regulated signatures ------

  insert_msg('Distance between the clusters')

  ## cosine distances convertible to cosine similarities

  bcg_reactome$clust_distances <- bcg_reactome$data %>%
    map(column_to_rownames, 'sample_id') %>%
    map2(., bcg_reactome$significant,
         ~select(.x, clust_id, all_of(.y))) %>%
    map2(., bcg_reactome$significant,
         ~subset_distance(.x,
                          variables = .y,
                          split_fct = 'clust_id',
                          dist_FUN = calculate_dist,
                          method = 'cosine'))

  ## summary heat maps of the distances

  bcg_reactome$dist_heat_maps <-  bcg_reactome$clust_distances %>%
    map(plot,
        'mean',
        cust_theme = globals$common_theme)

  bcg_reactome$dist_heat_maps <-
    list(x = bcg_reactome$dist_heat_maps,
         y = globals$cohort_labs[names(bcg_reactome$dist_heat_maps)],
         z = bcg_reactome$n_captions) %>%
    pmap(function(x, y, z) x +
           labs(title = y,
                subtitle = z,
                x = paste(y, 'clusters'),
                y = paste(y, 'clusters')) +
           scale_fill_gradient2(low = 'firebrick',
                                mid = 'white',
                                high = 'steelblue',
                                limits = bcg_reactome$clust_distances %>%
                                  map(summary) %>%
                                  map(~.x$mean) %>%
                                  reduce(c) %>%
                                  range,
                                midpoint = bcg_reactome$clust_distances %>%
                                  map(summary) %>%
                                  map(~.x$mean) %>%
                                  reduce(c) %>%
                                  range %>%
                                  mean,
                                name = 'mean cosine\ndistance'))

# Result tables -------

  insert_msg('Result tables')

  ## for the common regulated signatures

  bcg_reactome$result_tbl <-
    map2(bcg_reactome$stats,
         map(bcg_reactome$test,
             ~.x[c('variable', 'significance', 'eff_size')]),
         left_join, by = 'variable') %>%
    compress(names_to = 'cohort') %>%
    filter(variable %in% bcg_reactome$common_significant) %>%
    mutate(cohort = globals$cohort_labs[cohort],
           variable = ifelse(!stri_detect(variable, regex = '^Samp'),
                             exchange(variable,
                                      bcg_reactome$lexicon),
                             variable))

  bcg_reactome$result_tbl <- bcg_reactome$result_tbl %>%
    relocate(cohort, variable) %>%
    set_names(c('Cohort',
                'Variable',
                levels(bcg_reactome$data[[1]]$clust_id),
                'Significance',
                'Effect size'))

# END -----

  bcg_reactome$data <- NULL
  bcg_reactome$lexicon <- NULL
  bcg_reactome$cluster_descr <- NULL
  bcg_reactome$cluster_tcga <- NULL

  bcg_reactome <- compact(bcg_reactome)

  plan('sequential')

  insert_tail()
