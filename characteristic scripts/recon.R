# GSVA for the Recon metabolic subsystem gene signatures.
#
# Differences in ssGSEA scores between the hormonal clusters are investigated
# by one-way ANOVA with eta-square effect size metric.
#
# Differentially regulated signatures are defined by pFDR < 0.05
# and eta-square >= 0.14.

  insert_head()

# container ------

  bcg_recon <- list()

# analysis variables and data -------

  insert_msg('Analysis variables and data')

  bcg_recon$data <-
    map2(bcg_globals$assignment[names(recon$scores)],
         recon$scores,
         inner_join, by = 'sample_id')

  bcg_recon$lexicon <- recon$lexicon %>%
    filter(variable != 'Unassigned')

# N numbers -------

  insert_msg('N numbers')

  ## ready-to-use plot captions

  bcg_recon$n_captions <- bcg_recon$data %>%
    map(count, clust_id) %>%
    map(~map2_chr(.x[[1]], .x[[2]],
                  paste, sep = ': n = ')) %>%
    map(paste, collapse = ', ')

# descriptive stats --------

  insert_msg('Descriptive stats')

  bcg_recon$stats <- bcg_recon$data %>%
    map(fast_num_stats,
        split_factor = 'clust_id')

# Testing for differences between the hormonal clusters ------

  insert_msg('Testing')

  bcg_recon$test <- bcg_recon$data %>%
    map(test_anova,
        variables = bcg_recon$lexicon$variable,
        split_fct = 'clust_id',
        .parallel = FALSE) %>%
    map(~.x$anova)

# Formatting the results and significant effects -------

  insert_msg('Formatting the results and significant effects')

  bcg_recon$test <- bcg_recon$test %>%
    ## to get nicely formatted significances
    map(re_adjust, method = 'none') %>%
    map(mutate,
        variable = response,
        eff_size = paste('\u03B7\u00B2 =', signif(effect_size, 2)),
        ## axis labels to be used (optionally) at heat map Y axes
        ## significant effects are highlighted in bold
        axis_lab = exchange(variable,
                            bcg_recon$lexicon),
        axis_lab = ifelse(p_adjusted < 0.05 & effect_size >= 0.14,
                          html_bold(axis_lab), axis_lab))

  ## significant in single cohorts
  ## and common significant effect shared by the TCGA and GSE66420 cohorts

  bcg_recon$significant <- bcg_recon$test %>%
    map(filter,
        p_adjusted < 0.05,
        effect_size >= 0.14) %>%
    map(~.x$variable)

  bcg_recon$common_significant <-
    bcg_recon$significant[c("tcga", "gse99420")] %>%
    reduce(intersect)

# Classification of the common regulated signatures -------

  insert_msg('Classification of the common regulated signatures')

  ## done for the TCGA cohort and applied to the remaining collectives

  bcg_recon$classification_tcga <- bcg_recon$data$tcga %>%
    classify(variables = bcg_recon$common_significant,
             split_fct = 'clust_id') %>%
    .$classification

# Heat maps of ssGSEA scores of the common regulated signatures -------

  insert_msg('Heat maps of ssGSEA scores')

  bcg_recon$hm_plots <-
    list(data = bcg_recon$data,
         plot_title = paste('Recon metabolic signatures,',
                            globals$cohort_labs[names(bcg_recon$data)])) %>%
    pmap(heat_map,
         variables = bcg_recon$common_significant,
         split_fct = 'clust_id',
         normalize = FALSE,
         variable_classification = bcg_recon$classification_tcga,
         limits = c(-1, 1),
         midpoint = 0,
         oob = scales::squish,
         cust_theme = globals$common_theme +
           theme(axis.title.y = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.text.y = element_markdown(size = 7)),
         x_lab = 'cancer sample',
         fill_lab = 'ssGSEA score')

  ## Y axis labels: significant effects highlighted in bold

  bcg_recon$hm_plots <-
    list(x = bcg_recon$hm_plots,
         y = bcg_recon$test) %>%
    pmap(function(x, y) x +
           scale_y_discrete(labels = set_names(y$axis_lab,
                                               y$variable)))

# Distances between the clusters, regulated gene signatures -------

  insert_msg('Distances between the clusters, regulated signatures')

  ## cosine distances convertible to cosine similarities

  bcg_recon$clust_distances <- bcg_recon$data %>%
    map(column_to_rownames, 'sample_id') %>%
    map2(., bcg_recon$significant,
         ~select(.x, clust_id, all_of(.y))) %>%
    map2(., bcg_recon$significant,
         ~subset_distance(.x,
                          variables = .y,
                          split_fct = 'clust_id',
                          dist_FUN = calculate_dist,
                          method = 'cosine'))

  ## summary heat maps of the distances

  bcg_recon$dist_heat_maps <-  bcg_recon$clust_distances %>%
    map(plot,
        'mean',
        cust_theme = globals$common_theme)

  bcg_recon$dist_heat_maps <-
    list(x = bcg_recon$dist_heat_maps,
         y = globals$cohort_labs[names(bcg_recon$dist_heat_maps)],
         z = bcg_recon$n_captions) %>%
    pmap(function(x, y, z) x +
           labs(title = y,
                subtitle = z,
                x = paste(y, 'clusters'),
                y = paste(y, 'clusters')) +
           scale_fill_gradient2(low = 'firebrick',
                                mid = 'white',
                                high = 'steelblue',
                                limits = bcg_recon$clust_distances %>%
                                  map(summary) %>%
                                  map(~.x$mean) %>%
                                  reduce(c) %>%
                                  range,
                                midpoint = bcg_recon$clust_distances %>%
                                  map(summary) %>%
                                  map(~.x$mean) %>%
                                  reduce(c) %>%
                                  range %>%
                                  mean,
                                name = 'mean cosine\ndistance'))

# Result tables ------

  insert_msg('Result tables')

  ## for the common regulated gene signatures

  bcg_recon$result_tbl <-
    map2(bcg_recon$stats,
         map(bcg_recon$test, ~.x[c('variable', 'significance', 'eff_size')]),
         left_join, by = 'variable') %>%
    map(filter, variable %in% bcg_recon$common_significant) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = globals$cohort_labs[cohort],
           variable = ifelse(!stri_detect(variable, regex = '^Samp'),
                             exchange(variable,
                                      bcg_recon$lexicon),
                             variable)) %>%
    relocate(cohort) %>%
    set_names(c('Cohort',
                'Variable',
                levels(bcg_recon$data[[1]]$clust_id),
                'Significance',
                'Effect size'))

# END ------

  bcg_recon$data <- NULL
  bcg_recon$lexicon <- NULL

  bcg_recon <- compact(bcg_recon)

  insert_tail()
