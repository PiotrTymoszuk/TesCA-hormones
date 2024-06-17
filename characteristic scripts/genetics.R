# Differences in frequency of somatic mutations, gene amplifications and
# deletions between the hormonal clusters in the TCGA cohort.
# Statistical significance is determined by FDR-corrected chi-square test with
# Cramer's V effect size statistic.
# Features present in at least 2.5% of cancer samples are analyzed.

  insert_head()

# container ------

  bcg_genet <- list()

# analysis data ------

  insert_msg('Analysis data')

  bcg_genet$data[c("mutations", "deletions", "amplifications")] <-
    tcga[c("mutations", "deletions", "amplifications")] %>%
    map(filter, tissue == 'tumor')


  bcg_genet$variables <-
    expl_genet$top_features[names(bcg_genet$data)]

  bcg_genet$data <-
    map2(bcg_genet$data,
         bcg_genet$variables,
         ~.x[c('sample_id', .y)])

  bcg_genet$data <- bcg_genet$data %>%
    map(inner_join,
        bcg_globals$assignment$tcga,
        by = 'sample_id') %>%
    map(relocate,
        sample_id, clust_id)

# Descriptive stats ------

  insert_msg('Frequency of genetic features')

  bcg_genet$stats <- bcg_genet$data %>%
    map(column_to_rownames, 'sample_id') %>%
    map(count_binary,
        split_fct = 'clust_id')

# Testing for differences between the clusters -------

  insert_msg('Testing for differences between the clusters')

  bcg_genet$test <-
    list(x = bcg_genet$data,
         y = bcg_genet$variables,
         z = c(FALSE, TRUE, TRUE)) %>%
    pmap(function(x, y, z) x %>%
           compare_variables(variables = y,
                             split_factor = 'clust_id',
                             what = 'eff_size',
                             types = 'cramer_v',
                             exact = FALSE,
                             ci = FALSE,
                             pub_styled = FALSE,
                             adj_method = 'BH',
                             .parallel = z,
                             .paropts = furrr_options(seed = TRUE,
                                                      globals = c('bcg_genet')))) %>%
    map(mutate,
        eff_size = paste('V =', signif(2, estimate)),
        plot_cap = paste(eff_size, significance, sep = ', '))

  ## significant effects

  bcg_genet$significant <- bcg_genet$test %>%
    map(filter, p_adjusted < 0.05) %>%
    map(~.x$variable)

# Stack plots for all somatic mutations --------

  insert_msg('Stack plots for the mutations')

  bcg_genet$plots$mutations <-
    list(variable = bcg_genet$test$mutations$variable,
         plot_title = bcg_genet$test$mutations$variable %>%
           html_italic %>%
           paste(globals$cohort_labs["tcga"], sep = ' mutations, '),
         plot_subtitle = bcg_genet$test$mutations$plot_cap) %>%
    pmap(plot_variable,
         bcg_genet$data$mutations,
         split_factor = 'clust_id',
         type = 'stack',
         scale = 'percent',
         cust_theme = globals$common_theme +
           theme(plot.title = element_markdown()),
         x_lab = 'cluster',
         y_lab = '% of cluster',
         x_n_labs = TRUE) %>%
    map(~.x +
          scale_fill_manual(values = c('gray80', 'orangered3'),
                            labels = c('WT', 'mutated'),
                            name = '')) %>%
    set_names(bcg_genet$test$mutations$variable)

# Stack plots for selected deletions and gene gains -------

  insert_msg('Stack plots for selected deletions and apmlifications')

  ## plot metainformation

  bcg_genet$amp_del_vars <-
    bcg_genet$test[c("deletions", "amplifications")] %>%
    map(filter, p_value < 0.1) %>%
    map(~.x$variable)

  bcg_genet$ampl_del_colors <-
    list(deletions = c('gray80', 'steelblue'),
         amplifications = c('gray80', 'plum4'))

  bcg_genet$ampl_del_labels <-
    list(deletions = c('WT', 'deleted'),
         amplifications = c('WT', 'amplified'))

  bcg_genet$ampl_del_caps <-
    map2(bcg_genet$test[c("deletions", "amplifications")],
         bcg_genet$amp_del_vars,
         function(x, y) x %>%
           filter(variable %in% y) %>%
           mutate(variable = factor(variable, y)) %>%
           arrange(variable)) %>%
    map(~.x$plot_cap)

  ## plots

  for(i in c("deletions", "amplifications")) {

    bcg_genet$plots[[i]] <-
      list(variable = bcg_genet$amp_del_vars[[i]],
           plot_title = bcg_genet$amp_del_vars[[i]] %>%
             html_italic %>%
             paste(globals$cohort_labs["tcga"],
                   sep = paste0(' ', i, ', ')),
           plot_subtitle = bcg_genet$ampl_del_caps[[i]]) %>%
      pmap(plot_variable,
           bcg_genet$data[[i]],
           split_factor = 'clust_id',
           type = 'stack',
           scale = 'percent',
           cust_theme = globals$common_theme +
             theme(plot.title = element_markdown()),
           x_lab = 'cluster',
           y_lab = '% of cluster',
           x_n_labs = TRUE) %>%
      map(~.x +
            scale_fill_manual(values = bcg_genet$ampl_del_colors[[i]],
                              labels = bcg_genet$ampl_del_labels[[i]],
                              name = '')) %>%
      set_names(bcg_genet$amp_del_vars[[i]])

  }

# Result tables -------

  insert_msg('Result tables')

  bcg_genet$result_tbl <-
    map2(bcg_genet$stats,
         map(bcg_genet$test, ~.x[c('variable', 'significance', 'eff_size')]),
         left_join,
         by = 'variable') %>%
    map(mutate,
        percent = signif(percent, 3)) %>%
    map(select,
        variable,
        clust_id,
        percent,
        n, n_total,
        significance, eff_size) %>%
    map(set_names,
        c('Gene symbol',
          'Cluster',
          'Samples with alterations, percentage',
          'Samples with alterations, N',
          'Total samples',
          'Significance',
          'Effect size')) %>%
    compress(names_to = 'Feature type') %>%
    relocate(`Feature type`)

# Bar plot panels for the most frequent alterations in the clusters ----

  insert_msg('Bar plot panels')

  ## variables of interest: differences between the clusters of at least
  ## V >= 0.2. Significant effects upon FDR correction will be highlighted
  ## by bold colored font, raw significant effects will be highlighted with
  ## by colored font

  bcg_genet$panel_variables <- bcg_genet$test %>%
    map(filter, estimate >= 0.2) %>%
    map(~.x$variable)

  bcg_genet$panel_labs <-
    map2(bcg_genet$test,
         bcg_genet$panel_variables,
         ~filter(.x, variable %in% .y)) %>%
    map2(., c('#cd3700', '#4682B4', '#8b668b'),
         ~mutate(.x,
                 axis_lab = ifelse(p_value < 0.05,
                                   paste0("<i style = 'color:", .y, "'>",
                                          variable, '</i>'),
                                   paste0("<i style = 'color:#666666'>",
                                          variable, '</i>')),
                 axis_lab = ifelse(p_adjusted < 0.05,
                                   html_bold(axis_lab), axis_lab))) %>%
    map(~set_names(.x$axis_lab, .x$variable))

  bcg_genet$panel_stats <-
    map2(bcg_genet$stats,
         bcg_genet$panel_variables,
         ~filter(.x, variable %in% .y))

  ## classification of the alterations by cluster specificity

  bcg_genet$panel_classification <-
    map2(bcg_genet$data,
         bcg_genet$panel_variables,
         ~classify(.x,
                   variables = .y,
                   split_fct = 'clust_id')) %>%
    map(~.x$classification[c('variable', 'clust_id')]) %>%
    map(set_names, c('variable', 'class'))

  bcg_genet$panel_stats <-
    map2(bcg_genet$panel_stats,
         bcg_genet$panel_classification,
         left_join, by = 'variable') %>%
    map(mutate,
        clust_id = factor(clust_id,
                          levels(bcg_globals$assignment[[1]]$clust_id)))

  ## numbers of samples for the column labeller

  bcg_genet$panel_labeller <- bcg_genet$panel_stats %>%
    map(arrange, clust_id) %>%
    map(filter, !duplicated(clust_id)) %>%
    map(~map2_chr(.x[['clust_id']], .x[['n_total']],
                  paste, sep = '\nn = ')) %>%
    map(set_names, levels(bcg_genet$panel_stats[[1]]$clust_id))

  ## bar plots

  bcg_genet$panels <-
    list(x = bcg_genet$panel_stats,
         y = c('#cd3700', '#4682B4', '#8b668b'),
         v = bcg_genet$panel_labeller,
         w = bcg_genet$panel_labs,
         z = paste(c('Top mutations',
                     'Top deletions',
                     'Top amplifications'),
                   globals$cohort_labs["tcga"])) %>%
    pmap(function(x, y, v, w, z) x %>%
           ggplot(aes(x = percent,
                      y = reorder(variable, percent))) +
           geom_bar(stat = 'identity',
                    fill = y,
                    color = 'black') +
           facet_grid(class ~ clust_id,
                      scales = 'free_y',
                      space = 'free_y',
                      labeller = labeller(.cols = v)) +
           scale_y_discrete(labels = w) +
           globals$common_theme +
           theme(axis.text.y = element_markdown(),
                 axis.title.y = element_blank()) +
           labs(title = z,
                x = '% of cluster'))

# END -----

  bcg_genet <-
    bcg_genet[c("variables",
                "stats", "test", "significant",
                "plots", "result_tbl",
                "panel_labeller", "panels")]

  rm(i)

  insert_tail()
