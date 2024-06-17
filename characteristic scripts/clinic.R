# Demographic and clinical characteristic of the clusters.
#
# Statistical significance of differences between the clusters is determined
# by Kruskal-Wallis or chi-square test. Effect sizes: eta-square or Cramer's V.


  insert_head()

# container -----

  bcg_clinic <- list()

# analysis globals ------

  insert_msg('Analysis globals')

  ## analysis data
  ## re-coding the histological subtypes in the TCGA cohort

  bcg_clinic$data <- globals$cohort_expr %>%
    eval %>%
    map(~.x$clinic)

  bcg_clinic$data <- map2(bcg_globals$assignment[names(bcg_clinic$data)],
                          bcg_clinic$data,
                          inner_join,
                          by = 'sample_id') %>%
    map(safely_mutate,
        histology_icd = fct_collapse(histology_icd,
                                     TT = c('benign teratoma',
                                            'malignant teratoma',
                                            'teratocarcinoma'),
                                     TYST = c('yolk sac cancer'),
                                     MGCT = c('germinal mixed histology'),
                                     EMBCA = c('embryonal carcinoma'),
                                     SEM = c('seminoma'),
                                     TCCA = c('choriocarcinoma')),
        histology_icd = fct_relevel(histology_icd,
                                    'SEM',
                                    'MGCT',
                                    'EMBCA',
                                    'TT',
                                    'TYST')) %>%
    map(map_dfc,
        function(x) if(is.factor(x)) droplevels(x) else x)

  ## variable lexicon and variables,
  ## removal of survival-related variables
  ##
  ## vectors of variables present in the cohorts
  ## there are no patients who have received neoadjuvant therapy

  bcg_clinic$lexicon <- globals$clinic_lexicon %>%
    filter(class %in% c('demography', 'pathology', 'treatment'),
           variable != 'neoadjuvant') %>%
    mutate(test_type = ifelse(format == 'numeric',
                              'kruskal_etasq', 'cramer_v'),
           plot_type = ifelse(format == 'numeric',
                              'box', 'stack'),
           axis_lab = ifelse(is.na(unit),
                             '% of cluster', unit))

  bcg_clinic$variables <- bcg_clinic$data %>%
    map(names) %>%
    map(~filter(bcg_clinic$lexicon, variable %in% .x)) %>%
    map(~.x$variable)

  ## a vector with numeric variables, which will be used later
  ## for setting fill scales for the plots

  bcg_clinic$numeric_variables <- bcg_clinic$lexicon %>%
    filter(format == 'numeric') %>%
    .$variable

# Descriptive stats --------

  insert_msg('Descriptive stats')

  bcg_clinic$stats <-
    map2(bcg_clinic$data,
         bcg_clinic$variables,
         ~explore(.x,
                  variables = .y,
                  split_factor = 'clust_id',
                  what = 'table',
                  pub_styled = TRUE)) %>%
    map(format_stats)

# Testing for differences between the cohorts -------

  insert_msg('Testing for differences between the cohorts')

  bcg_clinic$test <-
    map2(bcg_clinic$data,
         bcg_clinic$variables,
         ~compare_variables(.x,
                            variables = .y,
                            split_factor = 'clust_id',
                            what = 'eff_size',
                            types = exchange(.y,
                                             bcg_clinic$lexicon,
                                             value = 'test_type'),
                            ci = FALSE,
                            exact = FALSE,
                            pub_styled = TRUE,
                            adj_method = 'BH')) %>%
    map(mutate, plot_cap = paste(eff_size, significance))

  ## significant differences

  bcg_clinic$significant <- bcg_clinic$test %>%
    map(filter, p_adjusted < 0.05) %>%
    map(~.x$variable)

# Plots ------

  insert_msg('Plots')

  for(i in names(bcg_clinic$data)) {

    bcg_clinic$plots[[i]] <-
      list(variable = bcg_clinic$test[[i]]$variable,
           plot_title = bcg_clinic$test[[i]]$variable %>%
             exchange(bcg_clinic$lexicon) %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = bcg_clinic$test[[i]]$plot_cap,
           y_lab = bcg_clinic$test[[i]]$variable %>%
             exchange(bcg_clinic$lexicon,
                      value = 'axis_lab'),
           type = bcg_clinic$test[[i]]$variable %>%
             exchange(bcg_clinic$lexicon,
                      value = 'plot_type')) %>%
      pmap(plot_variable,
           bcg_clinic$data[[i]],
           split_factor = 'clust_id',
           cust_theme = globals$common_theme,
           scale = 'percent',
           x_n_labs = TRUE,
           x_lab = 'cluster') %>%
      set_names(bcg_clinic$test[[i]]$variable)

    # additional styling

    bcg_clinic$plots[[i]][names(bcg_clinic$plots[[i]]) %in% bcg_clinic$numeric_variables] <-
      bcg_clinic$plots[[i]][names(bcg_clinic$plots[[i]]) %in% bcg_clinic$numeric_variables] %>%
      map(~.x + scale_fill_manual(values = globals$cluster_colors))

    bcg_clinic$plots[[i]][!names(bcg_clinic$plots[[i]]) %in% bcg_clinic$numeric_variables] <-
      bcg_clinic$plots[[i]][!names(bcg_clinic$plots[[i]]) %in% bcg_clinic$numeric_variables] %>%
      map(~.x + scale_fill_brewer(palette = 'Blues'))

    bcg_clinic$plots[[i]]$histology <-  bcg_clinic$plots[[i]]$histology +
      scale_fill_manual(values = globals$histo_colors)


  }

  ## formatting the plots with detailed histology subtyping

  for(i in c('tcga', 'gse3218')) {

    bcg_clinic$plots[[i]]$histology_icd <-
      bcg_clinic$plots[[i]]$histology_icd +
      scale_fill_manual(values = globals$histo_icd_colors)

  }

# Result table -----

  insert_msg('Result table')

  bcg_clinic$result_tbl <-
    map2(bcg_clinic$stats,
         map(bcg_clinic$test, ~.x[c('variable', 'significance', 'eff_size')]),
         left_join, by = 'variable') %>%
    map(format_tbl,
        lexicon = bcg_clinic$lexicon) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = globals$cohort_labs[cohort]) %>%
    relocate(cohort) %>%
    set_names(c('Cohort',
                'Variable',
                levels(bcg_clinic$data[[1]]$clust_id),
                'Significance',
                'Effect size'))

# END -----

  bcg_clinic$data <- NULL
  bcg_clinic$variables <- NULL
  bcg_clinic$lexicon <- NULL

  bcg_clinic <- compact(bcg_clinic)

  rm(i)

  insert_tail()
