# Comparison of numeric features of genetic stability (TMB, MSI scores) in the
# hormonal clusters of the TCGA cohort.
#
# Statistical significance is determined by Kruskal-Wallis test with eta-square
# effect size statistic.

  insert_head()

# container -----

  bcg_burdens <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  ## information provided by the study authors

  bcg_burdens$variables <-
    c("tmb", "mantis_msi_score", "sensor_msi_score")

  bcg_burdens$data$clinic <-
    tcga$clinic[, c("sample_id", bcg_burdens$variables)]

  ## counts of somatic mutations, deletions, and amplifications
  ## with the provided data sets

  bcg_burdens$data$panels <-
    tcga[c("mutations", "deletions", "amplifications")] %>%
    map(select, -tissue) %>%
    map(column_to_rownames, 'sample_id') %>%
    map(rowSums) %>%
    map2(., names(.),
         ~compress(.x,
                   names_to = 'sample_id',
                   values = .y)) %>%
    reduce(left_join, by = 'sample_id')

  ## the entire analysis data frame, appending with the
  ## cluster assignment

  bcg_burdens$data <- bcg_burdens$data %>%
    reduce(left_join, by = 'sample_id') %>%
    inner_join(bcg_globals$assignment$tcga, ., by = 'sample_id')

  ## update of teh variables and a variable lexicon

  bcg_burdens$variables <-
    c(bcg_burdens$variables,
      c('mutations', 'deletions', 'amplifications'))

  bcg_burdens$lexicon <-
    tibble(variable = bcg_burdens$variables,
           label = c('Total mutation burden',
                     'MANTIS MSI score',
                     'SENSOR MSI score',
                     'Mutation number',
                     'Gene deletion number',
                     'Gene amplification number'),
           y_label = c('alterations per MB',
                       'score',
                       'score',
                       'count per sample',
                       'count per sample',
                       'count per sample'),
           table_label = c('Total mutation burden, alterations per MB',
                           'MANTIS MSI score',
                           'SENSOR MSI score',
                           'Mutation number',
                           'Gene deletion number',
                           'Gene amplification number'))

# Descriptive stats --------

  insert_msg('Descriptive stats')

  bcg_burdens$stats <- bcg_burdens$data %>%
    select(-sample_id) %>%
    fast_num_stats(split_factor = 'clust_id')

# Tests -----

  insert_msg('Tests')

  bcg_burdens$test <- bcg_burdens$data %>%
    compare_variables(variables = bcg_burdens$variables,
                      split_factor = 'clust_id',
                      what = 'eff_size',
                      types = 'kruskal_etasq',
                      exact = FALSE,
                      ci = FALSE,
                      pub_styled = TRUE,
                      adj_method = 'BH') %>%
    mutate(plot_cap = paste(eff_size, significance, sep = ', '))

# Box plots --------

  insert_msg('Box plots')

  bcg_burdens$plots <-
    list(variable = bcg_burdens$lexicon$variable,
         plot_title = bcg_burdens$lexicon$label,
         plot_subtitle = bcg_burdens$test$plot_cap,
         y_lab = bcg_burdens$lexicon$y_label) %>%
    pmap(plot_variable,
         bcg_burdens$data,
         split_factor = 'clust_id',
         type = 'box',
         cust_theme = globals$common_theme,
         x_lab = 'cluster',
         x_n_labs = TRUE,
         point_hjitter = 0) %>%
    map(~.x +
          scale_fill_manual(values = globals$cluster_colors)) %>%
    set_names(bcg_burdens$lexicon$variable)

# Result table -------

  insert_msg('Result table')

  bcg_burdens$result_tbl <-
    left_join(bcg_burdens$stats,
              bcg_burdens$test[c("variable", "significance", "eff_size")],
              by = 'variable') %>%
    mutate(variable = ifelse(!stri_detect(variable, regex = '^Samp'),
                             exchange(variable,
                                      bcg_burdens$lexicon,
                                      value = 'table_label'),
                             variable)) %>%
    set_names(c('Variable',
                levels(bcg_burdens$data$clust_id),
                'Significance',
                'Effect size'))

# END -------

  bcg_burdens$variables <- NULL
  bcg_burdens$lexicon <- NULL
  bcg_burdens$data <- NULL

  bcg_burdens <- compact(bcg_burdens)

  insert_tail()
