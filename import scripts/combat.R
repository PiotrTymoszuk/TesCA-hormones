# ComBat batch adjustment of the expression data

  insert_head()

# container ------

  combat <- list()

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# common genes and expression -------

  insert_msg('Common genes and their expression estimates')

  combat$raw_expression <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    map(filter, tissue == 'tumor') %>%
    map(select, -tissue) %>%
    map(column_to_rownames, 'sample_id')

  combat$variables <- combat$raw_expression %>%
    map(names) %>%
    reduce(intersect)

  combat$raw_expression <- combat$raw_expression %>%
    map(select, all_of(combat$variables))

# Raw expression stats -------

  insert_msg('Raw expression stats')

  combat$stats$raw <- combat$raw_expression %>%
    future_map(distr_stats,
               .options = furrr_options(seed = TRUE)) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = factor(cohort, globals$analysis_cohorts))

# Batch correction ---------

  insert_msg('Batch correction')

  combat$multi_object <-
    multi_process(train = t(combat$raw_expression$tcga),
                  test = combat$raw_expression[c("gse3218", "gse99420")] %>%
                    map(t))

  ## retrieval of the batch-corrected expression data

  combat$expression <- c(list(combat$multi_object$train),
                         combat$multi_object$test) %>%
    set_names(names(combat$raw_expression)) %>%
    map(t) %>%
    map(as.data.frame) %>%
    map(rownames_to_column, 'sample_id') %>%
    map(as_tibble)

# Expression stats upon ComBat adjustment -------

  insert_msg('Adjusted expression stats')

  combat$stats$adjusted <- combat$expression %>%
    map(select, -sample_id) %>%
    future_map(distr_stats,
               .options = furrr_options(seed = TRUE)) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = factor(cohort, globals$analysis_cohorts))

# Plots of the expression stats before and after adjustment ------

  insert_msg('Plots of expression stats before and after adjustment')

  ## ready to use axis labels with N numbers of samples

  combat$n_numbers <- combat$expression %>%
    map_dbl(nrow) %>%
    map2_chr(names(.), .,
             ~paste(globals$cohort_labs[.x], .y, sep = '\nn = ')) %>%
    set_names(names(combat$expression))

  ## box plots of median log2-transformed gene expression values

  combat$plots <- combat$stats %>%
    map(blast, cohort) %>%
    map(map,
        ~tibble(median = median(.x$median),
                perc25 = quantile(.x$median, 0.25),
                perc75 = quantile(.x$median, 0.75),
                perc025 = quantile(.x$median, 0.025),
                perc975 = quantile(.x$median, 0.975))) %>%
    map(compress, names_to = 'cohort') %>%
    map(mutate,
        cohort = factor(cohort, globals$analysis_cohorts))

  combat$plots <- combat$plots %>%
    list(data = .,
         plot_title = c('Raw expression', 'ComBat-adjusted expression'),
         plot_subtitle = combat$stats %>%
           map(~.x$variable) %>%
           map(unique) %>%
           map_dbl(length) %>%
           paste('genes: n =', .)) %>%
    pmap(box_from_stats,
         x_variable = 'cohort',
         fill_variable = 'cohort',
         cust_theme = globals$common_theme,
         x_lab = NULL,
         y_lab = expression('median log'[2] * ' expression')) %>%
    map(~.x +
          scale_fill_manual(values = globals$cohort_colors) +
          scale_x_discrete(labels = combat$n_numbers) +
          guides(fill = 'none'))

# Caching the results -------

  insert_msg('Caching the results')

  combat$raw_expression <- NULL
  combat$variables <- NULL

  combat <- compact(combat)

  save(combat, file = './data/combat.RData')

# END -----

  plan('sequential')

  insert_tail()
