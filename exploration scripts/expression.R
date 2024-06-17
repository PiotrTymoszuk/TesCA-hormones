# Comparison of expression of the hormone-related genes between the cohort.
# Kruskal-Wallis test with eta-squared effect size metric is used for comparison
# of ComBat-adjusted log2-corrected mRNA levels between the cohorts.
# The genes passing variabilty and minimal expression criteria are used
# in the analysis.

  insert_head()

# container ------

  expl_expr <- list()

# analysis data ------

  insert_msg('Analysis data')

  expl_expr$variables <- expl_dist$top_variables

  expl_expr$data <- combat$expression %>%
    map(select, sample_id, all_of(expl_expr$variables)) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = globals$cohort_labs[cohort],
           cohort = factor(unname(cohort), unname(globals$cohort_labs))) %>%
    relocate(cohort)

# N numbers -------

  insert_msg('N numbers')

  expl_expr$n_numbers <- expl_expr$data %>%
    count(cohort) %>%
    column_to_rownames('cohort') %>%
    t %>%
    as_tibble %>%
    mutate(variable = 'Samples, N') %>%
    relocate(variable)

# Descriptive stats ------

  insert_msg('Descriptive stats')

  expl_expr$stats <- expl_expr$data %>%
    explore(variables = expl_expr$variables,
            split_factor = 'cohort',
            what = 'table',
            pub_styled = TRUE) %>%
    format_stats

# Tests --------

  insert_msg('Tests')

  expl_expr$test <- expl_expr$data %>%
    compare_variables(variables = expl_expr$variables,
                      split_factor = 'cohort',
                      what = 'eff_size',
                      types = 'kruskal_etasq',
                      exact = FALSE,
                      ci = FALSE,
                      adj_method = 'BH',
                      pub_styled = TRUE) %>%
    mutate(plot_cap = paste(eff_size, significance))

  ## significant differences

  expl_expr$significant <- expl_expr$test %>%
    filter(p_adjusted < 0.05) %>%
    .$variable

# Box plots -------

  insert_msg('Box plots')

  expl_expr$plots <-
    list(variable = expl_expr$test$variable,
         plot_title = html_italic(expl_expr$test$variable),
         plot_subtitle = expl_expr$test$plot_cap) %>%
    pmap(plot_variable,
         expl_expr$data,
         split_factor = 'cohort',
         type = 'box',
         cust_theme = globals$common_theme,
         x_n_labs = TRUE,
         y_lab = expression('log'[2] * ' expression')) %>%
    map(~.x +
          theme(plot.title = element_markdown()) +
          scale_fill_manual(values = unname(globals$cohort_colors))) %>%
    set_names(expl_expr$test$variable)

# Result table -------

  insert_msg('Result table')

  expl_expr$result_tbl <-
    left_join(expl_expr$stats,
              expl_expr$test[, c('variable', 'significance', 'eff_size')],
              by = 'variable') %>%
    format_tbl(lexicon = NULL,
               rm_complete = TRUE) %>%
    full_rbind(expl_expr$n_numbers, .) %>%
    set_names(c('Variable',
                levels(expl_expr$data$cohort),
                'Significance',
                'Effect size'))

# END --------

  expl_expr$data <- NULL
  expl_expr$n_numbers <- NULL

  expl_expr <- compact(expl_expr)

  insert_tail()
