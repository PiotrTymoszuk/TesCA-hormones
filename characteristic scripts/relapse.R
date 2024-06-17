# Frequency of relapses in the clusters. Statistical significance determined
# by Chi-square tests with Cramer's V effect size statistic.
#
# The analysis is done for the TCGA and GSE99420 cohorts

  insert_head()

# container -------

  bcg_relapse <- list()

# analysis data -------

  insert_msg('Analysis data')

  bcg_relapse$data <- list(tcga = tcga, gse99420 = gse99420) %>%
    map(~.x$clinic[c('sample_id', 'relapse_factor')]) %>%
    map2(bcg_globals$assignment[c("tcga", "gse99420")], .,
         inner_join, by = 'sample_id') %>%
    map(~filter(.x, complete.cases(.x)))

# Descriptive stats ------

  insert_msg('Descriptive stats')

  bcg_relapse$stats <- bcg_relapse$data %>%
    map(explore,
        variables = 'relapse_factor',
        split_factor = 'clust_id',
        what = 'table',
        pub_styled = TRUE) %>%
    map(format_stats)

# Tests -------

  insert_msg('Tests')

  bcg_relapse$test <- bcg_relapse$data %>%
    map(compare_variables,
        variables = 'relapse_factor',
        split_factor = 'clust_id',
        what = 'eff_size',
        types = 'cramer_v',
        exact = FALSE,
        ci = FALSE,
        pub_styled = TRUE) %>%
    map(mutate,
        plot_cap = paste(eff_size, significance, sep = ', '))

# Stack plots -------

  insert_msg('Stack plots')

  bcg_relapse$plots <-
    list(x = bcg_relapse$data,
         y = paste('Relapse,',
                   globals$cohort_labs[names(bcg_relapse$data)]),
         z = map(bcg_relapse$test, ~.x$plot_cap)) %>%
    pmap(function(x, y, z) x %>%
           plot_variable(variable = 'relapse_factor',
                         split_factor = 'clust_id',
                         type = 'stack',
                         scale = 'percent',
                         cust_theme = globals$common_theme,
                         plot_title = y,
                         plot_subtitle = z,
                         x_lab = 'cluster',
                         y_lab = '% of cluster',
                         x_n_labs = TRUE)) %>%
    map(~.x +
          scale_fill_manual(values = c(no = 'cornsilk',
                                       yes = 'orangered3'),
                            name = ''))

# Result table --------

  insert_msg('Result table')

  bcg_relapse$result_tbl <-
    map2(bcg_relapse$stats,
         map(bcg_relapse$test,
             ~.x[c('variable', 'significance', 'eff_size')]),
         left_join, by = 'variable') %>%
    compress(names_to = 'cohort') %>%
    format_tbl(lexicon = NULL) %>%
    mutate(variable = 'Relapse',
           cohort = globals$cohort_labs[cohort]) %>%
    relocate(variable, cohort) %>%
    set_names(c('Cohort', 'Variable',
                levels(bcg_relapse$data[[1]]$clust_id),
                'Significance', 'Effect size'))

# END --------

  bcg_relapse$data <- NULL

  bcg_relapse <- compact(bcg_relapse)

  insert_tail()
