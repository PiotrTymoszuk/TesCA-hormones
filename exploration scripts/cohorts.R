# Characteristic of the study cohorts

  insert_head()

# container -------

  expl_cohorts <- list()

# analysis data ------

  insert_msg('Analysis data')

  expl_cohorts$data <- globals$cohort_expr %>%
    eval %>%
    map(~.x$clinic) %>%
    map(filter,
        tissue == 'tumor')

  expl_cohorts$lexicon <- globals$clinic_lexicon %>%
    filter(class != 'index',
           !variable %in% c('death', 'progression', 'relapse')) %>%
    mutate(test_type = ifelse(format == 'numeric',
                              'kruskal_etasq',
                              'cramer_v'),
           y_lab = ifelse(format == 'factor',
                          '% of cohort', unit))

  ## cohort-specific variable vectors

  expl_cohorts$variables <- expl_cohorts$data %>%
    map(names) %>%
    map(~filter(expl_cohorts$lexicon,
                variable %in% .x)) %>%
    map(~.x$variable)

  expl_cohorts$n_cohorts <- expl_cohorts$variables %>%
    count_features %>%
    transmute(variable = as.character(element),
              n_cohorts = n)

  expl_cohorts$lexicon <-
    left_join(expl_cohorts$lexicon,
              expl_cohorts$n_cohorts,
              by = 'variable') %>%
    mutate(test_type = ifelse(n_cohorts == 1,
                              'none',
                              ifelse(test_type == 'kruskal_eta' & n_cohorts == 2,
                                     'wilcoxon_r', test_type)),
           plot_type = ifelse(test_type == 'none',
                              'none',
                              ifelse(format == 'numeric',
                                     'box', 'stack')))

  ## collapsing the data

  expl_cohorts$data <- expl_cohorts$data %>%
    map2_dfr(., names(.),
             ~mutate(.x,
                     cohort = .y,
                     cohort = globals$cohort_labs[cohort],
                     cohort = factor(cohort, unname(globals$cohort_labs))))

# Descriptive statistics ---------

  insert_msg('Descriptive stats')

  ## only for the TCGA and GSE99420 cohorts

  expl_cohorts$stats <- expl_cohorts$data %>%
    filter(cohort %in% c('TCGA', 'GSE99420')) %>%
    map_dfc(function(x) if(is.factor(x)) droplevels(x) else x) %>%
    explore(variables = reduce(expl_cohorts$variables, union),
            split_factor = 'cohort',
            what = 'table',
            pub_styled = TRUE) %>%
    map(mutate,
        statistic = ifelse(stri_detect(statistic, fixed = 'NaN'),
                           '', statistic)) %>%
    format_stats

# Testing for differences between the cohorts -------

  insert_msg('Testing for differences between the cohorts')

  expl_cohorts$test <- expl_cohorts$lexicon %>%
    filter(variable %in% c('histology', 'relapse_factor'))

  expl_cohorts$test <-
    map2_dfr(expl_cohorts$test$variable,
             expl_cohorts$test$test_type,
             ~compare_variables(expl_cohorts$data %>%
                                  select(cohort, all_of(.x)) %>%
                                  filter(cohort %in% c('TCGA', 'GSE99420'),
                                         complete.cases(.)) %>%
                                  map_dfc(function(x) if(is.factor(x)) droplevels(x) else x),
                                variables = .x,
                                split_factor = 'cohort',
                                what = 'eff_size',
                                types = .y,
                                ci = FALSE,
                                exact = FALSE,
                                pub_styled = TRUE)) %>%
    re_adjust(method = 'BH') %>%
    mutate(plot_cap = paste(eff_size, significance, sep = ', '))

# Plots --------

  insert_msg('Plots')

  expl_cohorts$plots <- expl_cohorts$lexicon %>%
    filter(variable %in% expl_cohorts$test$variable)

  expl_cohorts$plots <-
    list(variable = expl_cohorts$plots$variable,
         plot_title = expl_cohorts$plots$label,
         plot_subtitle = expl_cohorts$test$plot_cap,
         y_lab = expl_cohorts$plots$y_lab,
         type = expl_cohorts$plots$plot_type) %>%
    pmap(plot_variable,
         expl_cohorts$data,
         split_factor = 'cohort',
         scale = 'percent',
         x_n_labs = TRUE,
         cust_theme = globals$common_theme) %>%
    map(~.x +
          scale_fill_brewer(palette = 'Blues') +
          theme(axis.title.x = element_blank())) %>%
    set_names(expl_cohorts$plots$variable)

# Result table --------

  insert_msg('Result table')

  expl_cohorts$result_tbl <-
    left_join(expl_cohorts$stats,
              expl_cohorts$test[c('variable', 'significance', 'eff_size')],
              by = 'variable') %>%
    format_tbl(lexicon = expl_cohorts$lexicon) %>%
    set_names(c('Variable',
                c('TCGA', 'GSE99420'),
                'Significance',
                'Effect size'))

# END -------

  expl_cohorts$data <- NULL
  expl_cohorts$lexicon <- NULL
  expl_cohorts$variables <- NULL
  expl_cohorts$n_cohorts <- NULL

  expl_cohorts <- compact(expl_cohorts)

  insert_tail()
