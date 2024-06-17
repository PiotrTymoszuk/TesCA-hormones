# Response to therapy between the clusters: generally, this analysis is
# definitively under-powered and won't be shown in the manuscript -
# only on request!

  insert_head()

# container ------

  bcg_resp <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  ## best overall response

  bcg_resp$data$`best overall` <- tcga$treatment %>%
    filter(!is.na(response)) %>%
    select(sample_id, response) %>%
    blast(sample_id) %>%
    map(arrange, desc(response)) %>%
    map_dfr(~.x[1, ])

  ## data frames for particular therapy agents
  ## the 'other' category is removed

  bcg_resp$data <- tcga$treatment %>%
    filter(treatment_class != 'other',
           !is.na(response)) %>%
    blast(treatment_class, .skip = TRUE) %>%
    map(select, sample_id, response) %>%
    c(bcg_resp$data, .)

  bcg_resp$data <- bcg_resp$data %>%
    map(mutate,
        resposen = droplevels(response))

  ## appending with the cluster assignment

  bcg_resp$data <- bcg_resp$data %>%
    map(left_join,
        bcg_globals$assignment$tcga,
        by = 'sample_id') %>%
    map(relocate,
        sample_id, clust_id)

# Descriptive stats --------

  insert_msg('Descriptive stats')

  bcg_resp$stats <- bcg_resp$data %>%
    map(explore,
        variables = 'response',
        split_factor = 'clust_id',
        what = 'table',
        pub_styled = TRUE) %>%
    map_dfr(format_stats) %>%
    mutate(treatment_type = names(bcg_resp$data))

# Tests -------

  insert_msg('Tests')

  bcg_resp$test <- bcg_resp$data %>%
    map(compare_variables,
        variables = 'response',
        split_factor = 'clust_id',
        what = 'eff_size',
        types = 'cramer_v',
        exact = FALSE,
        ci = FALSE,
        pub_styled = TRUE) %>%
    compress(names_to = 'treatment_type') %>%
    re_adjust(method = 'BH') %>%
    mutate(plot_cap = paste(eff_size, significance, sep = ', '))

# Plots -------

  insert_msg('Plots')

  bcg_resp$plots <-
    list(x = bcg_resp$data,
         y = bcg_resp$test$treatment_type %>%
           stri_capitalize %>%
           paste('response,',
                 globals$cohort_labs["tcga"]),
         z = bcg_resp$test$plot_cap) %>%
    pmap(function(x, y, z) x %>%
           plot_variable(variable = 'response',
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
          scale_fill_manual(values = c('steelblue',
                                       'gray70',
                                       'orangered3'),
                            name = '')) %>%
    set_names(bcg_resp$test$treatment_type)

# Result table ------

  insert_msg('Result table')

  bcg_resp$result_tbl <-
    left_join(bcg_resp$stats %>%
                select(-variable),
              bcg_resp$test[c("treatment_type", "significance", "eff_size")],
              by = 'treatment_type') %>%
    format_tbl(lexicon = NULL) %>%
    relocate(treatment_type) %>%
    set_names(c('Treatment_type',
                levels(bcg_resp$data$`best overall`$clust_id),
                'Significance',
                'Effect size'))

# END ------

  bcg_resp$data <- NULL

  bcg_resp <- compact(bcg_resp)

  insert_tail()
