# Distribution of therapy agents in the hormonal clusters - done for TCGA
# cohort patients with available therapy information.
#
# Because there are multiple therapy agents per patient possible, we're
# investigating frequency of particular agents in the clusters of individuals
# with available therapy information.
#
# Statistical significance is determined by Chi-square test with Cramer's V
# effect size statistic.

  insert_head()

# container -----

  bcg_therapy <- list()

# analysis globals ------

  insert_msg('Analysis variables and data')

  ## dummy variables for treatment classes

  bcg_therapy$data <- tcga$treatment %>%
    select(sample_id, treatment_class)

  bcg_therapy$data <- levels(bcg_therapy$data$treatment_class) %>%
    set_names(levels(bcg_therapy$data$treatment_class)) %>%
    map(~filter(bcg_therapy$data, treatment_class %in% .x)) %>%
    map2(., names(.),
         ~mutate(.x, !!.y := 'yes')) %>%
    map(select, -treatment_class) %>%
    reduce(full_join, by = 'sample_id') %>%
    map_dfc(~ifelse(is.na(.x), 'no', .x))

  ## variables

  bcg_therapy$variables <- names(bcg_therapy$data[, -1])

  ## factors and appending with the hormonal cluster assignment

  bcg_therapy$data[, bcg_therapy$variables] <-
    bcg_therapy$data[, bcg_therapy$variables]  %>%
    map_dfc(factor, c('no', 'yes'))

  bcg_therapy$data <-
    left_join(bcg_therapy$data,
              bcg_globals$assignment$tcga,
              by = 'sample_id') %>%
    relocate(sample_id, clust_id)

# N numbers -------

  insert_msg('N numbers')

  bcg_therapy$n_numbers <- bcg_therapy$data %>%
    count(clust_id) %>%
    column_to_rownames('clust_id') %>%
    t %>%
    as_tibble %>%
    mutate(variable = 'Patients, N') %>%
    relocate(variable)

# Descriptive stats -------

  insert_msg('Descriptive stats')

  bcg_therapy$stats <- bcg_therapy$data %>%
    explore(variables = bcg_therapy$variables,
            split_factor = 'clust_id',
            what = 'table',
            pub_styled = TRUE) %>%
    format_stats

# Testing for differences of frequency between the clusters -----

  insert_msg('Testing')

  bcg_therapy$test <- bcg_therapy$data %>%
    compare_variables(variables = bcg_therapy$variables,
                      split_factor = 'clust_id',
                      what = 'eff_size',
                      types = 'cramer_v',
                      exact = FALSE,
                      ci = FALSE,
                      pub_styled = TRUE,
                      adj_method = 'BH') %>%
    mutate(plot_cap = paste(eff_size, significance, sep = ', '))

# Plots -------

  insert_msg('Plots')

  bcg_therapy$plots <-
    list(variable = bcg_therapy$test$variable,
         plot_title = bcg_therapy$test$variable %>%
           stri_capitalize %>%
           paste('treatment,', globals$cohort_labs["tcga"]),
         plot_subtitle = bcg_therapy$test$plot_cap) %>%
    pmap(plot_variable,
         bcg_therapy$data,
         split_factor = 'clust_id',
         scale = 'percent',
         type = 'stack',
         cust_theme = globals$common_theme,
         x_lab = 'cluster',
         y_lab = '% of cluster',
         x_n_labs = TRUE) %>%
    map(~.x +
          scale_fill_manual(values = c(no = 'steelblue',
                                       yes = 'orangered3'),
                            name = '')) %>%
    set_names(bcg_therapy$test$variable)

# Result table ------

  insert_msg('Result table')

  bcg_therapy$result_tbl <-
    left_join(bcg_therapy$stats,
              bcg_therapy$test[c('variable', 'significance', 'eff_size')],
              by = 'variable') %>%
    format_tbl(lexicon = NULL,
               rm_complete = TRUE) %>%
    full_rbind(bcg_therapy$n_numbers, .) %>%
    set_names(c('Variable',
                levels(bcg_therapy$data$clust_id),
                'Significance',
                'Effect size'))

# END -------

  bcg_therapy$data <- NULL
  bcg_therapy$n_numbers <- NULL

  bcg_therapy <- compact(bcg_therapy)

  insert_tail()
