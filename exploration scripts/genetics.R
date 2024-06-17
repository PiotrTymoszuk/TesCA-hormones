# General frequency of somatic mutations and CNA in the TCGA cohort.

  insert_head()

# container ------

  expl_genet <- list()

# parallel backend -----

  insert_msg('Parallel backend')

  plan('multisession')

# analysis data and variables ------

  insert_msg('Analysis variables and data')

  expl_genet$data[c("mutations",
                    "deletions",
                    "amplifications")] <- tcga[c("mutations",
                                                 "deletions",
                                                 "amplifications")] %>%
    map(filter, tissue == 'tumor') %>%
    map(select, -tissue)

  expl_genet$variables <- expl_genet$data %>%
    map(select, -sample_id) %>%
    map(names)

# Frequency ------

  insert_msg('Descriptive stats')

  expl_genet$stats <- expl_genet$data %>%
    map(column_to_rownames, 'sample_id') %>%
    future_map(count_binary,
               .options = furrr_options(seed = TRUE))

# Top variables: features present in at least 2.5% of samples -------

  insert_msg('Top most frequent genetic features')

  ## they will be used later for characteristic of the hormonal clusters

  expl_genet$top_features <- expl_genet$stats %>%
    map(filter, percent >= 2.5) %>%
    map(~.x$variable)

# END ------

  expl_genet$data <- NULL
  expl_genet$variables <- NULL

  expl_genet <- compact(expl_genet)

  plan('sequential')

  insert_tail()
