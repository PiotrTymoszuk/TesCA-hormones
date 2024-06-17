# Selection of the clustering algorithm in the TCGA training cohort.
# We're working with the regularized HTK algorithm.
# Lambdas are selected by minimum of BIC, while the k-number of clusters is
# chosen by comparing cluster perofrmance stats in the training data
# and cross-validation

  insert_head()

# container -------

  clust_dev <- list()

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis globals -------

  insert_msg('Analysis globals')

  clust_dev$data <- clust_globals$data$tcga %>%
    center_data

# Creating the clustering objects --------

  insert_msg('Clustering objects')

  clust_dev$algos <-
    list(k = 2:8) %>%
    pmap(htk_cluster,
         data = clust_dev$data,
         lambdas = seq(0, 0.5, by = 0.025),
         select_stat = 'BIC',
         standardize = FALSE) %>%
    set_names(paste0('k', 2:8))

# Performance stats -------

  insert_msg('Performance stats')

  ## in the training data set

  clust_dev$stats$train <- clust_dev$algos %>%
    map(summary) %>%
    compress(names_to = 'k') %>%
    mutate(k = stri_extract(k, regex = '\\d+'),
           k = as.numeric(k),
           accuracy = NA,
           error = NA)

  ## and in cross-validation

  clust_dev$stats$cv <- clust_dev$algos %>%
    future_map(cv,
               nfolds = 5,
               kNN = 9,
               active_variables = TRUE,
               simple_vote = FALSE,
               .parallel = FALSE,
               .options = furrr_options(seed = 1245))

  clust_dev$stats$cv <- clust_dev$stats$cv %>%
    map(summary) %>%
    map(select, ends_with('mean')) %>%
    map(~set_names(.x,
                   stri_replace(names(.x),
                                fixed = '_mean',
                                replacement = ''))) %>%
    compress(names_to = 'k') %>%
    mutate(k = stri_extract(k, regex = '\\d+'),
           k = as.numeric(k))

  ## merging the stats frame

  clust_dev$stats <- clust_dev$stats %>%
    compress(names_to = 'dataset') %>%
    mutate(dataset = factor(dataset, c('train', 'cv')),
           neighborhood_error = 1 - frac_np)

  ## appending with lambdas,
  ## numbers of active variables and fractions of active variables

  clust_dev$lambdas_vars <-
    tibble(k = stri_extract(names(clust_dev$algos), regex = '\\d+'),
           lambda = map_dbl(clust_dev$algos,
                            ~.x$lambdas),
           n_active_variables = map_dbl(clust_dev$algos,
                                        ~length(.x$active_variables))) %>%
    mutate(k = as.numeric(k),
           n_variables = ncol(clust_dev$data),
           frac_active_variables = n_active_variables/n_variables)

  clust_dev$stats <-
    left_join(clust_dev$stats,
              clust_dev$lambdas_vars,
              by = 'k')

# Plots of stats in the training and CV setting -------

  insert_msg('Plots of the performance stats')

  clust_dev$stat_plots <- clust_dev$stats %>%
    pivot_longer(cols = c(sil_width,
                          frac_misclassified,
                          frac_var,
                          neighborhood_error),
                 names_to = 'statistic',
                 values_to = 'value') %>%
    blast(dataset) %>%
    list(x = .,
         y = c('entire data set',
               '5-fold cross-validation')) %>%
    pmap(function(x, y) x %>%
           ggplot(aes(x = k,
                      y = value,
                      color = statistic)) +
           geom_vline(xintercept = 4,
                      linetype = 'dashed') +
           geom_line() +
           geom_point(shape = 16,
                      size = 2) +
           expand_limits(y = 0) +
           #scale_y_continuous(limits = c(0, 1)) +
           scale_color_manual(values = c(frac_misclassified = 'orangered3',
                                         neighborhood_error = 'steelblue',
                                         sil_width = 'darkolivegreen4',
                                         frac_var = 'plum4'),
                              labels = c(frac_misclassified = 'misclassification rate',
                                         neighborhood_error = 'neighborhood error',
                                         sil_width = 'silhouette width',
                                         frac_var = 'explained variance')) +
           globals$common_theme +
           labs(title = 'Cluster number choice',
                subtitle = y,
                x = 'cluster number, k',
                y = 'statistic value'))

# Caching the results ------

  insert_msg('Caching the results')

  clust_dev$data <- NULL
  clust_dev$lambdas_vars <- NULL

  clust_dev <- compact(clust_dev)

  save(clust_dev, file = './cache/clust_dev.RData')

# END ------

  plan('sequential')

  insert_tail()
