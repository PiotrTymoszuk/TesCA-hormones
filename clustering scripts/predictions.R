# Tuning and training of the cluster classifier (ranger's Random Forest),
# which assigns cancer samples to the hormonal clusters based on expression
# of the hormone-related genes.
# The clusters are established in the TCGA cohort by regularized, unsupervised
# HTK algorithm witk k = 4 clusters.
#
# The rationale behind using Random Forest for the cluster predictions instead
# of kNN classifiers (as usual!): upon normalization/Z-score, expression values
# of the clustering factors are located on completely different scales.
# This may result from different composition of histological types between the
# cohorts. Hence, we're working with a cluster-prediction algorithm that works
# with the genuine ComBat expression levels.

  insert_head()

# container --------

  clust_pred <- list()

# analysis data and cluster assignment in the training cohort -------

  insert_msg('Data and cluster assignment in the training cohort')

  clust_pred$data <- clust_globals$data %>%
    map(rownames_to_column, 'sample_id')

  clust_pred$variables <- clust_globals$variables

  clust_pred$train_assignment <- clust_dev$algos$k4 %>%
    rename(c('1' = '#1',
             '2' = '#2',
             '3' = '#3',
             '4' = '#4')) %>%
    extract('assignment') %>%
    set_names(c('sample_id', 'clust_id'))

  clust_pred$data$tcga <-
    left_join(clust_pred$train_assignment,
              clust_pred$data$tcga,
              by = 'sample_id')

  clust_pred$data <- clust_pred$data %>%
    map(column_to_rownames, 'sample_id')

# tuning of the cluster assignment classifier --------

  insert_msg('Tuning of the cluster assignment classifier')

  set.seed(12345)

  clust_pred$tune_object <-
    tune_rf(clust_pred$data$tcga,
            formula = clust_id ~ .,
            tune_grid = expand.grid(mtry = 2:floor(length(clust_pred$variables)/2),
                                    min.node.size = c(1, 2, 5),
                                    splitrule = c('gini', 'extratrees'),
                                    stringsAsFactors = FALSE),
            num.trees = 1000)

  ## plot of OOB classification errors during of the tuning process

  clust_pred$tune_plot <-
    clust_pred$tune_object$stats %>%
    plot_ranger_tuning(split_by_node_size = TRUE,
                       plot_title = paste('Cluster assignment classifier,',
                                          'tuning,',
                                          globals$cohort_labs["tcga"]))

# Training of the RF classifier and predictions --------

  insert_msg('Training and predictions')

  set.seed(12345)

  clust_pred$model <-
    ranger(formula = clust_id ~ .,
           data = clust_pred$data$tcga,
           num.trees = 1000,
           mtry = clust_pred$tune_object$best_tune$mtry,
           min.node.size = clust_pred$tune_object$best_tune$min.node.size,
           splitrule = clust_pred$tune_object$best_tune$splitrule)

  clust_pred$assignment <- clust_pred$data %>%
    map(~predict(clust_pred$model,
                 data = .x)) %>%
    map(~.x$predictions) %>%
    map2(clust_pred$data, .,
         ~tibble(observation = rownames(.x),
                 clust_id = .y))

# Cluster objects --------

  insert_msg('Cluster objects and cluster assignments')

  ## 'sumofsquares' distance between the observations, i.e.
  ## the distance metric used for development of the genuine clusters
  ## with the HTK algorithm

  clust_pred$clust_objects <-
    list(x = map(clust_globals$data, center_data),
         y = clust_pred$assignment) %>%
    pmap(function(x, y) list(data = quo(!!x),
                             dist_mtx = calculate_dist(x, 'sumofsquares'),
                             dist_method = 'sumofsquares',
                             clust_obj = NULL,
                             clust_fun = 'prediction',
                             clust_assignment = y,
                             dots = list2())) %>%
    map(clust_analysis)

  ## re-naming the cluster assignment data frames:
  ## this will facilitate further use in the analyses

  clust_pred$assignment <- clust_pred$assignment %>%
    map(set_names, c('sample_id', 'clust_id'))

# END ------

  clust_pred$data <- NULL
  clust_pred$variables <- NULL
  clust_pred$train_assignment <- NULL

  clust_pred <- compact(clust_pred)

  insert_tail()
