# Assessment of quality of the hormonal clusters

  insert_head()

# container -------

  clust_eval <- list()

# analysis globals --------

  insert_msg('Analysis globals')

  clust_eval$clust_objects <- clust_pred$clust_objects

# Cluster sizes -------

  insert_msg('Cluster sizes')

  ## cluster frequencies

  clust_eval$clust_sizes <- clust_eval$clust_objects %>%
    map(ngroups) %>%
    map(arrange, desc(clust_id)) %>%
    map(mutate,
        n_total = sum(n),
        percent = n/sum(n) * 100,
        plot_pos = cumsum(percent) - 0.5 * percent) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = factor(cohort, rev(globals$analysis_cohorts)),
           axis_label = paste(globals$cohort_labs[as.character(cohort)],
                              n_total,
                              sep = '\nn = '))

  ## ready to use legend labels and captions with numbers of
  ## observations in the clusters

  clust_eval$clust_size_labs <- clust_eval$clust_objects %>%
    map(ngroups) %>%
    map(~map2_chr(.x[[1]], .x[[2]], paste, sep = ': n = ')) %>%
    map(set_names, levels(clust_eval$clust_sizes$clust_id))

  ## testing for differences in cluster sizes

  clust_eval$clust_size_test <- clust_eval$clust_objects %>%
    map(extract, 'assignment') %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = factor(cohort, globals$analysis_cohorts)) %>%
    compare_variables(variables = 'clust_id',
                      split_factor = 'cohort',
                      what = 'eff_size',
                      types = 'cramer_v',
                      exact = FALSE,
                      ci = FALSE,
                      pub_styled = TRUE) %>%
    mutate(plot_cap = paste(eff_size, significance, sep = ', '))

  ## stack plots with cluster sizes

  clust_eval$clust_size_plot <- clust_eval$clust_sizes %>%
    ggplot(aes(x = percent,
               y = reorder(axis_label, as.numeric(cohort)),
               fill = clust_id)) +
    geom_bar(stat = 'identity',
             color = 'black',
             position = position_stack()) +
    geom_label(aes(label = signif(percent, 2),
                   x = plot_pos),
               size = 2.5,
               show.legend = FALSE) +
    scale_fill_manual(values = globals$cluster_colors,
                      name = 'Hormonal cluster') +
    globals$common_theme +
    theme(axis.title.y = element_blank()) +
    labs(title = 'Hormonal cluster sizes',
         subtitle = clust_eval$clust_size_test$plot_cap,
         x = '% of cancer samples')

# Numeric stats of cluster structure performance ---------

  insert_msg('Numeric stats')

  clust_eval$stats <- clust_eval$clust_objects %>%
    map(summary) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = factor(cohort, rev(globals$analysis_cohorts)),
           neighborhood_error = 1 - frac_np)

  ## bar plots of the clustering stats

  clust_eval$stat_plots <-
    list(x = c('sil_width',
               'frac_misclassified',
               'frac_var',
               'neighborhood_error'),
         y = c('Cluster separation',
               'Misclassification rate',
               'Explained variance',
               'Neighborhood error'),
         z = c('silhouette width',
               'fraction with negative silhouette width ',
               'between-cluster to total sum of squares',
               'fraction of neighbors in different clusters')) %>%
    pmap(function(x, y, z, w) clust_eval$stats %>%
           ggplot(aes(x = .data[[x]],
                      y = cohort,
                      fill = cohort)) +
           geom_bar(stat = 'identity',
                    color = 'black') +
           geom_text(aes(label = signif(.data[[x]], 2)),
                     size = 2.5,
                     hjust = 1.3,
                     color = 'white') +
           scale_fill_manual(values = globals$cohort_colors)  +
           scale_y_discrete(labels = globals$cohort_labs) +
           guides(fill = 'none') +
           globals$common_theme +
           theme(axis.title.y = element_blank()) +
           labs(title = y,
                x = z)) %>%
    set_names(c('sil_width',
                'frac_misclassified',
                'frac_var',
                'neighborhood_error'))

# Pairwise homologous distances --------

  insert_msg('Pairwise homologous distances')

  clust_eval$homo_distances <- clust_eval$clust_objects %>%
    map(cross_distance)

  ## heat map of pairwise squared euclidean distances

  clust_eval$dist_heat_maps <- clust_eval$clust_objects %>%
    map(plot,
        'heat_map',
        cust_theme = globals$common_theme)

  clust_eval$dist_heat_maps <-
    list(x = clust_eval$dist_heat_maps,
         y = globals$cohort_labs[names(clust_eval$dist_heat_map)],
         z = clust_eval$clust_size_labs) %>%
    pmap(function(x, y, z) x +
           theme(axis.text = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks = element_blank(),
                 plot.tag = element_blank()) +
           labs(title = y,
                subtitle = paste(z, collapse = ', ')) +
           scale_fill_gradient2(low = 'firebrick',
                                mid = 'white',
                                high = 'steelblue',
                                limits = c(0, 200),
                                midpoint = 100,
                                oob = scales::squish,
                                name = 'sum-of-squares\ndistance'))

  ## summary heat maps with mean distances

  clust_eval$homo_heat_maps <- clust_eval$homo_distances %>%
    map(plot,
        'mean',
        cust_theme = globals$common_theme)

  clust_eval$homo_heat_maps <-
    list(x = clust_eval$homo_heat_map,
         y = globals$cohort_labs[names(clust_eval$homo_heat_map)],
         z = clust_eval$clust_size_labs) %>%
    pmap(function(x, y, z) x +
           scale_fill_gradient2(low = 'firebrick',
                                mid = 'white',
                                high = 'steelblue',
                                limits = clust_eval$homo_distances %>%
                                  map(summary) %>%
                                  map(~.x$mean) %>%
                                  reduce(c) %>%
                                  range,
                                midpoint = clust_eval$homo_distances %>%
                                  map(summary) %>%
                                  map(~.x$mean) %>%
                                  reduce(c) %>%
                                  range %>%
                                  mean,
                                name = 'mean sum-of-squares\ndistance') +
           labs(title = y,
                subtitle = paste(z, collapse = ', '),
                x = paste(y, 'clusters'),
                y = paste(y, 'clusters')))

# Cross-distances to the training cohort --------

  insert_msg('Cross distance to the training cohort')

  plan('multisession')

  clust_eval$hetero_distances <-
    clust_eval$clust_objects[c("gse3218", "gse99420")] %>%
    map(cross_distance,
        x = clust_eval$clust_objects$gse99420)

  plan('sequential')

  ## heat maps of the mean cross-distances to the training cohort's clusters

  clust_eval$hetero_heat_maps <- clust_eval$hetero_distances %>%
    map(plot,
        'mean',
        cust_theme = globals$common_theme)

  clust_eval$hetero_heat_maps <-
    list(x = clust_eval$hetero_heat_maps,
         y = globals$cohort_labs[names(clust_eval$hetero_heat_maps)]) %>%
    pmap(function(x, y) x +
           scale_fill_gradient2(low = 'firebrick',
                                mid = 'white',
                                high = 'steelblue',
                                limits = clust_eval$hetero_distances %>%
                                  map(summary) %>%
                                  map(~.x$mean) %>%
                                  reduce(c) %>%
                                  range,
                                midpoint = clust_eval$hetero_distances %>%
                                  map(summary) %>%
                                  map(~.x$mean) %>%
                                  reduce(c) %>%
                                  range %>%
                                  mean,
                                name = 'mean sum-of-squares\ndistance') +
           labs(title = y,
                subtitle = 'distances between training and test cohort clusters',
                x = paste(globals$cohort_labs["tcga"], 'clusters'),
                y = paste(y, 'clusters')))

# END ------

  clust_eval$clust_objects <- NULL

  clust_eval <- compact(clust_eval)

  insert_tail()
