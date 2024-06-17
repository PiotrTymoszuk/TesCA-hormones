# Cosine distances between the clusters in respect to differentially
# regulated genes.

  insert_head()

# container ------

  bcg_expdist <- list()

# parallel backend ------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis globals ------

  insert_msg('Analysis globals')

  ## genes found to be differentially regulated between the clusters
  ## specific for each cohort

  bcg_expdist$variables <- bcg_dge$significant %>%
    transpose %>%
    map(reduce, union)

  ## expression data, variables of interest, appending with the cluster
  ## assignment of the samples

  bcg_expdist$data <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression)

  bcg_expdist$data <-
    map2(bcg_expdist$data,
         bcg_expdist$variables[names(bcg_expdist$data)],
         ~select(.x, sample_id, all_of(.y)))

  bcg_expdist$data <-
    map2(bcg_globals$assignment[names(bcg_expdist$data)],
         bcg_expdist$data,
         inner_join, by = 'sample_id') %>%
    map(column_to_rownames, 'sample_id')

# N numbers of samples in the clusters ------

  insert_msg('N numbers')

  ## in form of ready-to-use plot subtitles

  bcg_expdist$n_captions <- bcg_expdist$data %>%
    map(count, clust_id) %>%
    map(~map2_chr(.x[[1]], .x[[2]], paste, sep = ': n = ')) %>%
    map(paste, collapse = ', ')

# Calculation of the cosine distances ------

  insert_msg('Cosine distances')

  bcg_expdist$clust_distances <-
    list(data = bcg_expdist$data,
         variables = bcg_expdist$variables) %>%
    future_pmap(subset_distance,
                split_fct = 'clust_id',
                dist_FUN = calculate_dist,
                method = 'cosine',
                .options = furrr_options(seed = TRUE))

# Summary heat maps -------

  insert_msg('Summary heat maps')

  bcg_expdist$dist_heat_maps <-
    bcg_expdist$clust_distances %>%
    map(plot,
        'mean',
        cust_theme = globals$common_theme)

  ## styling

  bcg_expdist$dist_heat_maps <-
    list(x = bcg_expdist$dist_heat_maps,
         y = globals$cohort_labs[names(bcg_expdist$dist_heat_maps)],
         z = bcg_expdist$n_captions) %>%
    pmap(function(x, y, z) x +
           labs(title = y,
                subtitle = z,
                x = paste(y, 'clusters'),
                y = paste(y, 'clusters')) +
           scale_fill_gradient2(low = 'firebrick',
                                mid = 'white',
                                high = 'steelblue',
                                limits = bcg_expdist$clust_distances %>%
                                  map(summary) %>%
                                  map(~.x$mean) %>%
                                  reduce(c) %>%
                                  range,
                                midpoint = bcg_expdist$clust_distances %>%
                                  map(summary) %>%
                                  map(~.x$mean) %>%
                                  reduce(c) %>%
                                  range %>%
                                  mean,
                                name = 'mean cosine\ndistance'))

# END ------

  bcg_expdist$data <- NULL

  bcg_expdist <- compact(bcg_expdist)

  plan('sequential')

  insert_tail()
