# Distribution of expression of the hormone-related genes in cancer tissues.
# Identification of outliers by Mahalanobis distance.
#
# ComBat-corrected, log2-transformed mRNA levels are used in the analysis

  insert_head()

# container -------

  expl_dist <- list()

# analysis variables and data -------

  insert_msg('Analysis vaiables and data')

  expl_dist$variables <- globals$genes

  expl_dist$data <- combat$expression %>%
    map(select, sample_id, all_of(expl_dist$variables)) %>%
    map(column_to_rownames, 'sample_id')

# Numbers of observations in the cohorts -------

  insert_msg('Number of observations')

  expl_dist$n_numbers <- expl_dist$data %>%
    map_dbl(nrow) %>%
    map2_chr(names(.), .,
             ~paste(globals$cohort_labs[.x], .y, sep = '\nn = ')) %>%
    set_names(names(expl_dist$data))

# Distribution stats ------

  insert_msg('Distribution statistics')

  ## appended with the gene classification

  expl_dist$stats <- expl_dist$data %>%
    map(distr_stats) %>%
    map(left_join,
        set_names(globals$gene_lexicon[c('gene_symbol', 'class')],
                  c('variable', 'class')),
        by = 'variable')

  ## genes to be used in subsequent analysis: CYP11B1, FSHB, and CYP11B2
  ## are excluded due to low expression and variability

  expl_dist$top_variables <-
    expl_dist$variables[!expl_dist$variables %in% c('CYP11B1', 'CYP11B2', 'FSHB')]

# Distribution plots -------

  insert_msg('Distribution plots')

  ## Gini coefficients versus median
  ## Variance versus mean

  expl_dist$gini_median_plots <-
    list(x = expl_dist$stats,
         y = globals$cohort_labs[names(expl_dist$stats)],
         z = map_dbl(expl_dist$data, nrow)) %>%
    pmap(function(x, y, z) x %>%
           ggplot(aes(x = median,
                      y = gini_coef,
                      fill = class)) +
           geom_point(shape = 21,
                      size = 2,
                      color = 'black') +
           geom_text_repel(aes(label = variable),
                           size = 2.5,
                           fontface = 'italic') +
           scale_fill_manual(values = globals$gene_class_colors,
                             name = 'Gene classification') +
           globals$common_theme +
           labs(title = y,
                subtitle = paste('n =', z),
                x = expression('median log'[2] * ' expression'),
                y = 'Gini coefficient'))

  expl_dist$variance_mean_plots <-
    list(x = expl_dist$stats,
         y = globals$cohort_labs[names(expl_dist$stats)],
         z = map_dbl(expl_dist$data, nrow)) %>%
    pmap(function(x, y, z) x %>%
           ggplot(aes(x = mean,
                      y = var,
                      fill = class)) +
           geom_point(shape = 21,
                      size = 2,
                      color = 'black') +
           geom_text_repel(aes(label = variable),
                           size = 2.5,
                           fontface = 'italic') +
           scale_fill_manual(values = globals$gene_class_colors,
                             name = 'Gene classification') +
           globals$common_theme +
           labs(title = y,
                subtitle = paste('n =', z),
                x = expression('mean log'[2] * ' expression'),
                y = expression('variance log'[2] * ' expression')))

# Normality: Shapiro-Wilk test --------

  insert_msg('Shapiro-Wilk test')

  expl_dist$normality <- expl_dist$data %>%
    map(explore,
        variables = expl_dist$variables,
        what = 'normality',
        pub_styled = TRUE)

# Potential outliers --------

  insert_msg('Outliers')

  ## Mahalanobis distances and their distributions

  expl_dist$mahalanobis <- expl_dist$data %>%
    map(~mahalanobis(.x, center = colMeans(.x), cov = cov(.x))) %>%
    map(compress,
        names_to = 'sample_id',
        values_to = 'distance') %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = factor(cohort, globals$analysis_cohorts))

  ## box plots of the distances

  expl_dist$mahalanobis_plot <-  expl_dist$mahalanobis %>%
    ggplot(aes(x = cohort,
               y = distance,
               fill = cohort)) +
    geom_boxplot(outlier.color = NA,
                 alpha = 0.25) +
    geom_point(size = 2,
               shape = 21,
               alpha = 0.75,
               position = position_jitter(width = 0.1,
                                          height = 0)) +
    scale_fill_manual(values = globals$cohort_colors) +
    scale_x_discrete(labels = expl_dist$n_numbers) +
    guides(fill = 'none') +
    globals$common_theme +
    theme(axis.title.x = element_blank()) +
    labs(title = 'Candidate outliers',
         y = 'Mahalanobis distance from centroid')

# END ------

  expl_dist$data <- NULL
  expl_dist$variables <- NULL

  insert_tail()
