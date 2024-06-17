# Co-expression of the hormone-related genes investigated by pairwise Spearman's
# correlations.
#
# The analysis is done with Combat-corrected log2-transformed expression levels.
# significant correlations: pFDR < 0.05 and abs(rho) >= 0.3.
# The genes passing variability and minimal expression criteria were used in the
# analysis.

  insert_head()

# container ---------

  expl_corr <- list()

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis data ------

  insert_msg('Analysis data')

  expl_corr$variables <- expl_dist$top_variables

  expl_corr$data <- combat$expression %>%
    map(select,
        sample_id, all_of(expl_corr$variables)) %>%
    map(column_to_rownames, 'sample_id') %>%
    map(center_data)

  expl_corr$pairs <- expl_corr$variables %>%
    combn(m = 2, simplify = FALSE)

# Correlation analysis -------

  insert_msg('Correlation analysis')

  for(i in names(expl_corr$data)) {

    expl_corr$test[[i]] <- expl_corr$pairs %>%
      future_map_dfr(~correlate_variables(expl_corr$data[[i]],
                                          variables = .x,
                                          type = 'spearman',
                                          ci = TRUE,
                                          pub_styled = FALSE),
                     .options = furrr_options(seed = TRUE))

  }

  ## formatting the results

  expl_corr$test <- expl_corr$test %>%
    map(re_adjust, method = 'BH') %>%
    map(mutate,
        correlation = ifelse(p_adjusted >= 0.05,
                             'ns',
                             ifelse(estimate > 0.3, 'positive',
                                    ifelse(estimate < -0.3,
                                           'negative', 'ns'))),
        correlation = factor(correlation, c('positive', 'negative', 'ns')),
        pair_id = paste(variable1, variable2, sep = '_'))

# Significant correlations --------

  insert_msg('Significant correlations')

  ## in single cohorts

  expl_corr$significant <- expl_corr$test %>%
    map(filter, correlation %in% c('positive', 'negative')) %>%
    map(blast, correlation) %>%
    transpose %>%
    map(map, ~.x$pair_id)

  ## shared by at least two cohorts

  expl_corr$common_significant <- expl_corr$significant %>%
    map(shared_features, m = 2) %>%
    map(as.character)

# Visualization of the common significant correlations with bubble plots -------

  insert_msg('Bubble plots')

  ## plotting data

  expl_corr$bubble_data <- expl_corr$test %>%
    map(filter,
        pair_id %in% reduce(expl_corr$common_significant, union)) %>%
    map(select,
        variable1, variable2,
        estimate, correlation)

  ## bubble plots

  expl_corr$bubble_plots <-
    list(x = expl_corr$bubble_data,
         y = globals$cohort_labs[names(expl_corr$bubble_data)],
         z = expl_corr$data %>% map_dbl(nrow)) %>%
    pmap(function(x, y, z) x %>%
           ggplot(aes(x = variable1,
                      y = variable2,
                      size = abs(estimate),
                      fill = correlation,
                      color = correlation)) +
           geom_point(shape = 21,
                      color = 'black') +
           #geom_text(aes(label = signif(estimate, 2)),
            #         size = 2.5,
             #        hjust = 0.5,
              #       vjust = -1.2,
               #      show.legend = FALSE) +
           scale_fill_manual(values = c(positive = 'firebrick',
                                        negative = 'steelblue',
                                        ns = 'gray60'),
                             name = 'Correlation') +
           scale_color_manual(values = c(positive = 'firebrick',
                                         negative = 'steelblue',
                                         ns = 'gray60'),
                              name = 'Correlation') +
           scale_size_area(max_size = 6,
                           limits = c(0, 1),
                           name = expression('abs(', rho, ')')) +
           scale_x_discrete(limits = globals$gene_lexicon$gene_symbol) +
           scale_y_discrete(limits = globals$gene_lexicon$gene_symbol) +
           guides(x = guide_axis(angle = 45)) +
           globals$common_theme +
           theme(axis.title = element_blank(),
                 axis.text = element_text(face = 'italic')) +
           labs(title = y,
                subtitle = paste('samples: n =', z)))

# END -------

  expl_corr$data <- NULL
  expl_corr$bubble_data <- NULL
  expl_corr$variables <- NULL
  expl_corr$pairs <- NULL

  expl_corr <- compact(expl_corr)

  plan('sequential')

  insert_tail()
