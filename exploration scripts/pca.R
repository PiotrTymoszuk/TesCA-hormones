# Principal component analysis and clustering tendency
# measured by Hopking statistic in the hormone related gene expression data.
#
# Z-scores of ComBat-adjusted log2-transformed expression of the genes passing
# variability and minimal expressin criteria are included in the analysis.

  insert_head()

# container -------

  expl_pca <- list()

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis variables and data ----------

  insert_msg('Analysis variables and data')

  expl_pca$variables <- expl_dist$top_variables

  expl_pca$data <- combat$expression %>%
    map(select,
        sample_id, all_of(expl_pca$variables)) %>%
    map(column_to_rownames, 'sample_id') %>%
    map(center_data)

# PCA objects -------

  insert_msg('PCA objects')

  ## six dimensions: as evident from the scree plots, first 3 - 4 do
  ## actually play the role

  expl_pca$pca_objects <- expl_pca$data %>%
    map(reduce_data,
        kdim = 6,
        red_fun = 'pca')

# Scree, score and loadings plots -----

  insert_msg('PCA plots')

  ## scree plots

  expl_pca$scree_plots <-
    list(x = expl_pca$pca_objects,
         segment_color = globals$cohort_colors[names(expl_pca$pca_objects)]) %>%
    pmap(plot,
         type = 'scree',
         cust_theme = globals$common_theme)

  ## score plots

  expl_pca$score_plots <-
    list(x = expl_pca$pca_objects,
         point_color = globals$cohort_colors[names(expl_pca$pca_objects)]) %>%
    pmap(plot,
         type = 'scores',
         cust_theme = globals$common_theme) %>%
    map(~.x +
          geom_vline(xintercept = 0,
                     linetype = 'dashed') +
          geom_hline(yintercept = 0,
                     linetype = 'dashed'))

  ## loadings plots, first two dimensions

  expl_pca$loadings_plots <-
    list(x = expl_pca$pca_objects,
         point_color = globals$cohort_colors[names(expl_pca$pca_objects)]) %>%
    pmap(plot,
         type = 'loadings',
         cust_theme = globals$common_theme,
         txt_type = 'text') %>%
    map(~.x +
          geom_vline(xintercept = 0,
                     linetype = 'dashed') +
          geom_hline(yintercept = 0,
                     linetype = 'dashed'))

  ## additional styling

  for(i in c('scree_plots', 'score_plots', 'loadings_plots')) {

    expl_pca[[i]] <-
      list(x = expl_pca[[i]],
           y = globals$cohort_labs[names(expl_pca[[i]])],
           z = map_dbl(expl_pca$data, nrow)) %>%
      pmap(function(x, y, z) x +
             labs(title = y,
                  subtitle = paste('samples: n =', z)) +
             theme(plot.tag = element_blank()))

  }

# Key genes for the principal components --------

  insert_msg('Key genes for the principal components')

  ## loadings for the principal components, appending with the gene class

  expl_pca$loadings <- expl_pca$pca_objects %>%
    map(extract, 'loadings') %>%
    map(left_join,
        set_names(globals$gene_lexicon[c('gene_symbol', 'class')],
                  c('variable', 'class')),
        by = 'variable') %>%
    map(relocate, variable, class)

  ## variances associated with the PCs

  expl_pca$variance <- expl_pca$pca_objects %>%
    map(clustTools::var)

  ## bar plots with the loadings

  for(i in names(expl_pca$loadings)) {

    expl_pca$pc_plots[[i]] <-
      list(x = names(expl_pca$loadings[[i]])[-2:-1],
           y = paste0('PC',
                      1:(ncol(expl_pca$loadings[[i]]) - 2)),
           z = expl_pca$variance[[i]]$perc_var) %>%
      pmap(function(x, y, z) expl_pca$loadings[[i]] %>%
             ggplot(aes(x = .data[[x]],
                        y = reorder(variable, .data[[x]]),
                        fill = class)) +
             geom_bar(stat = 'identity',
                      color = 'black') +
             scale_fill_manual(values = set_names(globals$gene_lexicon$color,
                                                  globals$gene_lexicon$class),
                               name = 'Gene classification') +
             globals$common_theme +
             theme(axis.title.y = element_blank(),
                   axis.text.y = element_text(face = 'italic')) +
             labs(title = paste(y, globals$cohort_labs[i], sep = ', '),
                  subtitle = paste0('fraction of variance: ', signif(z, 2), '%'),
                  x = 'Loading value')) %>%
      set_names(names(expl_pca$loadings[[i]])[-2:-1])

  }

# UMAP --------

  insert_msg('Two-dimensional UMAP')

  expl_pca$umap_objects <- expl_pca$data %>%
    map(reduce_data,
        kdim = 2,
        distance_method = 'cosine',
        red_fun = 'umap',
        random_state = 12345)

  ## UMAP embedding plots

  expl_pca$umap_plots <-
    list(x = expl_pca$umap_objects,
         point_color = globals$cohort_colors[names(expl_pca$umap_objects)]) %>%
    pmap(plot,
         cust_theme = globals$common_theme)

  expl_pca$umap_plots <-
    list(x = expl_pca$umap_plots,
         y = globals$cohort_labs[names(expl_pca$umap_plots)],
         z = map_dbl(expl_pca$data, nrow)) %>%
    pmap(function(x, y, z) x +
           labs(title = y,
                subtitle = paste('samples: n =', z)) +
           theme(plot.tag = element_blank()))

# Spontaneous clustering tendency -------

  insert_msg('Clustering tendency')

  set.seed(12345)

  expl_pca$clust_tendency <- expl_pca$data %>%
    future_map(~get_clust_tendency(.x,
                                   n = floor(0.2 * nrow(.x))),
               .options = furrr_options(seed = TRUE))

# END ------

  expl_pca$data <- NULL
  expl_pca$variables <- NULL

  expl_pca <- compact(expl_pca)

  plan('sequential')

  rm(i)

  insert_tail()
