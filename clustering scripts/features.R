# Differences in expression levels of the cluster-defining genes between
# the clusters.
#
# ComBat-adjusted log2-transformed expression values are compared by
# Kruskal-Wallis tests corrected for multiple testing with the false-discovery
# rate method. Effect size stat: eta-squared

  insert_head()

# container --------

  clust_ft <- list()

# parallel backend ---------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis globals -------

  insert_msg('Analysis globals')

  clust_ft$variables <- clust_globals$variables

  clust_ft$data <- clust_globals$data %>%
    map(rownames_to_column, 'sample_id') %>%
    map2(clust_pred$assignment, .,
         inner_join, by = 'sample_id')

  ## appending with the sample histology,
  ## consecutive observation numbering, which will be used in heat maps
  ## coding for the sample's histology

  clust_ft$clinic <- globals$cohort_expr %>%
    eval %>%
    map(~.x$clinic[c('sample_id', 'histology')]) %>%
    map(mutate,
        histology = ifelse(is.na(histology),
                           'not assigned', as.character(histology)),
        histology = factor(histology, c('seminoma', 'NSGCT', 'not assigned')))

  clust_ft$data <-
    map2(clust_ft$data,
         clust_ft$clinic,
         left_join, by = 'sample_id') %>%
    map(~mutate(.x, observation = paste0('obs_', 1:nrow(.x)))) %>%
    map(relocate,
        sample_id, observation, clust_id, histology)

# N numbers --------

  insert_msg('N numbers')

  clust_ft$n_numbers <- clust_ft$data %>%
    map(count, clust_id) %>%
    map(column_to_rownames, 'clust_id') %>%
    map(t) %>%
    map(as_tibble) %>%
    map(mutate, variable = 'Samples, N') %>%
    map(relocate, variable)

  ## ready-to-use legend labels,
  ## to be used in legends of box plots

  clust_ft$legend_labs <- clust_ft$n_numbers %>%
    map(select, -variable) %>%
    map(~map2_chr(names(.x), unlist(.x[1, ]),
                  paste, sep = ': n = '))

# Descriptive stats --------

  insert_msg('Descriptive stats')

  clust_ft$stats <- clust_ft$data %>%
    map(select, clust_id, all_of(clust_ft$variables)) %>%
    map(fast_num_stats,
        split_factor = 'clust_id')

# Tests ------

  insert_msg('Tests')

  clust_ft$test <- clust_ft$data %>%
    future_map(compare_variables,
               variables = clust_ft$variables,
               split_factor = 'clust_id',
               what = 'eff_size',
               types = 'kruskal_etasq',
               exact = FALSE,
               ci = FALSE,
               pub_styled = TRUE,
               adj_method = 'BH',
               .options = furrr_options(seed = TRUE)) %>%
    map(mutate,
        plot_cap = paste(eff_size, significance),
        gene_lab = html_italic(variable),
        gene_lab = ifelse(p_adjusted < 0.05,
                          html_bold(gene_lab), gene_lab))

  ## significant effects,
  ## common significant effects shared by the TCGA and GSE99420 cohorts

  clust_ft$significant <- clust_ft$test %>%
    map(filter, p_adjusted < 0.05) %>%
    map(~.x$variable)

  clust_ft$common_significant <-
    clust_ft$significant[c("tcga", "gse99420")] %>%
    reduce(intersect)

# Result table --------

  insert_msg('Result table')

  clust_ft$result_tbl <-
    map2(clust_ft$stats,
         map(clust_ft$test, ~.x[c('variable', 'significance', 'eff_size')]),
         left_join, by = 'variable') %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = globals$cohort_labs[cohort]) %>%
    relocate(cohort) %>%
    set_names(c('Cohort',
                'Variable',
                levels(clust_ft$data[[1]]$clust_id),
                'Significance',
                'Effect size'))

# Box plots of mRNA levels for single genes --------

  insert_msg('Box plots')

  for(i in names(clust_ft$data)) {

    clust_ft$plots[[i]] <-
      list(variable = clust_ft$test[[i]]$variable,
           plot_title = clust_ft$test[[i]]$variable %>%
             html_italic %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = clust_ft$test[[i]]$plot_cap) %>%
      pmap(plot_variable,
           clust_ft$data[[i]],
           split_factor = 'clust_id',
           type = 'box',
           cust_theme = globals$common_theme +
             theme(plot.title = element_markdown()),
           y_lab = expression('log'[2] * ' expression'),
           x_n_labs = FALSE) %>%
      map(~.x +
            theme(plot.tag = element_blank()) +
            scale_fill_manual(values = globals$cluster_colors,
                              labels = clust_ft$legend_labs[[i]])) %>%
      set_names(clust_ft$test[[i]]$variable)

  }

# Classification of the variables by specificity for the clusters -------

  insert_msg('Classification of the variables for the clusters')

  ## done in the TCGA training cohort

  clust_ft$classification_tcga <- clust_ft$data$tcga %>%
    classify(variables = clust_ft$variables,
             split_fct = 'clust_id') %>%
    .$classification %>%
    mutate(gene_symbol = variable) %>%
    left_join(globals$gene_lexicon[c('gene_symbol', 'class')],
              by = 'gene_symbol')

# Heat map plots of normalized expression levels -------

  insert_msg('Heat maps of expression Z-scores')

  clust_ft$hm_plots <-
    list(data = clust_ft$data,
         plot_title = globals$cohort_labs[names(clust_ft$data)]) %>%
    pmap(heat_map,
         variables = clust_ft$variable,
         split_fct = 'clust_id',
         normalize = TRUE,
         variable_classification = clust_ft$classification_tcga,
         cust_theme = globals$common_theme +
           theme(axis.title.y = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.text.y = element_markdown()),
         limits = c(-3, 3),
         midpoint = 0,
         oob = scales::squish,
         x_lab = 'cancer sample')

  ## styling: highlighting significant genes

  clust_ft$hm_plots <-
    list(x = clust_ft$hm_plots,
         y = clust_ft$test) %>%
    pmap(function(x, y) x +
           scale_y_discrete(labels = set_names(y$gene_lab, y$variable)))

# Rug heat maps, vertical: gene classification --------

  insert_msg('Rug heat map, gene classification')

  clust_ft$hm_gene <- clust_ft$classification_tcga %>%
    ggplot(aes(x = 'gene\nclass',
               y = reorder(variable, delta_auc),
               fill = class)) +
    facet_grid(clust_id ~ .,
               scales = 'free',
               space = 'free') +
    geom_tile(color = 'black') +
    scale_fill_manual(values = globals$gene_class_colors,
                      name = 'Gene\nclassification') +
    globals$common_theme +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_blank())

# Rug heat maps, horizontal: sample histology -------

  insert_msg('Rug heat maps, sample histology')

  ## plotting data

  clust_ft$hm_sample_data <- clust_ft$hm_plots %>%
    map(~.x$data[c('observation', 'clust_id')]) %>%
    map2(.,
         map(clust_ft$data, ~.x[c('observation', 'histology')]),
         left_join, by = 'observation')

  clust_ft$hm_sample_data <- clust_ft$hm_sample_data %>%
    map(~mutate(.x,
                observation = factor(observation, unique(.x$observation)))) %>%
    map(filter,
        !duplicated(observation)) %>%
    map(arrange, clust_id, observation)

  ## heat maps

  clust_ft$hm_sample <- clust_ft$hm_sample_data %>%
    map(~ggplot(.x,
                aes(x = observation,
                    y = 'histology',
                    fill = histology)) +
          facet_grid(. ~ clust_id,
                     scales = 'free',
                     space = 'free') +
          geom_tile() +
          scale_fill_manual(values = globals$histo_colors,
                            name = 'Histology') +
          globals$common_theme +
          theme(axis.title.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank(),
                axis.line = element_blank(),
                panel.grid.major = element_blank()) +
          labs(x = 'cancer sample'))

# END ------

  clust_ft$data <- NULL
  clust_ft$clinic <- NULL
  clust_ft$hm_sample_data <- NULL

  clust_ft <- compact(clust_ft)

  rm(i)

  plan('sequential')

  insert_tail()
