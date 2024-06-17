# Expression of steroid hormone genes

  insert_head()

# container ------

  bcg_gfr <- list()

# analysis globals ------

  insert_msg('Analysis globals')

  ## genes of interest and their classification

  bcg_gfr$genes <- c('EGFR',
                    paste0('ERBB', 2:4),

                    c('EGF', 'TGFA', 'HBEGF',
                      'AREG', 'BTC', 'EPGN',
                      'EREG', paste0('NRG', 1:4)),

                    paste0('FGFR', 1:4),
                    paste0('FGF', 1:23))

  ## significant effects in single cohorts

  bcg_gfr$significant <- bcg_dge$significant %>%
    map(map, ~.x[.x %in% bcg_gfr$genes])

  ## common significant genes

  bcg_gfr$common_significant <- bcg_dge$common_significant %>%
    map(~.x[.x %in% bcg_gfr$genes])

  ## regulation estimates for the genes of interest,
  ## appended with ANOVA results (to be shown in box plots)

  bcg_gfr$estimates <- bcg_dge$dev_test %>%
    map(filter,
        gene_symbol %in% bcg_gfr$genes) %>%
    map(select,
        clust_id,
        gene_symbol, entrez_id,
        deviation_center, lower_ci, upper_ci,
        regulation)

  bcg_gfr$estimates <-
    map2(bcg_gfr$estimates,
         map(bcg_dge$anova,
             ~.x[c('gene_symbol', 'effect_size', 'p_adjusted')]),
         left_join,
         by = 'gene_symbol') %>%
    map(re_adjust, 'p_adjusted', method = 'none') %>%
    map(mutate,
        plot_cap = paste('\u03B7\u00B2 =', signif(effect_size, 2)),
        plot_cap = paste(plot_cap, significance, sep = ', '))

# Plotting variables and expression values for heat maps ------

  insert_msg('Plotting variables and expression values')

  ## plotting variables: genes found to be differentially
  ## regulated in at least two cohorts

  bcg_gfr$plot_variables <- bcg_gfr$common_significant %>%
    unlist %>%
    unname %>%
    unique

  bcg_gfr$expression <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    map(select, sample_id, any_of(bcg_gfr$plot_variables))

  bcg_gfr$expression <-
    map2(bcg_globals$assignment[names(bcg_gfr$expression)],
         bcg_gfr$expression,
         left_join,
         by = 'sample_id')

# Classification of the genes by cluster specificity in the TCGA cohort ------

  insert_msg('Classification of the genes in the TCGA cohort')

  bcg_gfr$classification_tcga <- bcg_gfr$expression$tcga %>%
    classify(variables = bcg_gfr$plot_variables,
             split_fct = 'clust_id') %>%
    .$classification

# Heat map plots -------

  insert_msg('Heat maps of Z-scores')

  bcg_gfr$hm_plots <-
    list(data = bcg_gfr$expression,
         plot_title = paste('Hormone receptor genes,',
                            globals$cohort_labs[names(bcg_gfr$expression)])) %>%
    pmap(heat_map,
         variables = bcg_gfr$plot_variables,
         split_fct = 'clust_id',
         normalize = TRUE,
         variable_classification = bcg_gfr$classification_tcga,
         limits = c(-3, 3),
         oob = scales::squish,
         midpoint = 0,
         x_lab = 'cancer sample',
         cust_theme = globals$common_theme +
           theme(axis.title.y = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.text.y = element_markdown()))

# Box plots for single genes --------

  insert_msg('Box plots for single genes')

  bcg_gfr$box_subtitles <- bcg_gfr$estimates %>%
    map(filter,
        gene_symbol %in% bcg_gfr$plot_variables,
        !duplicated(gene_symbol)) %>%
    map(mutate,
        gene_symbol = factor(gene_symbol, bcg_gfr$plot_variables)) %>%
    map(arrange, gene_symbol) %>%
    map(~.x$plot_cap)

  ## box plots: ANOVA results are displayed in the plot subtitles

  for(i in names(bcg_gfr$expression)) {

    bcg_gfr$box_plots[[i]] <-
      list(variable = bcg_gfr$plot_variables,
           plot_title = bcg_gfr$plot_variables %>%
             html_italic %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = bcg_gfr$box_subtitles[[i]]) %>%
      pmap(plot_variable,
           bcg_gfr$expression[[i]],
           split_factor = 'clust_id',
           type = 'box',
           cust_theme = globals$common_theme +
             theme(plot.title = element_markdown()),
           x_n_labs = TRUE,
           y_lab = expression('log'[2] * ' expression'),
           x_lab = 'cluster') %>%
      map(~.x +
            scale_fill_manual(values = globals$cluster_colors)) %>%
      set_names(bcg_gfr$plot_variables)

  }

# END ------

  bcg_gfr <- bcg_gfr[c("genes",
                       "significant", "common_significant",
                       "plot_variables",
                       "hm_plots", "box_plots")]

  rm(i)

  insert_tail()
