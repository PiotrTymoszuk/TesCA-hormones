# Expression of steroid hormone genes

  insert_head()

# container ------

  bcg_hr <- list()

# analysis globals ------

  insert_msg('Analysis globals')

  ## genes of interest and their classification

  bcg_hr$genes <- c('ESR1', 'ESR2',
                    'ESRRA', 'ESRRB', 'ESRRG',
                    'GPER1',
                    'PGR',
                    'PGRMC1', 'PGRMC2',
                    'AR',
                    'NR3C1', 'NR3C2',
                    'FSHR', 'LHCGR', 'PRLR')

  bcg_hr$classification <-
    c('ESR1' = 'estrogen',
      'ESR2' = 'estrogen',
      'ESRRA' = 'estrogen',
      'ESRRB' = 'estrogen',
      'ESRRG' = 'estrogen',
      'GPER1' = 'estrogen',
      'PGR' = 'progesterone',
      'PGRMC1' = 'progesterone',
      'PGRMC2' = 'progesterone',
      'AR' = 'androgen',
      'NR3C1' = 'corticosteroid',
      'NR3C2' = 'corticosteroid',
      'FSHR' = 'gonadotropin',
      'LHCGR' = 'gonadotropin',
      'PRLR' = 'gonadotropin') %>%
    compress(names_to = 'gene_symbol',
             values_to = 'class')

  ## significant effects in single cohorts

  bcg_hr$significant <- bcg_dge$significant %>%
    map(map, ~.x[.x %in% bcg_hr$genes])

  ## common significant genes

  bcg_hr$common_significant <- bcg_dge$common_significant %>%
    map(~.x[.x %in% bcg_hr$genes])

  # regulation estimates for the genes of interest,
  ## appended with ANOVA results (to be shown in box plots)

  bcg_hr$estimates <- bcg_dge$dev_test %>%
    map(filter,
        gene_symbol %in% bcg_hr$genes) %>%
    map(select,
        clust_id,
        gene_symbol, entrez_id,
        deviation_center, lower_ci, upper_ci,
        regulation)

  bcg_hr$estimates <-
    map2(bcg_hr$estimates,
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

  ## plotting variables: all genes of interest

  bcg_hr$expression <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    map(select, sample_id, any_of(bcg_hr$genes))

  bcg_hr$expression <-
    map2(bcg_globals$assignment[names(bcg_hr$expression)],
         bcg_hr$expression,
         left_join,
         by = 'sample_id')

  ## heat map Y axis scales: significant effects highlighted with
  ## bold font

  bcg_hr$hm_scales <- bcg_hr$estimates %>%
    map(mutate,
        axis_lab = paste(html_italic(gene_symbol),
                         plot_cap,
                         sep = '<br>'),
        axis_lab = ifelse(stri_detect(significance, regex = '^ns'),
                          axis_lab,
                          html_bold(axis_lab))) %>%
    map(~set_names(.x$axis_lab, .x$gene_symbol))

# Classification of the genes by cluster specificity in the TCGA cohort ------

  insert_msg('Classification of the genes in the TCGA cohort')

  bcg_hr$classification_tcga <- bcg_hr$expression$tcga %>%
    classify(variables = bcg_hr$genes,
             split_fct = 'clust_id') %>%
    .$classification

# Heat map plots -------

  insert_msg('Heat maps of Z-scores')

  bcg_hr$hm_plots <-
    list(data = bcg_hr$expression,
         plot_title = paste('Hormone receptor genes,',
                            globals$cohort_labs[names(bcg_hr$expression)])) %>%
    pmap(heat_map,
         variables = bcg_hr$genes,
         split_fct = 'clust_id',
         normalize = TRUE,
         variable_classification = bcg_hr$classification_tcga,
         limits = c(-3, 3),
         oob = scales::squish,
         midpoint = 0,
         x_lab = 'cancer sample',
         cust_theme = globals$common_theme +
           theme(axis.title.y = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.text.y = element_markdown())) %>%
    map2(., bcg_hr$hm_scales,
         ~.x + scale_y_discrete(labels = .y))

# Box plots for single genes --------

  insert_msg('Box plots for single genes')

  bcg_hr$box_subtitles <- bcg_hr$estimates %>%
    map(filter,
        !duplicated(gene_symbol)) %>%
    map(mutate,
        gene_symbol = factor(gene_symbol, bcg_hr$genes)) %>%
    map(arrange, gene_symbol) %>%
    map(~.x$plot_cap)

  ## box plots: ANOVA results are displayed in the plot subtitles

  for(i in names(bcg_hr$expression)) {

    bcg_hr$box_plots[[i]] <-
      list(variable = bcg_hr$genes,
           plot_title = bcg_hr$genes %>%
             html_italic %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = bcg_hr$box_subtitles[[i]]) %>%
      pmap(plot_variable,
           bcg_hr$expression[[i]],
           split_factor = 'clust_id',
           type = 'box',
           cust_theme = globals$common_theme +
             theme(plot.title = element_markdown()),
           x_n_labs = TRUE,
           y_lab = expression('log'[2] * ' expression'),
           x_lab = 'cluster') %>%
      map(~.x +
            scale_fill_manual(values = globals$cluster_colors)) %>%
      set_names(bcg_hr$genes)

  }

# END ------

  rm(i)

  insert_tail()
