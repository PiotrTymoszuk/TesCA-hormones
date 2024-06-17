# differential expression of estrogen-responsive genes in the hormonal clusters
#
# Sources of the genes:
# Cancer Res (2023) 83 (16): 2656–2674.
# https://doi.org/10.1158/0008-5472.CAN-23-0539
#
# Table S8

  insert_head()

# container -------

  bcg_estro <- list()

# genes of interest ------

  insert_msg('Genes of interest, regulation estimates, and significant features')

  ## genes of interest: some of them are aliases, which need
  ## to be converted to official gene symbols, solving some conflicts per hand

  bcg_estro$gene_lexicon <-
    read_excel('./data/signatures/6780365/can-23-0539_supplementary_table_s8_suppst8.xlsx',
               skip = 1,
               range = 'A2:E89') %>%
    set_names(c('alias', 'mode', 'percent_up', 'percent_down', 'beta_score')) %>%
    mutate(gene_symbol = mapIds(org.Hs.eg.db,
                                keys = alias,
                                keytype = 'ALIAS',
                                column = 'SYMBOL'),
           gene_symbol = ifelse(alias == 'AR',
                                'AR', gene_symbol),
           axis_lab = ifelse(mode == 'Up-regulated',
                             paste0("<i style = 'color:#ff4500'>",
                                    gene_symbol, "</i>"),
                             paste0("<i style = 'color:#8b668b'>",
                                    gene_symbol, "</i>")))

  bcg_estro$genes <- bcg_estro$gene_lexicon$gene_symbol %>%
    unique

  # significant effects in single cohorts

  bcg_estro$significant <- bcg_dge$significant %>%
    map(map, ~.x[.x %in% bcg_estro$genes])

  ## common significant genes

  bcg_estro$common_significant <- bcg_dge$common_significant %>%
    map(~.x[.x %in% bcg_estro$genes])

  # regulation estimates for the genes of interest,
  ## appended with ANOVA results (to be shown in box plots)

  bcg_estro$estimates <- bcg_dge$dev_test %>%
    map(filter,
        gene_symbol %in% bcg_estro$genes) %>%
    map(select,
        clust_id,
        gene_symbol, entrez_id,
        deviation_center, lower_ci, upper_ci,
        regulation)

  bcg_estro$estimates <-
    map2(bcg_estro$estimates,
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

  ## plotting variables: genes found to be regulated in at least one cluster
  ## in both the TCGA and GSE99420 cohort
  ## we're restricting this list to the top differentially regulated genes
  ## defined by effect size in ANOVA of at least 0.26
  ## in the TCGA and GSE99420 cohorts

  bcg_estro$plot_variables <- bcg_estro$estimates[c("tcga", "gse99420")] %>%
    map(filter,
        effect_size >= 0.26,
        gene_symbol %in% unique(unname(unlist(bcg_estro$common_significant)))) %>%
    map(~.x$gene_symbol) %>%
    reduce(intersect)

  bcg_estro$expression <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    map(select, sample_id, any_of(bcg_estro$genes))

  bcg_estro$expression <-
    map2(bcg_globals$assignment[names(bcg_estro$expression)],
         bcg_estro$expression,
         left_join,
         by = 'sample_id')

  ## data set-specific variable vectors

  bcg_estro$plot_variables <- bcg_estro$expression %>%
    map(names) %>%
    map(~bcg_estro$plot_variables[bcg_estro$plot_variables %in% .x])

# Classification of the genes of interest by specificity for the clusters ------

  insert_msg('Classification of the common significnt genes')

  ## done for the TCGA cohort and applied to the remaining collectives

  bcg_estro$classification_tcga <- bcg_estro$expression$tcga %>%
    classify(variables = bcg_estro$plot_variables$tcga,
             split_fct = 'clust_id') %>%
    .$classification

# Heat maps of expression Z scores ------

  insert_msg('Heat maps of expression Z-scores')

  bcg_estro$hm_plots <-
    list(data = bcg_estro$expression,
         plot_title = paste('Estrogen-responsive genes,',
                            globals$cohort_labs[names(bcg_estro$expression)]),
         variables = bcg_estro$plot_variables) %>%
    pmap(heat_map,
         split_fct = 'clust_id',
         normalize = TRUE,
         variable_classification = bcg_estro$classification_tcga,
         limits = c(-3, 3),
         oob = scales::squish,
         midpoint = 0,
         x_lab = 'cancer sample',
         cust_theme = globals$common_theme +
           theme(axis.title.y = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.text.y = element_markdown())) %>%
    map(~.x +
          scale_y_discrete(labels = set_names(bcg_estro$gene_lexicon$axis_lab,
                                              bcg_estro$gene_lexicon$gene_symbol)))

# Box plots for top up- and downregulated genes in the references data base -------

  insert_msg('Box plots')

  ## and 'canonical' ER-resonsive genes

  bcg_estro$box_variables <- bcg_estro$gene_lexicon %>%
    group_by(mode) %>%
    top_n(n = 10, abs(beta_score)) %>%
    ungroup %>%
    .$gene_symbol

  bcg_estro$box_variables <- c(bcg_estro$box_variables,
                               'AREG', 'PGR') %>%
    unique

  bcg_estro$box_subtitles <- bcg_estro$estimates %>%
    map(filter,
        gene_symbol %in% bcg_estro$box_variables,
        !duplicated(gene_symbol)) %>%
    map(mutate,
        gene_symbol = factor(gene_symbol, bcg_estro$box_variables)) %>%
    map(arrange, gene_symbol) %>%
    map(~set_names(.x$plot_cap,
                   .x$gene_symbol))

  ## box plots: ANOVA results are displayed in the plot subtitles

  for(i in names(bcg_estro$expression)) {

    bcg_estro$box_plots[[i]] <-
      list(variable = names(bcg_estro$box_subtitles[[i]]),
           plot_title = names(bcg_estro$box_subtitles[[i]]) %>%
             exchange(bcg_estro$gene_lexicon,
                      key = 'gene_symbol',
                      value = 'axis_lab') %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = bcg_estro$box_subtitles[[i]]) %>%
      pmap(plot_variable,
           bcg_estro$expression[[i]],
           split_factor = 'clust_id',
           type = 'box',
           cust_theme = globals$common_theme +
             theme(plot.title = element_markdown()),
           x_n_labs = TRUE,
           y_lab = expression('log'[2] * ' expression'),
           x_lab = 'cluster') %>%
      map(~.x +
            scale_fill_manual(values = globals$cluster_colors)) %>%
      set_names(names(bcg_estro$box_subtitles[[i]]))

  }

# END -----

  rm(i)

  insert_tail()
