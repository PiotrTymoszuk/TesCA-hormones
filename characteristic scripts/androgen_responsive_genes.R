# differential expression of androgen-responsive genes in the hormonal clusters
#
# Sources of the genes:
# EMBO J. 2011 Jul 6; 30(13): 2719–2733.
# doi: 10.1038/emboj.2011.158
#
# Table S7

  insert_head()

# container -------

  bcg_andro <- list()

# genes of interest ------

  insert_msg('Genes of interest, regulation estimates, and significant features')

  ## genes of interest:
  ## we're selecting genes with sustainable regulation sign, i.e.
  ## the sme mode of regulation at the early and late time points

  bcg_andro$gene_lexicon <-
    read_excel('./data/signatures/emboj2011158s8.xls',
               skip = 1)[1:5] %>%
    set_names(c('entrez_id',
                'gene_symbol',
                'early_regulation',
                'late_regulation',
                'auto_corr')) %>%
    filter(early_regulation == late_regulation) %>%
    mutate(mode = factor(early_regulation, c('up', 'down'))) %>%
    filter(complete.cases(.))

  bcg_andro$gene_lexicon <- bcg_andro$gene_lexicon %>%
    mutate(axis_lab = ifelse(mode == 'up',
                             paste0("<i style = 'color:#ff4500'>",
                                    gene_symbol, "</i>"),
                             paste0("<i style = 'color:#8b668b'>",
                                    gene_symbol, "</i>")))

  bcg_andro$genes <- bcg_andro$gene_lexicon$gene_symbol %>%
    unique

  # significant effects in single cohorts

  bcg_andro$significant <- bcg_dge$significant %>%
    map(map, ~.x[.x %in% bcg_andro$genes])

  ## common significant genes

  bcg_andro$common_significant <- bcg_dge$common_significant %>%
    map(~.x[.x %in% bcg_andro$genes])

  # regulation estimates for the genes of interest,
  ## appended with ANOVA results (to be shown in box plots)

  bcg_andro$estimates <- bcg_dge$dev_test %>%
    map(filter,
        gene_symbol %in% bcg_andro$genes) %>%
    map(select,
        clust_id,
        gene_symbol, entrez_id,
        deviation_center, lower_ci, upper_ci,
        regulation)

  bcg_andro$estimates <-
    map2(bcg_andro$estimates,
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

  bcg_andro$plot_variables <- bcg_andro$estimates[c("tcga", "gse99420")] %>%
    map(filter,
        effect_size >= 0.26,
        gene_symbol %in% unique(unname(unlist(bcg_andro$common_significant)))) %>%
    map(~.x$gene_symbol) %>%
    reduce(intersect)

  bcg_andro$expression <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    map(select, sample_id, any_of(bcg_andro$genes))

  bcg_andro$expression <-
    map2(bcg_globals$assignment[names(bcg_andro$expression)],
         bcg_andro$expression,
         left_join,
         by = 'sample_id')

  ## data set-specific variable vectors

  bcg_andro$plot_variables <- bcg_andro$expression %>%
    map(names) %>%
    map(~bcg_andro$plot_variables[bcg_andro$plot_variables %in% .x])

# Classification of the genes of interest by specificity for the clusters ------

  insert_msg('Classification of the common significnt genes')

  ## done for the TCGA cohort and applied to the remaining collectives

  bcg_andro$classification_tcga <- bcg_andro$expression$tcga %>%
    classify(variables = bcg_andro$plot_variables$tcga,
             split_fct = 'clust_id') %>%
    .$classification

# Heat maps of expression Z scores ------

  insert_msg('Heat maps of expression Z-scores')

  bcg_andro$hm_plots <-
    list(data = bcg_andro$expression,
         plot_title = paste('AR-responsive genes,',
                            globals$cohort_labs[names(bcg_andro$expression)]),
         variables = bcg_andro$plot_variables) %>%
    pmap(heat_map,
         split_fct = 'clust_id',
         normalize = TRUE,
         variable_classification = bcg_andro$classification_tcga,
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
          scale_y_discrete(labels = set_names(bcg_andro$gene_lexicon$axis_lab,
                                              bcg_andro$gene_lexicon$gene_symbol)))

# Box plots for top up- and downregulated genes in the references data base -------

  insert_msg('Box plots')

  bcg_andro$box_variables <- bcg_andro$gene_lexicon %>%
    group_by(mode) %>%
    top_n(n = 10, abs(auto_corr)) %>%
    ungroup %>%
    .$gene_symbol

  bcg_andro$box_subtitles <- bcg_andro$estimates %>%
    map(filter,
        gene_symbol %in% bcg_andro$box_variables,
        !duplicated(gene_symbol)) %>%
    map(mutate,
        gene_symbol = factor(gene_symbol, bcg_andro$box_variables)) %>%
    map(arrange, gene_symbol) %>%
    map(~set_names(.x$plot_cap,
                   .x$gene_symbol))

  ## box plots: ANOVA results are displayed in the plot subtitles

  for(i in names(bcg_andro$expression)) {

    bcg_andro$box_plots[[i]] <-
      list(variable = names(bcg_andro$box_subtitles[[i]]),
           plot_title = names(bcg_andro$box_subtitles[[i]]) %>%
             exchange(bcg_andro$gene_lexicon,
                      key = 'gene_symbol',
                      value = 'axis_lab') %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = bcg_andro$box_subtitles[[i]]) %>%
      pmap(plot_variable,
           bcg_andro$expression[[i]],
           split_factor = 'clust_id',
           type = 'box',
           cust_theme = globals$common_theme +
             theme(plot.title = element_markdown()),
           x_n_labs = TRUE,
           y_lab = expression('log'[2] * ' expression'),
           x_lab = 'cluster') %>%
      map(~.x +
            scale_fill_manual(values = globals$cluster_colors)) %>%
      set_names(names(bcg_andro$box_subtitles[[i]]))

  }

# END -----

  rm(i)

  insert_tail()
