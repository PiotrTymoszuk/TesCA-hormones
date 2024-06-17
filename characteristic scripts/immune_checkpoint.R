# Differential gene expression of genes associated with immune checkpoint
# and T cell exhaustion in the hormonal clusters.
#
# Differential expression compared with the cohort mean is evaluated. P values
# refer to whole-genome effects, i.e. are corrected by FDR at the whole-genome
# scale

  insert_head()

# container ------

  bcg_tex <- list()

# published genes -------

  insert_msg('Published genes')

  ## signatures of bcg_texion
  ## sources:
  ## https://doi.org/10.1016/j.ebiom.2022.104207
  ## https://doi.org/10.1158/1078-0432.CCR-22-0275
  ## https://doi.org/10.1158/1078-0432.CCR-17-1846
  ## https://doi.org/10.1002/path.5435
  ## There is a discrepant gene VSIR (Entrez ID: 64115)
  ## which has not official symbol in the org.Hs database
  ## but corresponds to the C10orf54 item

  bcg_tex$genes <- c('BTLA',
                     'CTLA4',
                     'CD274',
                     'PDCD1',
                     'PDCD1LG2',
                     'HAVCR2',
                     'IDO1',
                     'SIGLEC7',
                     'TIGIT',
                     'LAG3',
                     'VSIR',
                     'CXCR4',
                     'CD160',
                     'ENTPD1',
                     'CD244',
                     'EOMES',
                     'CD274',
                     'CD200',
                     'NOS2',
                     'TOX',
                     'TOX2',
                     'C10orf54')

  bcg_tex$genes <-
    read_excel('./data/signatures/path5435-sup-0002-supptables.xlsx',
               sheet = 'S5',
               skip = 3) %>%
    .$gene %>%
    c(bcg_tex$gene) %>%
    unique

  bcg_tex$genes <- bcg_tex$genes[bcg_tex$genes != 'VSIR']

  ## gene labels

  bcg_tex$genes <- ifelse(bcg_tex$genes == 'C10orf54',
                          'VSIR', bcg_tex$genes) %>%
    set_names(bcg_tex$genes) %>%
    compress(names_to = 'gene_symbol',
             values_to = 'gene_label') %>%
    filter(!duplicated(gene_label))

# Estimates of differential expression ------

  insert_msg('Estimates of differential expression')

  ## significant effects in single cohorts

  bcg_tex$significant <- bcg_dge$significant %>%
    map(map, ~.x[.x %in% bcg_tex$genes$gene_symbol])

  ## common significant genes

  bcg_tex$common_significant <- bcg_dge$common_significant %>%
    map(~.x[.x %in% bcg_tex$genes$gene_symbol])

  # regulation estimates for the genes of interest,
  ## appended with ANOVA results (to be shown in box plots)

  bcg_tex$estimates <- bcg_dge$dev_test %>%
    map(filter,
        gene_symbol %in% bcg_tex$genes$gene_symbol) %>%
    map(select,
        clust_id,
        gene_symbol, entrez_id,
        deviation_center, lower_ci, upper_ci,
        regulation)

  bcg_tex$estimates <-
    map2(bcg_tex$estimates,
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

  bcg_tex$plot_variables <- bcg_tex$common_significant %>%
    unlist %>%
    unname %>%
    unique

  bcg_tex$expression <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    map(select, sample_id, any_of(bcg_tex$genes$gene_symbol))

  bcg_tex$expression <-
    map2(bcg_globals$assignment[names(bcg_tex$expression)],
         bcg_tex$expression,
         left_join,
         by = 'sample_id')

  ## data set-specific variable vectors

  bcg_tex$plot_variables <- bcg_tex$expression %>%
    map(names) %>%
    map(~bcg_tex$plot_variables[bcg_tex$plot_variables %in% .x])

# Classification of the genes of interest by specificity for the clusters ------

  insert_msg('Classification of the common significnt genes')

  ## done for the TCGA cohort and applied to the remaining collectives

  bcg_tex$classification_tcga <- bcg_tex$expression$tcga %>%
    classify(variables = bcg_tex$plot_variables$tcga,
             split_fct = 'clust_id') %>%
    .$classification

# Heat maps of expression Z scores ------

  insert_msg('Heat maps of expression Z-scores')

  bcg_tex$hm_plots <-
    list(data = bcg_tex$expression,
         plot_title = paste('IC- and TEx- associated genes,',
                            globals$cohort_labs[names(bcg_tex$expression)]),
         variables = bcg_tex$plot_variables) %>%
    pmap(heat_map,
         split_fct = 'clust_id',
         normalize = TRUE,
         variable_classification = bcg_tex$classification_tcga,
         limits = c(-3, 3),
         oob = scales::squish,
         midpoint = 0,
         x_lab = 'cancer sample',
         cust_theme = globals$common_theme +
           theme(axis.title.y = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.text.y = element_text(face = 'italic')))

# Forest plots of the regulation estimates -------

  insert_msg('Forest plots')

  bcg_tex$estimates <- bcg_tex$estimates %>%
    map(left_join,
        set_names(bcg_tex$classification_tcga[c("variable", "clust_id")],
                  c('gene_symbol', 'class')),
        by = 'gene_symbol')

  bcg_tex$forest_plots <-
    list(x = bcg_tex$estimates %>%
           map(filter, gene_symbol %in% reduce(bcg_tex$plot_variables, union)),
         y =  paste('IC- and TEx- associated genes,',
                    globals$cohort_labs[names(bcg_tex$expression)])) %>%
    pmap(function(x, y) x %>%
           ggplot(aes(x = deviation_center,
                      y = gene_symbol,
                      color = regulation)) +
           facet_grid(class ~ clust_id,
                      space = 'free_y',
                      scales = 'free_y') +
           geom_vline(xintercept = 0,
                      linetype = 'dashed') +
           geom_errorbarh(aes(xmin = lower_ci,
                              xmax = upper_ci),
                          height = 0) +
           geom_point(shape = 16,
                      size = 2) +
           scale_color_manual(values = c(upregulated = 'firebrick',
                                         downregulated = 'steelblue',
                                         ns = 'gray70'),
                              name = 'Regulation\nvs cohort mean') +
           globals$common_theme +
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_text(face = 'italic'),
                 strip.background.y = element_blank(),
                 strip.text.y = element_blank()) +
           labs(title = y,
                x = expression('log'[2] * ' fold-regulation vs cohort mean')))

# Box plots for selected, clinically relevant genes ------

  insert_msg('Box plots for selected genes')

  bcg_tex$box_variables <- c('CD274', 'CXCR4', 'CTLA4', 'HAVCR2',
                             'IDO1', 'PDCD1', 'PDCD1LG2',
                             'TIGIT', 'NOS2')

  bcg_tex$box_subtitles <- bcg_tex$estimates %>%
    map(filter,
        gene_symbol %in% bcg_tex$box_variables,
        !duplicated(gene_symbol)) %>%
    map(mutate,
        gene_symbol = factor(gene_symbol, bcg_tex$box_variables)) %>%
    map(arrange, gene_symbol) %>%
    map(~.x$plot_cap)

  ## box plots: ANOVA results are displayed in the plot subtitles

  for(i in names(bcg_tex$expression)) {

    bcg_tex$box_plots[[i]] <-
      list(variable = bcg_tex$box_variables,
           plot_title = bcg_tex$box_variables %>%
             html_italic %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = bcg_tex$box_subtitles[[i]]) %>%
      pmap(plot_variable,
           bcg_tex$expression[[i]],
           split_factor = 'clust_id',
           type = 'box',
           cust_theme = globals$common_theme +
             theme(plot.title = element_markdown()),
           x_n_labs = TRUE,
           y_lab = expression('log'[2] * ' expression'),
           x_lab = 'cluster') %>%
      map(~.x +
            scale_fill_manual(values = globals$cluster_colors)) %>%
      set_names(bcg_tex$box_variables)

  }

# END -----

  rm(i)

  bcg_tex$estimates <- NULL
  bcg_tex$expression <- NULL
  bcg_tex$box_variables <- NULL
  bcg_tex$box_subtitles <- NULL

  bcg_tex <- compact(bcg_tex)

  insert_tail()
