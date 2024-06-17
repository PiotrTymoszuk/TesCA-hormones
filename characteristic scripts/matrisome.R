# Differential expression of matrisome genes

  insert_head()

# container ------

  bcg_ecm <- list()

# analysis globals: matrisome genes and their DGE estimates ------

  insert_msg('Matrisome genes')

  ## regulation estimates, appending with the ANOVA results
  ## The ANOVA results will be used for identification of the top
  ## differentially expressed genes in matrisome categories and clusters

  bcg_ecm$estimates <- bcg_dge$dev_test %>%
    map(select,
        clust_id,
        gene_symbol,
        entrez_id,
        deviation_center,
        lower_ci,
        upper_ci,
        regulation) %>%
    map(as.data.frame) %>%
    map(matriannotate,
        gene.column = 'gene_symbol',
        species = 'human') %>%
    map(filter,
        `Annotated Matrisome Category` != 'Non-matrisome') %>%
    map(mutate,
        gene_symbol = `Annotated Gene`) %>%
    map(as_tibble)

  bcg_ecm$estimates <-
    map2(bcg_ecm$estimates,
         map(bcg_dge$anova,
             ~.x[c('gene_symbol', 'effect_size', 'p_adjusted')]),
         left_join,
         by = 'gene_symbol')

  ## matrisome genes and their classification

  bcg_ecm$genes <- bcg_ecm$estimates %>%
    map(~.x$gene_symbol) %>%
    reduce(union)

  bcg_ecm$classification <- bcg_ecm$estimates %>%
    map_dfr(~.x[c('gene_symbol',
                  'Annotated Matrisome Category')]) %>%
    set_names(c('variable', 'ecm_class')) %>%
    filter(!duplicated(variable)) %>%
    mutate(ecm_class = factor(ecm_class,
                              c('Collagens',
                                'ECM Glycoproteins',
                                'Proteoglycans',
                                'ECM-affiliated Proteins',
                                'ECM Regulators',
                                'Secreted Factors')))

  # significant effects in single cohorts

  bcg_ecm$significant <- bcg_dge$significant %>%
    map(map, ~.x[.x %in% bcg_ecm$genes])

  ## common significant genes

  bcg_ecm$common_significant <- bcg_dge$common_significant %>%
    map(~.x[.x %in% bcg_ecm$genes])

# Plotting variables and expression values for heat maps ------

  insert_msg('Plotting variables and expression values')

  ## plotting variables: genes found to be regulated in at least one cluster
  ## in both the TCGA and GSE99420 cohort

  bcg_ecm$plot_variables <- bcg_ecm$estimates[c("tcga", "gse99420")] %>%
    map(filter,
        gene_symbol %in% unique(unname(unlist(bcg_ecm$common_significant)))) %>%
    map(~.x$gene_symbol) %>%
    reduce(intersect)

  bcg_ecm$expression <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    map(select, sample_id, any_of(bcg_ecm$plot_variables))

  bcg_ecm$expression <-
    map2(bcg_globals$assignment[names(bcg_ecm$expression)],
         bcg_ecm$expression,
         left_join,
         by = 'sample_id')

  ## data set-specific variable vectors

  bcg_ecm$plot_variables <- bcg_ecm$expression %>%
    map(names) %>%
    map(~bcg_ecm$plot_variables[bcg_ecm$plot_variables %in% .x])

# Classification of the genes of interest by specificity for the clusters ------

  insert_msg('Classification of the common significnt genes')

  ## done for the TCGA cohort and applied to the remaining collectives

  bcg_ecm$classification_tcga <- bcg_ecm$expression$tcga %>%
    classify(variables = bcg_ecm$plot_variables$tcga,
             split_fct = 'clust_id') %>%
    .$classification

  ## appending the classification with the ECM division information

  bcg_ecm$classification_tcga <-
    left_join(bcg_ecm$classification_tcga,
              bcg_ecm$classification,
              by = 'variable') %>%
    arrange(clust_id, desc(ecm_class), delta_auc)

  bcg_ecm$classification_tcga <-
    bcg_ecm$classification_tcga %>%
    mutate(variable = factor(variable,
                             bcg_ecm$classification_tcga$variable))

# Scales for the heat maps: top regulated genes -------

  insert_msg('Scales for the heat maps')

  ## top regulated genes per cluster in the TCGA cohort
  ## note that the label order does not correspond to
  ## the actual gene location on the Y scale.
  ## The purpose is just to show the top cluster markers

  bcg_ecm$cluster_scale <- bcg_ecm$classification_tcga %>%
    blast(clust_id) %>%
    map2_dfr(., c(77, 13, 21, 11),
             ~top_n(.x, n = .y, auc)) %>%
    .$variable %>%
    as.character

  set.seed(12345)

  bcg_ecm$cluster_scale <- bcg_ecm$classification_tcga %>%
    blast(clust_id) %>%
    map2(., globals$cluster_hex_colors,
         ~mutate(.x,
                 axis_lab = ifelse(variable %in% bcg_ecm$cluster_scale,
                                   as.character(variable),
                                   NA),
                 axis_lab = space_evenly(axis_lab),
                 axis_lab = ifelse(is.na(axis_lab),
                                   '',
                                   paste0("<i style = 'color:", .y, "'>",
                                          axis_lab, "</i>")))) %>%
    map(~set_names(.x$axis_lab, as.character(.x$variable))) %>%
    reduce(c)

# Heat maps of expression Z scores: main plots ------

  insert_msg('Heat map plots')

  ## the genes are classified by cluster specificity and matrisome division

  bcg_ecm$hm_plots <-
    list(data = bcg_ecm$expression,
         plot_title = paste('Matrisome genes,',
                            globals$cohort_labs[names(bcg_ecm$expression)]),
         variables = bcg_ecm$plot_variables) %>%
    pmap(heat_map,
         split_fct = 'clust_id',
         normalize = TRUE,
         variable_classification = bcg_ecm$classification_tcga,
         limits = c(-3, 3),
         oob = scales::squish,
         midpoint = 0,
         x_lab = 'cancer sample',
         cust_theme = globals$common_theme +
           theme(axis.title.y = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.text.y = element_markdown())) %>%
    map(~.x +
          scale_y_discrete(labels = bcg_ecm$cluster_scale) +
          guides(y = guide_axis(n.dodge = 4)))

# Rug heat maps with classification of matrisome genes ------

  insert_msg('Rug heat maps')

  bcg_ecm$hm_rug <-
    bcg_ecm$classification_tcga %>%
    ggplot(aes(x = 'gene',
               y = variable,
               fill = ecm_class)) +
    facet_grid(clust_id ~ .,
               scales = 'free',
               space = 'free') +
    geom_tile() +
    scale_fill_viridis_d(name = 'Gene\nclassification') +
    theme_void() +
    theme(strip.text.y = element_blank(),
          legend.text = globals$common_text,
          legend.title = globals$common_text)

# END ----

  bcg_ecm <- bcg_ecm[c("genes",
                       "classification", "classification_tcga",
                       "significant", "common_significant",
                       "plot_variables",
                       "cluster_scale",
                       "hm_plots", "hm_rug")]

  insert_tail()
