# differential expression of testis- and cancer-specific antigens.
#
# Sources of the genes:
# Wang, C., Gu, Y., Zhang, K. et al.
# Systematic identification of genes with a cancer-testis expression pattern
# in 19 cancer types. Nat Commun 7, 10499 (2016).
# https://doi.org/10.1038/ncomms10499
# Supplementary Data 3A.
#
# http://www.cta.lncc.br
#
# The gene list is appended with the clinically relevant
# testicular cancer markers.

  insert_head()

# container -------

  bcg_cta <- list()

# genes of interest ------

  insert_msg('Genes of interest, regulation estimates, and significant features')

  bcg_cta$genes$wang <-
    read_xlsx(path = './data/signatures/41467_2016_BFncomms10499_MOESM726_ESM.xlsx',
              sheet = 2,
              skip = 1) %>%
    .$Description

  bcg_cta$genes$ctpedia <-
    read_html('./data/signatures/CTpedia.htm') %>%
    html_table %>%
    .[[1]] %>%
    .[['Family member']]

  bcg_cta$genes$ctpedia <-
    stri_split_fixed(bcg_cta$genes$ctpedia,
                     pattern = '/',
                     simplify = TRUE)[, 1]

  bcg_cta$genes <- bcg_cta$genes %>%
    reduce(union)

  bcg_cta$genes <- c(bcg_cta$genes,
                     'AFP',
                     'CGA',
                     paste('CGB', c(1, 2, 3, 5, 7, 8)),
                     'LDHA', 'LDHB', 'LDHC')

  # significant effects in single cohorts

  bcg_cta$significant <- bcg_dge$significant %>%
    map(map, ~.x[.x %in% bcg_cta$genes])

  ## common significant genes

  bcg_cta$common_significant <- bcg_dge$common_significant %>%
    map(~.x[.x %in% bcg_cta$genes])

  # regulation estimates for the genes of interest,
  ## appended with ANOVA results (to be shown in box plots)

  bcg_cta$estimates <- bcg_dge$dev_test %>%
    map(filter,
        gene_symbol %in% bcg_cta$genes) %>%
    map(select,
        clust_id,
        gene_symbol, entrez_id,
        deviation_center, lower_ci, upper_ci,
        regulation)

  bcg_cta$estimates <-
    map2(bcg_cta$estimates,
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

  bcg_cta$plot_variables <- bcg_cta$estimates[c("tcga", "gse99420")] %>%
    map(filter,
        effect_size >= 0.26,
        gene_symbol %in% unique(unname(unlist(bcg_cta$common_significant)))) %>%
    map(~.x$gene_symbol) %>%
    reduce(intersect)

  bcg_cta$expression <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    map(select, sample_id, any_of(bcg_cta$genes))

  bcg_cta$expression <-
    map2(bcg_globals$assignment[names(bcg_cta$expression)],
         bcg_cta$expression,
         left_join,
         by = 'sample_id')

  ## data set-specific variable vectors

  bcg_cta$plot_variables <- bcg_cta$expression %>%
    map(names) %>%
    map(~bcg_cta$plot_variables[bcg_cta$plot_variables %in% .x])

# Classification of the genes of interest by specificity for the clusters ------

  insert_msg('Classification of the common significnt genes')

  ## done for the TCGA cohort and applied to the remaining collectives

  bcg_cta$classification_tcga <- bcg_cta$expression$tcga %>%
    classify(variables = bcg_cta$plot_variables$tcga,
             split_fct = 'clust_id') %>%
    .$classification

# Heat maps of expression Z scores ------

  insert_msg('Heat maps of expression Z-scores')

  bcg_cta$hm_plots <-
    list(data = bcg_cta$expression,
         plot_title = paste('CTA and cancer marker genes,',
                            globals$cohort_labs[names(bcg_cta$expression)]),
         variables = bcg_cta$plot_variables) %>%
    pmap(heat_map,
         split_fct = 'clust_id',
         normalize = TRUE,
         variable_classification = bcg_cta$classification_tcga,
         limits = c(-3, 3),
         oob = scales::squish,
         midpoint = 0,
         x_lab = 'cancer sample',
         cust_theme = globals$common_theme +
           theme(axis.title.y = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.text.y = element_text(face = 'italic'))) %>%
    map(~.x +
          guides(y = guide_axis(n.dodge = 4)))

# Forest plots for selected CTA families --------

  insert_msg('Forest plots')

  ## variables of interest: MAGE, PRAME, TSPY, TEX, and HORMAD families

  bcg_cta$forest_variables <- bcg_cta$common_significant %>%
    unlist %>%
    unname %>%
    unique

  bcg_cta$forest_variables <-
    bcg_cta$forest_variables[stri_detect(bcg_cta$forest_variables,
                                         regex = '^(MAGE|PRAME|TSPY|TEX|HORMAD)')] %>%
    tibble(gene_symbol = .) %>%
    mutate(class = stri_extract(gene_symbol,
                                regex = '^(MAGE|PRAME|TSPY|TEX|HORMAD)'))

  bcg_cta$estimates <- bcg_cta$estimates %>%
    map(left_join,
        bcg_cta$forest_variables,
        by = 'gene_symbol')

  ## Forest plots

  bcg_cta$forest_plots <-
    list(x = bcg_cta$estimates %>%
           map(filter, gene_symbol %in% bcg_cta$forest_variables$gene_symbol),
         y = paste('CTA genes,',
                   globals$cohort_labs[names(bcg_cta$estimates)])) %>%
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

# Box plots for selected, clinically relevant cancer markers -------

  insert_msg('Box plots')

  bcg_cta$box_variables <- c('AFP',
                             'CGA',
                             paste0('CGB', c(1, 2, 3, 5, 7, 8)),
                             'LDHA', 'LDHB', 'LDHC')

  bcg_cta$box_subtitles <- bcg_cta$estimates %>%
    map(filter,
        gene_symbol %in% bcg_cta$box_variables,
        !duplicated(gene_symbol)) %>%
    map(mutate,
        gene_symbol = factor(gene_symbol, bcg_cta$box_variables)) %>%
    map(arrange, gene_symbol) %>%
    map(~set_names(.x$plot_cap,
                   .x$gene_symbol))

  ## box plots: ANOVA results are displayed in the plot subtitles

  for(i in names(bcg_cta$expression)) {

    bcg_cta$box_plots[[i]] <-
      list(variable = names(bcg_cta$box_subtitles[[i]]),
           plot_title = names(bcg_cta$box_subtitles[[i]]) %>%
             html_italic %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = bcg_cta$box_subtitles[[i]]) %>%
      pmap(plot_variable,
           bcg_cta$expression[[i]],
           split_factor = 'clust_id',
           type = 'box',
           cust_theme = globals$common_theme +
             theme(plot.title = element_markdown()),
           x_n_labs = TRUE,
           y_lab = expression('log'[2] * ' expression'),
           x_lab = 'cluster') %>%
      map(~.x +
            scale_fill_manual(values = globals$cluster_colors)) %>%
      set_names(names(bcg_cta$box_subtitles[[i]]))

  }

# END -----

  rm(i)

  bcg_cta <- bcg_cta[c("genes", "significant", "common_significant",
                       "plot_variables", "hm_plots",
                       "forest_variables", "forest_plots",
                       "box_variables", "box_plots")]

  insert_tail()
