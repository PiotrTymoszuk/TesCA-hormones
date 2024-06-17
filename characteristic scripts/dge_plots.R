# Plots for results of differential gene expression

  insert_head()

# container ------

  bcg_dgeplots <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  ## Y axis labels with the total N numbers of genes

  bcg_dgeplots$n_axis_labs <- bcg_dge$dge_numbers %>%
    filter(cohort != 'common',
           !duplicated(cohort))

  bcg_dgeplots$n_axis_labs <-
    map2(bcg_dgeplots$n_axis_labs$cohort,
         bcg_dgeplots$n_axis_labs$n_total,
         ~paste(globals$cohort_labs[.x],
                .y,
                sep = '\nn = ')) %>%
    set_names(bcg_dgeplots$n_axis_labs$cohort)

  ## cluster names

  bcg_dgeplots$cluster_levels <-
    levels(bcg_dge$dev_test[[1]]$clust_id)

  ## common differentially expressed genes and their expression values
  ## (i.e. differentially regulated genes shared by TCGA and GSE99420)

  bcg_dgeplots$common_genes <- bcg_dge$common_significant %>%
    reduce(union)

  bcg_dgeplots$expression <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    map(select, sample_id, any_of(bcg_dgeplots$common_genes))

  bcg_dgeplots$expression <-
    map2(bcg_globals$assignment[names(bcg_dgeplots$expression)],
         bcg_dgeplots$expression,
         inner_join, by = 'sample_id') %>%
    map(column_to_rownames, 'sample_id')

  ## top common regulated genes in the clusters
  ## to be displayed at the heat map axes

  bcg_dgeplots$top_common_genes <- bcg_dge$dev_test[c("tcga", "gse99420")] %>%
    map_dfr(filter,
            gene_symbol %in% bcg_dgeplots$common_genes) %>%
    select(clust_id, gene_symbol, deviation_center) %>%
    group_by(clust_id, gene_symbol) %>%
    summarise(mean_deviation = mean(deviation_center)) %>%
    ungroup %>%
    group_by(clust_id) %>%
    top_n(n = 10, mean_deviation) %>%
    ungroup

# Numbers of differentially regulated genes --------

  insert_msg('Numbers of differentially regulated genes')

  ## bar plots

  bcg_dgeplots$n_numbers <- bcg_dge$dge_numbers %>%
    filter(cohort != 'common') %>%
    mutate(percent = n/n_total * 100,
           percent_plot = ifelse(regulation == 'downregulated',
                                 -percent, percent)) %>%
    ggplot(aes(x = percent_plot,
               y = factor(cohort, rev(globals$analysis_cohorts)),
               fill = regulation)) +
    geom_vline(xintercept = 0,
               linetype = 'dashed') +
    geom_bar(stat = 'identity',
             color = 'black') +
    facet_grid(. ~ clust_id) +
    scale_x_continuous(labels = function(x) abs(x)) +
    scale_y_discrete(labels = bcg_dgeplots$n_axis_labs) +
    scale_fill_manual(values = c(upregulated = 'firebrick',
                                 downregulated = 'steelblue',
                                 name = 'Regulation\nvs cohort mean')) +
    globals$common_theme +
    theme(axis.title.y = element_blank()) +
    labs(title = 'Differentially regulated genes',
         subtitle = 'Cluster versus cohort mean',
         x = '% of total genes')

# Volcano plots ---------

  insert_msg('Volcano plots')

  ## log2 fold-regulation of regulation of gene expression in particular clusters
  ## as compared with the cohort mean and significance

  ## volcano plots: only ANOVA-significant genes are presented

  for(i in bcg_dgeplots$cluster_levels) {

    bcg_dgeplots$volcano_plots[[i]] <-
      list(data = map(bcg_dge$dev_test,
                      filter,
                      anova_significant == 'yes',
                      clust_id %in% i,
                      p_adjusted > 0),
           plot_title = paste(i,
                              globals$cohort_labs[names(bcg_dge$dev_test)],
                              sep = ' vs cohort mean, ')) %>%
      pmap(plot_volcano,
           regulation_variable = 'deviation_center',
           p_variable = 'p_adjusted',
           regulation_level = 0,
           x_lab = "log<sub>2</sub> fold-regulation vs cohort mean",
           y_lab = expression('-log'[10] * ' pFDR'),
           top_regulated = 20,
           label_variable = 'gene_symbol',
           label_type = 'text',
           txt_size = 2.3,
           txt_face = 'italic',
           cust_theme = globals$common_theme +
             theme(axis.title.x = element_markdown())) %>%
      map(~.x +
            labs(subtitle = .x$labels$tag) +
            theme(plot.tag = element_blank()))

  }

# Forest plots with the top regulated genes in each cluster ------

  insert_msg('Forest plots for top regulated genes')

  for(i in bcg_dgeplots$cluster_levels) {

    bcg_dgeplots$top_plots[[i]] <-
      list(data = map(bcg_dge$dev_test,
                      filter,
                      clust_id %in% i,
                      regulation %in% c('upregulated', 'downregulated')),
           plot_title = paste(i,
                              globals$cohort_labs[names(bcg_dge$dev_test)],
                              sep = ' vs cohort mean, ')) %>%
      pmap(plot_top,
           regulation_variable = 'deviation_center',
           label_variable = 'gene_symbol',
           p_variable = 'p_adjusted',
           top_regulated = 20,
           lower_ci_variable = 'lower_ci',
           upper_ci_variable = 'upper_ci',
           fill_title = 'Regulation\nvs cohort mean',
           cust_theme = globals$common_theme +
             theme(axis.title.x = element_markdown(),
                   axis.text.y = element_text(face = 'italic')),
           x_lab = "log<sub>2</sub> fold-regulation vs cohort mean")

  }

# Heat maps of Z-scores of expression of regulated genes --------

  insert_msg('Heat maps of expression of common regulated genes')

  ## classification of the common regulated genes in the TCGA cohort

  bcg_dgeplots$classification_tcga <-
    classify(bcg_dgeplots$expression$tcga,
             variables = bcg_dgeplots$common_genes,
             split_fct = 'clust_id')

  bcg_dgeplots$classification_tcga <-
    bcg_dgeplots$classification_tcga$classification %>%
    mutate(gene_label = ifelse(variable %in% bcg_dgeplots$top_common_genes$gene_symbol,
                               variable, NA))

  ## heat map plots

  bcg_dgeplots$hm_plots <-
    list(x = bcg_dgeplots$expression,
         y = paste('Common DGE,',
                   globals$cohort_labs[names(bcg_dgeplots$expression)])) %>%
    pmap(function(x, y) x %>%
           heat_map(variables = bcg_dgeplots$common_genes[bcg_dgeplots$common_genes %in% names(x)],
                    split_fct = 'clust_id',
                    normalize = TRUE,
                    variable_classification = bcg_dgeplots$classification_tcga %>%
                      filter(variable %in% names(x)),
                    plot_title = y,
                    cust_theme = globals$common_theme +
                      theme(axis.text.x = element_blank(),
                            axis.ticks = element_blank(),
                            axis.title.y = element_blank()),
                    limits = c(-3, 3),
                    midpoint = 0,
                    oob = scales::squish,
                    x_lab = 'cancer sample',
                    y_lab = 'gene'))

  ## appending with symbols of the top regulated genes per cluster.
  ## I'm re-shuffling the symbols within the cluster, to make them visible
  ## in the plots

  bcg_dgeplots$top_scale <- bcg_dgeplots$classification_tcga %>%
    blast(clust_id) %>%
    map_dfr(mutate,
            axis_label = space_evenly(gene_label),
            axis_label = ifelse(is.na(axis_label), '', axis_label),
            clust_color = globals$cluster_hex_colors[clust_id],
            axis_label = ifelse(axis_label != '',
                                paste0("<i style = 'color:",
                                       clust_color, "'>",
                                       axis_label,
                                       '</i>'),
                                axis_label))

  bcg_dgeplots$top_scale <-
    set_names(bcg_dgeplots$top_scale$axis_label,
              bcg_dgeplots$top_scale$variable)

  bcg_dgeplots$hm_plots <- bcg_dgeplots$hm_plots %>%
    map(~.x +
          guides(y = guide_axis(n.dodge = 4)) +
          scale_y_discrete(labels = bcg_dgeplots$top_scale) +
          theme(axis.text.y = element_markdown(size = 7,
                                               face = 'italic')))

# Forest plots of the common top regulated genes in the clusters --------

  insert_msg('Top common genes')

  ## only for the TCGA and GSE99420 cohorts

  ## top shared variables per cluster and their regulation estimates

  bcg_dgeplots$cmm_variables <-
    list(bcg_dge$common_significant[c("#1.upregulated", "#1.downregulated")],
         bcg_dge$common_significant[c("#2.upregulated", "#2.downregulated")],
         bcg_dge$common_significant[c("#3.upregulated", "#3.downregulated")],
         bcg_dge$common_significant[c("#4.upregulated", "#4.downregulated")]) %>%
    map(reduce, union) %>%
    set_names(levels(bcg_globals$assignment[[1]]$clust_id))

  bcg_dgeplots$top_cmm_estimates <-
    bcg_dge$dev_test[c("tcga", "gse99420")] %>%
    map(blast, clust_id) %>%
    transpose %>%
    map(compress, names_to = 'cohort') %>%
    map2(., bcg_dgeplots$cmm_variables,
         ~filter(.x, variable %in% .y)) %>%
    map(select,
        cohort,
        variable, regulation,
        deviation_center,
        lower_ci, upper_ci) %>%
    map(group_by, cohort, regulation) %>%
    map(top_n, n = 30, abs(deviation_center)) %>%
    map(ungroup)

  bcg_dgeplots$cmm_variables <- bcg_dgeplots$top_cmm_estimates %>%
    map(blast, regulation) %>%
    map(map, blast, cohort) %>%
    map(map, map, ~.x$variable) %>%
    map(map, reduce, intersect) %>%
    map(reduce, union) %>%
    map(~.x[!.x %in% globals$genes])

  bcg_dgeplots$top_cmm_estimates <-
    map2(bcg_dgeplots$top_cmm_estimates,
         bcg_dgeplots$cmm_variables,
         ~filter(.x, variable %in% .y)) %>%
    map(mutate, cohort = factor(cohort, c('tcga', 'gse99420')))

  ## Forest plots

  bcg_dgeplots$cmm_top_plots <-
    list(x = bcg_dgeplots$top_cmm_estimates,
         y = paste('Top regulated genes, cluster',
                   names(bcg_dgeplots$top_cmm_estimates))) %>%
    pmap(function(x, y) x %>%
           ggplot(aes(x = deviation_center,
                      y = reorder(variable, deviation_center),
                      color = regulation)) +
           facet_grid(regulation ~ cohort,
                      labeller = labeller(.cols = globals$cohort_labs),
                      space = 'free_y',
                      scale = 'free_y') +
           geom_vline(xintercept = 0,
                      linetype = 'dashed') +
           geom_errorbarh(aes(xmin = lower_ci,
                              xmax = upper_ci),
                          height = 0) +
           geom_point(size = 2,
                      shape = 16) +
           scale_color_manual(values = c(upregulated = 'firebrick',
                                         downregulated = 'steelblue',
                                         ns = 'gray60'),
                              name = 'Regulation\nvs cohort mean') +
           globals$common_theme +
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_text(face = 'italic')) +
           labs(title = y,
                x = expression('log'[2] * ' fold-regulation vs cohort mean')))

# END ------

  bcg_dgeplots <-
    bcg_dgeplots[c('n_numbers',
                   'volcano_plots',
                   'top_plots',
                   'top_scale',
                   'hm_plots',
                   'cmm_top_plots')]

  insert_tail()
