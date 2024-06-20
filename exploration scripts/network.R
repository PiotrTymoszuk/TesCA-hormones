# Co-expression analysis for the hormone-related genes done with a graph of
# Spearman's correlations (edges: correlations rho > 0.3).
#
# CombBat-adjusted log2-transformed expression of the genes passing variability
# and minimal expression criteria is included in the analysis.
#
# Potential hubs of hormone metabolism are identified by a visual analysis
# of graph plots, and comparison of the gene's degree, betweenness
# and hub scores.

  insert_head()

# container -------

  expl_net <- list()

# analysis variables and data ------

  insert_msg('Analysis variables and data')

  expl_net$variables <- expl_dist$top_variables

  expl_net$data <- combat$expression %>%
    map(select,
        sample_id, all_of(expl_net$variables)) %>%
    map(column_to_rownames, 'sample_id') %>%
    map(center_data)

# Similarity graphs and minimum spanning trees -----------

  insert_msg('Similarity graphs')

  expl_net$graph_objects <- expl_net$data %>%
    map(simil_graph,
        method = 'spearman',
        corr_cutoff = 0.3,
        diag = FALSE)

  expl_net$mst_objects <- expl_net$graph_objects %>%
    map(mst)

# Graph statistics ---------

  insert_msg('Graph statistics')

  expl_net$stats <- expl_net$graph_objects %>%
    map(get_graph_stats,
        directed = FALSE,
        normalized = TRUE)

  ## top most important nodes to be presented in the plots

  expl_net$top_genes <- expl_net$stats %>%
    map(~list(slice_max(.x, hub_score, n = 15),
              slice_max(.x, betweenness, n = 15),
              slice_max(.x, degree, n = 15))) %>%
    map(map, ~.x$variable) %>%
    map(reduce, union)

  expl_net$stats <-
    map2(expl_net$stats,
         expl_net$top_genes,
         ~mutate(.x,
                 plot_label = ifelse(variable %in% .y,
                                     variable, NA)))

  ## plots of degree, hub scores, and betweenness

  expl_net$stat_plots <-
    list(x = expl_net$stats,
         y = globals$cohort_labs[names(expl_net$stats)]) %>%
    pmap(function(x, y) x %>%
           ggplot(aes(x = degree,
                      y = hub_score,
                      fill = class,
                      size = betweenness)) +
           geom_point(color = 'black',
                      shape = 21) +
           geom_text_repel(aes(label = plot_label),
                           size = 2.3,
                           fontface = 'italic') +
           scale_fill_manual(values = globals$gene_class_colors,
                             name = 'Gene\nclassification') +
           scale_size_area(max_size = 4,
                           limits = expl_net$stats %>%
                             map(~.x$betweenness) %>%
                             reduce(c) %>%
                             range,
                           name = 'Normalized\nbetweenness') +
           globals$common_theme +
           labs(title = y,
                subtitle = paste('genes: n =', nrow(x)),
                x = 'degree',
                y = 'normalized hub score'))


# Visualization of the graphs --------

  insert_msg('Visualiztion of the graphs and minimal spanning trees')

  set.seed(1234)

  expl_net$graph_plots <-
    list(graph_object = expl_net$graph_objects,
         plot_title = globals$cohort_labs[names(expl_net$graph_objects)],
         select_nodes = expl_net$top_genes) %>%
    pmap(plot_simil_graph,
         txt_size = 2,
         cust_theme = globals$common_theme +
           theme(panel.grid.major = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.text.x = element_blank(),
                 axis.text.y = element_blank(),
                 axis.title = element_blank()))

  expl_net$mst_plots <-
    list(graph_object = expl_net$mst_objects,
         plot_title = globals$cohort_labs[names(expl_net$mst_objects)]) %>%
    pmap(plot_simil_graph,
         txt_size = 2,
         cust_theme = globals$common_theme +
           theme(panel.grid.major = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.text.x = element_blank(),
                 axis.text.y = element_blank(),
                 axis.title = element_blank()))

# END ------

  expl_net$data <- NULL
  expl_net$variables <- NULL

  expl_net <- compact(expl_net)

  insert_tail()
