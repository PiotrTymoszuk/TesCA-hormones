# Co-regulation networks of proteins in the hormonal clusters of the TCGA
# cohort.

  insert_head()

# container -------

  bcg_pronet <- list()

# analysis data -------

  insert_msg('Analysis data')

  bcg_pronet$lexicon <- tcga$protein_annotation %>%
    transmute(variable = variable,
              class = category,
              color = globals$protein_colors[category])

  ## caspase 3 and caspase 9 are removed because they contain only NA

  bcg_pronet$data <-
    inner_join(bcg_globals$assignment$tcga,
               tcga$protein,
               by = 'sample_id') %>%
    filter(tissue == 'tumor') %>%
    select(-tissue,
           -`CASP3|Caspase-3`, -`CASP9|Caspase-9`,
           -`ARID1A|ARID1A`, -`BRCA2|BRCA2`)%>%
    blast(clust_id, .skip = TRUE) %>%
    map(column_to_rownames, 'sample_id')

  bcg_pronet$clust_levels <- levels(bcg_pronet$data$clust_id)

# analysis globals -------

  insert_msg('Analysis globals')

  ## color scales for functional categories and communities
  ## of the proteins

  bcg_pronet$colors <-
    list(class = globals$protein_colors,
         comm_id = c('community 1' ='steelblue',
                     'community 2' = 'orangered3',
                     'community 3' = 'gold3',
                     'community 4' = 'pink2',
                     'community 5' = 'aquamarine3',
                     'community 6' = 'firebrick',
                     'community 7' = 'plum4',
                     'other' = 'gray70'))

# Similarity graphs and community detection -----------

  insert_msg('Similarity graphs')

  bcg_pronet$graph_objects <- bcg_pronet$data %>%
    map(simil_graph,
        method = 'pearson',
        corr_cutoff = 0.5,
        diag = FALSE,
        attr_lexicon = bcg_pronet$lexicon)

  ## community assignment by maximal modularity,
  ## communities with less than 10 proteins are merged

  bcg_pronet$communities <- bcg_pronet$graph_objects %>%
    map(cluster_edge_betweenness)

  bcg_pronet$comm_assignment <- bcg_pronet$communities %>%
    map(membership) %>%
    map(compress,
        names_to = 'variable',
        values_to = 'comm_id') %>%
    map(~mutate(.x,
                comm_id = factor(paste('community', comm_id)),
                comm_id = fct_lump_min(comm_id, 10, other_level = 'other'),
                index = 1:nrow(.x)))

  ## appending the vertices with their community assignment information

  for(i in names(bcg_pronet$graph_objects)) {

    V(bcg_pronet$graph_objects[[i]])$comm_id <-
      bcg_pronet$comm_assignment[[i]]$comm_id

  }

# Graph statistics ---------

  insert_msg('Graph statistics')

  ## appended with the community assignment

  bcg_pronet$stats <- bcg_pronet$graph_objects %>%
    map(get_graph_stats,
        directed = FALSE,
        normalized = TRUE) %>%
    map(arrange, -hub_score, -betweenness, -degree)

  ## top most important proteins to be labelled in the plots

  bcg_pronet$top_proteins <- bcg_pronet$stats %>%
    map(~list(slice_max(.x, hub_score, n = 15),
              slice_max(.x, betweenness, n = 15),
              slice_max(.x, degree, n = 15))) %>%
    map(map, ~.x$variable) %>%
    map(reduce, union)

  bcg_pronet$stats <-
    map2(bcg_pronet$stats,
         bcg_pronet$top_proteins,
         ~mutate(.x,
                 plot_label = ifelse(variable %in% .y,
                                     variable, NA)))

  ## plots of degree and betweenness
  ## point color codes for protein function or community assignment

  for(i in c('class', 'community')) {

    bcg_pronet$stat_plots[[i]] <-
      list(x = bcg_pronet$stats,
           y = paste0('Protein networks, cluster ',
                      names(bcg_pronet$stats),
                      ', ', globals$cohort_labs["tcga"])) %>%
      pmap(function(x, y) x %>%
             ggplot(aes(x = degree,
                        y = hub_score,
                        fill = .data[[i]],
                        size = betweenness)) +
             geom_point(color = 'black',
                        shape = 21) +
             geom_text_repel(aes(label = protein_labeller(variable)),
                             size = 2.3,
                             fontface = 'italic') +
             scale_size_area(max_size = 4,
                             limits = bcg_pronet$stats %>%
                               map(~.x$betweenness) %>%
                               reduce(c) %>%
                               range,
                             name = 'Normalized\nbetweenness') +
             globals$common_theme +
             labs(title = y,
                  subtitle = paste('proteins: n =', nrow(x)),
                  x = 'degree',
                  y = 'normalized hub score'))

  }

  ## color scales

  bcg_pronet$stat_plots$class <-
    bcg_pronet$stat_plots$class %>%
    map(~.x +
          scale_fill_manual(values = set_names(bcg_pronet$lexicon$color,
                                               bcg_pronet$lexicon$class),
                            name = 'Gene\nclassification'))

  bcg_pronet$stat_plots$community <-
    bcg_pronet$stat_plots$community %>%
    map(~.x +
          scale_fill_manual(values = bcg_pronet$colors$comm_id,
                            name = 'Protein\ncommunity'))

# Statistics for communities --------

  insert_msg('Node stats for the communities')

  bcg_pronet$community_stats <- bcg_pronet$comm_assignment %>%
    map(blast, comm_id) %>%
    map(map, ~.x$index)

  for(i in names(bcg_pronet$graph_objects)) {

    bcg_pronet$community_stats[[i]] <- bcg_pronet$community_stats[[i]] %>%
      map(~subgraph(bcg_pronet$graph_objects[[i]],
                    .x)) %>%
      map(get_graph_stats)

    bcg_pronet$top_community_proteins[[i]] <- bcg_pronet$community_stats[[i]] %>%
      map(~list(slice_max(.x, hub_score, n = 10),
                slice_max(.x, betweenness, n = 10),
                slice_max(.x, degree, n = 10))) %>%
      map(map, ~.x$variable) %>%
      map(reduce, union)

  }

# GO enrichment analysis for the 'non-other' communities ------

  insert_msg('Community GO enrichment analysis')

  bcg_pronet$community_gos <- bcg_pronet$comm_assignment %>%
    map(community_go)

# Visualization of the graphs --------

  insert_msg('Visualiztion of the graphs')

  ## with color coding of functional protein category
  ## and community assignment

  bcg_pronet$plot_proteins <-
    list(class = bcg_pronet$top_proteins,
         comm_id = bcg_pronet$top_community_proteins %>%
           map(~.x[names(.x) != 'other']) %>%
           map(reduce, union))

  for(i in c('class', 'comm_id')) {

    bcg_pronet$graph_plots[[i]] <-
      list(graph_object = bcg_pronet$graph_objects,
           plot_title = paste0('Protein networks, cluster ',
                               names(bcg_pronet$graph_objects),
                               ', ', globals$cohort_labs["tcga"]),
           select_nodes = bcg_pronet$plot_proteins[[i]]) %>%
      pmap(plot_simil_graph,
           node_color_var = i,
           txt_size = 1.8,
           node_palette = bcg_pronet$colors[[i]],
           edge_title = "Pearson's r",
           node_title = 'Protein\nclassification',
           node_txt_face = 'plain',
           cust_theme = globals$common_theme +
             theme(panel.grid.major = element_blank(),
                   axis.line = element_blank(),
                   axis.ticks = element_blank(),
                   axis.text.x = element_blank(),
                   axis.text.y = element_blank(),
                   axis.title = element_blank()),
           max.overlaps = 30,
           na.rm = TRUE,
           weighting_order = 3,
           label_fun = protein_labeller)


  }

  set.seed(1234)

# END ------

  insert_tail()
