# Plots of common (i.e.. TCGA and GSE99420) differentially regulated signaling
# pathways identified by SPIA

  insert_head()

# container ------

  bcg_spiaplots <- list()

# plotting globals -------

  insert_msg('Plotting globals')

  ## pathways to be displayed

  bcg_spiaplots$pathways <- bcg_spia$common_significant %>%
    unlist %>%
    unname %>%
    unique

  ## classification of the pathways

  bcg_spiaplots$pathway_classification <-
    c('Calcium signaling pathway' = 'class1',
      'Pathogenic Escherichia coli infection' = 'class1',

      'Focal adhesion' = 'class2',
      'ECM-receptor interaction' = 'class2',

      'Pathways in cancer' = 'class3',
      'Melanoma' = 'class3') %>%
    compress(names_to = 'Name',
             values_to = 'class')

  ## plotting data

  bcg_spiaplots$data <- bcg_spia$test %>%
    transpose %>%
    map(compress,
        names_to = 'clust_id') %>%
    map(mutate,
        clust_id = factor(clust_id,
                          levels(bcg_globals$assignment[[1]]$clust_id))) %>%
    map(filter, Name %in% bcg_spiaplots$pathways) %>%
    map(select,
        clust_id, Name, tA, regulation)

  bcg_spiaplots$data <- bcg_spiaplots$data %>%
    map(left_join,
        bcg_spiaplots$pathway_classification,
        by = 'Name')

# Bubble plots ------

  insert_msg('Bubble plots')

  bcg_spiaplots$bubble_plots <-
    list(x = bcg_spiaplots$data,
         y = paste('SPIA signaling pathways,',
                   globals$cohort_labs[names(bcg_spiaplots$data)])) %>%
    pmap(function(x, y) x %>%
           ggplot(aes(x = clust_id,
                      y = Name,
                      fill = regulation,
                      color = regulation,
                      size = abs(tA))) +
           facet_grid(class ~ .,
                      space = 'free',
                      scales = 'free') +
           geom_point(shape = 21,
                      color = 'black') +
           geom_text(aes(label = signif(tA, 2)),
                     size = 2.5,
                     hjust = 0.5,
                     vjust = -1.2,
                     show.legend = FALSE) +
           scale_fill_manual(values = c(activated = 'firebrick',
                                        inhibited = 'steelblue',
                                        ns = 'gray70'),
                             name = 'Regulation\nvs cohort mean') +
           scale_color_manual(values = c(activated = 'firebrick',
                                         inhibited = 'steelblue',
                                         ns = 'gray70'),
                              name = 'Regulation\nvs cohort mean') +
           scale_size_area(max_size = 5,
                           limits = bcg_spiaplots$data %>%
                             map(~.x$tA) %>%
                             reduce(c) %>%
                             abs %>%
                             range,
                           name = 'abs(tA)') +
           globals$common_theme +
           theme(axis.title.y = element_blank(),
                 strip.background = element_blank(),
                 strip.text = element_blank()) +
           labs(title = y,
                subtitle = 'Common regulated pathways vs cohort mean',
                x = 'cluster'))

# END ------

  bcg_spiaplots <-
    bcg_spiaplots[c("pathways", "pathway_classification", "bubble_plots")]

  insert_tail()
