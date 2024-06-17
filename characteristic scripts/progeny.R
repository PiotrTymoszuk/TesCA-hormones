# Differential regulation of PROGENy signaling pathways

  insert_head()

# container ------

  bcg_progeny <- list()

  # parallel backend ------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis globals ------

  insert_msg('Analysis globals')

  ## the collecTRI model

  if(file.exists('./data/progeny.RData')) {

    load('./data/progeny.RData')

  } else {

    progeny <- get_progeny()

    save(progeny, file = './data/progeny.RData')

  }

  ## t stats of the differential gene expression
  ## in the cluster as compared with the cohort mean

  bcg_progeny$t_stats <- bcg_dge$dev_test %>%
    map(filter,
        !is.na(gene_symbol),
        !is.na(t),
        !is.infinite(t)) %>%
    map(blast, clust_id) %>%
    map(map, select, gene_symbol, t) %>%
    map(map, column_to_rownames, 'gene_symbol') %>%
    map(map, as.matrix)

# Modeling -------

  insert_msg('Modeling')

  for(i in names(bcg_progeny$t_stats)) {

    bcg_progeny$test[[i]] <- bcg_progeny$t_stats[[i]] %>%
      future_map(run_mlm,
                 network = progeny,
                 .source = 'source',
                 .target = 'target',
                 .mor = 'weight',
                 minsize = 5,
                 .options = furrr_options(seed = TRUE))

    ## formatting the results

    bcg_progeny$test[[i]] <- bcg_progeny$test[[i]] %>%
      map(re_adjust, method = 'BH') %>%
      map(mutate,
          regulation = ifelse(p_adjusted >= 0.05, 'ns',
                              ifelse(score > 0, 'activated',
                                     ifelse(score < 0,
                                            'inhibited', 'ns'))),
          regulation = factor(regulation,
                              c('activated', 'inhibited', 'ns')))

  }

  bcg_progeny$test <- transpose(bcg_progeny$test)

# significantly modulated signaling pathways ------

  insert_msg('Significant pathways')

  ## in single cohorts

  bcg_progeny$significant <- bcg_progeny$test %>%
    map(map, filter, regulation %in% c('activated', 'inhibited')) %>%
    map(map, blast, regulation) %>%
    map(transpose) %>%
    map(map, map, ~.x$source)

  ## common ones: shared by the TCGA and GSE99420 cohort

  bcg_progeny$common_significant <-
    bcg_progeny$significant %>%
    map(map, ~.x[c('tcga', 'gse99420')]) %>%
    map(map, reduce, intersect)

# Bubble plots for the common significant signaling pathways --------

  insert_msg('Bubble plots for the common significant pathways')

  ## plotting data, classification of the pathways based
  ## on the pathway functions

  bcg_progeny$plot_variables <- bcg_progeny$common_significant %>%
    unlist %>%
    unname %>%
    unique

  bcg_progeny$plot_data <- bcg_progeny$test %>%
    transpose %>%
    map(compress, names_to = 'clust_id') %>%
    map(mutate,
        clust_id = factor(clust_id,
                             levels(bcg_globals$assignment[[1]]$clust_id))) %>%
    map(filter,
        source %in% bcg_progeny$plot_variables) %>%
    map(select, clust_id, source, score, regulation)

  bcg_progeny$plot_classification <-
    c('Hypoxia' = 'class1',
      'VEGF' = 'class1',

      'WNT' = 'class2',
      'TGFb' = 'class2',
      'EGFR' = 'class2',
      'MAPK' = 'class2',

      'JAK-STAT' = 'class3',
      'p53' = 'class3',

      'Estrogen' = 'class4',
      'Androgen' = 'class4') %>%
    compress(names_to = 'source',
             values_to = 'class')

  bcg_progeny$plot_data <- bcg_progeny$plot_data %>%
    map(left_join,
        bcg_progeny$plot_classification,
        by = 'source') %>%
    map(mutate,
        source = factor(source,
                        rev(bcg_progeny$plot_classification$source)))

  ## bubble plots

  bcg_progeny$bubble_plots <-
    list(x = bcg_progeny$plot_data,
         y = paste('PROGENy signaling pathways,',
                   globals$cohort_labs[names(bcg_progeny$plot_data)])) %>%
    pmap(function(x, y) x %>%
           ggplot(aes(x = clust_id,
                      y = source,
                      fill = regulation,
                      color = regulation,
                      size = abs(score))) +
           facet_grid(class ~ .,
                      scales = 'free',
                      space = 'free') +
           geom_point(shape = 21,
                      color = 'black') +
           geom_text(aes(label = signif(score, 2)),
                     size = 2.5,
                     hjust = 0.5,
                     vjust = -1.2,
                     show.legend = FALSE) +
           scale_fill_manual(values = c(activated = 'firebrick',
                                        inhibited = 'steelblue',
                                        ns = 'gray70'),
                             name = 'Regulation\nvc cohort mean') +
           scale_color_manual(values = c(activated = 'firebrick',
                                         inhibited = 'steelblue',
                                         ns = 'gray70'),
                              name = 'Regulation\nvc cohort mean') +
           scale_size_area(max_size = 5,
                           limits = bcg_progeny$plot_data %>%
                             map(~.x$score) %>%
                             reduce(c) %>%
                             abs %>%
                             range,
                           name = 'abs(LM score)') +
           globals$common_theme +
           theme(axis.title.y = element_blank(),
                 strip.background = element_blank(),
                 strip.text = element_blank()) +
           labs(title = y,
                subtitle = 'Common regulated pathways vs cohort mean',
                x = 'cluster'))

# END -----

  bcg_progeny <- bcg_progeny[c("test",
                               "significant", "common_significant",
                               "bubble_plots")]

  rm(i)

  plan('sequential')

  insert_tail()
