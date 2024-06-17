# differential regulation of transcriptional regulons between the
# hormonal clusters done with the decoupler's univariable linear modeling tool.
# Significant modulation of a regulon is assumed for pFDR < 0.05
# and abs(LM score) > 0

  insert_head()

# container -------

  bcg_collectri <- list()

# parallel backend ------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis globals ------

  insert_msg('Analysis globals')

  ## the collecTRI model

  if(file.exists('./data/collectri.RData')) {

    load('./data/collectri.RData')

  } else {

    collectri <- get_collectri(split_complexes = FALSE)

    save(collectri, file = './data/collectri.RData')

  }

  ## t stats of the differential gene expression
  ## in the cluster as compared with the cohort mean

  bcg_collectri$t_stats <- bcg_dge$dev_test %>%
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

  for(i in names(bcg_collectri$t_stats)) {

    bcg_collectri$test[[i]] <- bcg_collectri$t_stats[[i]] %>%
      future_map(run_ulm,
                 network = collectri,
                 .source = 'source',
                 .target = 'target',
                 .mor = 'mor',
                 minsize = 5,
                 .options = furrr_options(seed = TRUE))

    ## formatting the results

    bcg_collectri$test[[i]] <- bcg_collectri$test[[i]] %>%
      map(re_adjust, method = 'BH') %>%
      map(mutate,
          regulation = ifelse(p_adjusted >= 0.05, 'ns',
                              ifelse(score > 0, 'activated',
                                     ifelse(score < 0,
                                            'inhibited', 'ns'))),
          regulation = factor(regulation,
                              c('activated', 'inhibited', 'ns')))

  }

  bcg_collectri$test <- transpose(bcg_collectri$test)

# significant regulons ------

  insert_msg('Significant regulons')

  ## in single cohorts

  bcg_collectri$significant <- bcg_collectri$test %>%
    map(map, filter, regulation %in% c('activated', 'inhibited')) %>%
    map(map, blast, regulation) %>%
    map(transpose) %>%
    map(map, map, ~.x$source)

  ## common ones: shared by the TCGA and GSE99420 cohort

  bcg_collectri$common_significant <-
    bcg_collectri$significant %>%
    map(map, ~.x[c('tcga', 'gse99420')]) %>%
    map(map, reduce, intersect)

  ## numbers of the common regulons and ready-to-use plot subtitles

  bcg_collectri$common_n <- bcg_collectri$common_significant %>%
    map(map_dbl, length) %>%
    map(~paste0('activated: n = ', .x[1],
                ', inhibited: n = ', .x[2]))

# Bar plots for the common significant regulons -------

  insert_msg('Bubble plots')

  ## variables to be presented in the plots
  ## and data: TCGA and GSE99420 cohorts

  bcg_collectri$plot_variables <-
    bcg_collectri$common_significant %>%
    map(reduce, union)

  bcg_collectri$plot_data <- bcg_collectri$test %>%
    map(compress, names_to = 'cohort') %>%
    map(filter, cohort %in% c('tcga', 'gse99420')) %>%
    map(mutate,
        cohort = factor(cohort, globals$analysis_cohorts),
        cohort = droplevels(cohort)) %>%
    map(select,
        cohort, source, score, regulation)

    bcg_collectri$plot_data <-
      map2(bcg_collectri$plot_data,
           bcg_collectri$plot_variables,
           ~filter(.x, source %in% .y))

    ## bar plots

    bcg_collectri$bar_plots <-
      list(x = bcg_collectri$plot_data ,
           y = paste('Regulons, cluster',
                     names(bcg_collectri$plot_data)),
           z =  bcg_collectri$common_n,
           w = c(1, 2, 2, 1)) %>%
      pmap(function(x, y, z, w) x %>%
             ggplot(aes(x = score,
                        y = reorder(source, score),
                        fill = regulation)) +
             facet_grid(. ~ cohort,
                        labeller = as_labeller(globals$cohort_labs)) +
             geom_vline(xintercept = 0,
                        linetype = 'dashed') +
             geom_bar(stat = 'identity') +
             scale_fill_manual(values = c(activated = 'firebrick',
                                          inhibited = 'steelblue'),
                               name = 'activity status\vs cohort mean') +
             guides(y = guide_axis(n.dodge = w)) +
             globals$common_theme +
             theme(axis.title.y = element_blank()) +
             labs(title = y,
                  subtitle = z,
                  x = 'LM score'))

# Bar plots with the top modulated regulons ------

  insert_msg('Bar plots for the top regulons')

  ## the plots are generated only for the TCGA and GSE99420 cohort

  ## plotting data: top 10 up- and downregulated
  ## transcriptional regulons per cluster shared by both cohorts

  bcg_collectri$top_data <-
    map2(bcg_collectri$test,
         bcg_collectri$plot_variables,
         function(test, vars) test %>%
           map(filter, source %in% vars)) %>%
    map(~.x[c('tcga', 'gse99420')]) %>%
    map(map, group_by, sign(score)) %>%
    map(map, top_n, n = 10, abs(score)) %>%
    map(map, ungroup)

  bcg_collectri$top_variables <- bcg_collectri$top_data %>%
    map(map, ~.x$source) %>%
    map(reduce, intersect)

  bcg_collectri$top_data <- bcg_collectri$top_data %>%
    map(compress, names_to = 'cohort') %>%
    map2(., bcg_collectri$top_variables,
         ~filter(.x, source %in% .y)) %>%
    compress(names_to = 'clust_id') %>%
    mutate(cohort = factor(cohort, c('tcga', 'gse99420')),
           clust_id = factor(clust_id,
                             levels(bcg_globals$assignment[[1]]$clust_id)))

  ## bar plot with LM scores

  bcg_collectri$top_plot <- bcg_collectri$top_data %>%
    ggplot(aes(x = score,
               y = reorder(source, score),
               fill = regulation)) +
    facet_grid(clust_id ~ cohort,
               scales = 'free',
               space = 'free_y',
               labeller = labeller(.cols = globals$cohort_labs)) +
    geom_vline(xintercept = 0,
               linetype = 'dashed') +
    geom_bar(stat = 'identity',
             color = 'black') +
    scale_fill_manual(values = c(activated = 'firebrick',
                                 inhibited = 'steelblue'),
                      name = 'Regulation\nvs cohort mean') +
    globals$common_theme +
    theme(axis.title.y = element_blank()) +
    labs(title = 'Top collecTRI regulons',
         subtitle = 'Common modulated regulons vs cohort mean',
         x = 'LM score')

# END ------

  bcg_collectri <-
      bcg_collectri[c("test",
                      "significant", "common_significant",
                      "common_n", "bar_plots", "top_plot")]

  rm(i)

  plan('sequential')

  insert_tail()
