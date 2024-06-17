# Differences in non-malignant cell infiltration between the hormonal clusters.
#
# Statistical significance is determined by Kruskal-Wallis test with eta-square
# effect size statistic. Biologically relevant differences are considered for
# pFDR < 0.05 and eta-square >= 0.14.

  insert_head()

# container -----

  bcg_infil <- list()

# parallel backend -----

  insert_msg('Parallel backend')

  plan('multisession')

# analysis globals ------

  insert_msg('Analysis globals')

  ## analysis data, appending with the cluster assignment

  bcg_infil$data <- list(quantiseq = quantiseq,
                         xcell = xcell,
                         mcp = mcp) %>%
    map(~.x$infiltration)

  for(i in names(bcg_infil$data)) {

    bcg_infil$data[[i]] <-
      map2(bcg_globals$assignment[names(bcg_infil$data[[i]])],
           bcg_infil$data[[i]],
           inner_join, by = 'sample_id') %>%
      map(~filter(.x, complete.cases(.x)))

  }

  ## variables, scores and progenitor cells are removed

  bcg_infil$lexicons <- list(quantiseq = quantiseq,
                             xcell = xcell,
                             mcp = mcp) %>%
    map(~.x$lexicon) %>%
    map(filter,
        !stri_detect(variable, regex = 'progenitor$'),
        !stri_detect(variable, regex = 'score$'),
        !stri_detect(variable, fixed = 'stem'))

  ## labels of the deconvolution algorithms

  bcg_infil$algo_labs <-
    c(quantiseq = 'QuanTIseq',
         xcell = 'xCell',
         mcp = 'MCP')

  bcg_infil$algo_scale_labs <-
    c(quantiseq = 'fraction of tumor cells, QuanTIseq',
      xcell = 'fraction of tumor cells, xCell',
      mcp = 'cell count, MCP Counter')

# Descriptive stats -------

  insert_msg('Descriptive stats')

  for(i in names(bcg_infil$data)) {

    bcg_infil$stats[[i]] <- bcg_infil$data[[i]] %>%
      map(select, clust_id, all_of(bcg_infil$lexicon[[i]]$variable)) %>%
      map(fast_num_stats,
          split_factor = 'clust_id')

  }

# Testing for differences between the clusters ------

  insert_msg('Testing for differences between the clusters')

  for(i in names(bcg_infil$data)) {

    bcg_infil$test[[i]] <- bcg_infil$data[[i]] %>%
      future_map(compare_variables,
                 variables = bcg_infil$lexicons[[i]]$variable,
                 split_factor = 'clust_id',
                 what = 'eff_size',
                 types = 'kruskal_etasq',
                 ci = FALSE,
                 exact = FALSE,
                 pub_styled = FALSE,
                 adj_method = 'BH',
                 .options = furrr_options(seed = TRUE)) %>%
      map(mutate,
          eff_size = paste('\u03B7\u00B2 =', signif(estimate, 2)),
          plot_cap = paste(eff_size, significance, sep = ', '),
          ## axis labels to be presented in heat maps,
          ## significant effects highlighted in bold
          axis_lab = exchange(variable,
                              bcg_infil$lexicons[[i]]),
          axis_lab = ifelse(p_adjusted < 0.05,
                            html_bold(axis_lab), axis_lab))

    ## significant effects,
    ## common significant effects shared by the TCGA and GSE99420

    bcg_infil$significant[[i]] <- bcg_infil$test[[i]] %>%
      map(filter,
          p_adjusted < 0.05) %>%
      map(~.x$variable)

    bcg_infil$common_significant[[i]] <-
      bcg_infil$significant[[i]][c('tcga', 'gse99420')] %>%
      reduce(intersect)

  }

  ## top significant variables: significant effects for at least two
  ## algorithms

  bcg_infil$top_significant <- bcg_infil$common_significant %>%
    shared_features(m = 2) %>%
    as.character

# Box plots of single variables --------

  insert_msg('Box plots for single variables')

  ## N numbers of samples in the clusters
  ## to be presented in legends of the box plots and used later
  ## in box plot panels

  for(i in names(bcg_infil$data)) {

    bcg_infil$legend_labs[[i]] <- bcg_infil$data[[i]] %>%
      map(count, clust_id) %>%
      map(~map2_chr(.x[[1]], .x[[2]],
                    paste, sep = ': n = '))

  }

  ## box plots

  for(i in names(bcg_infil$data)) {

    for(j in names(bcg_infil$data[[i]])) {

      bcg_infil$plots[[i]][[j]] <-
        list(variable = bcg_infil$test[[i]][[j]]$variable,
             plot_title = bcg_infil$test[[i]][[j]]$variable %>%
               exchange(bcg_infil$lexicons[[i]]) %>%
               paste(globals$cohort_labs[j], sep = ', '),
             plot_subtitle = bcg_infil$test[[i]][[j]]$plot_cap) %>%
        future_pmap(plot_variable,
                    bcg_infil$data[[i]][[j]],
                    split_factor = 'clust_id',
                    type = 'box',
                    cust_theme = globals$common_theme,
                    y_lab = bcg_infil$algo_scale_labs[i],
                    x_lab = 'cluster',
                    x_n_labs = TRUE,
                    point_hjitter = 0,
                    .options = furrr_options(seed = TRUE)) %>%
        map(~.x +
              scale_fill_manual(values = globals$cluster_colors,
                                labels = bcg_infil$legend_labs[[i]][[j]])) %>%
        set_names(bcg_infil$test[[i]][[j]]$variable)


    }

    bcg_infil$plots[[i]] <- transpose(bcg_infil$plots[[i]])

  }

# Result tables -------

  insert_msg('Result tables')

  for(i in names(bcg_infil$data)) {

    bcg_infil$result_tbl[[i]] <-
      map2(bcg_infil$stats[[i]],
           map(bcg_infil$test[[i]],
               ~.x[c('variable', 'significance', 'eff_size')]),
           left_join, by = 'variable') %>%
      compress(names_to = 'cohort') %>%
      mutate(cohort = globals$cohort_labs[cohort]) %>%
      relocate(cohort) %>%
      set_names(c('Cohort',
                  'Variable',
                  levels(bcg_infil$data[[i]][[1]]$clust_id),
                  'Significance',
                  'Effect size'))

  }

# Classification of the common significant variables by the cluster specificity ------

  insert_msg('Classification of the common significant variables')

  ## in the TCGA cohort

  bcg_infil$classification_tcga <- bcg_infil$data %>%
    map(~.x$tcga) %>%
    list(data = .,
         variables = bcg_infil$common_significant) %>%
    pmap(classify,
         split_fct = 'clust_id') %>%
    map(~.x$classification)

# Heat maps for the common significant infiltration estimates -------

  insert_msg('Heat maps for the common significant variables')

  for(i in names(bcg_infil$data)) {

    bcg_infil$hm_plots[[i]] <-
      list(data = bcg_infil$data[[i]],
           plot_title = bcg_infil$algo_labs[i] %>%
             paste(globals$cohort_labs[names(bcg_infil$data[[i]])],
                   sep = ', ')) %>%
      pmap(heat_map,
           variables = bcg_infil$common_significant[[i]],
           split_fct = 'clust_id',
           normalize = TRUE,
           variable_classification = bcg_infil$classification_tcga[[i]],
           cust_theme = globals$common_theme +
             theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.y = element_blank(),
                   axis.text.y = element_markdown()),
           limits = c(-3, 3),
           oob = scales::squish,
           midpoint = 0,
           x_lab = 'cancer sample')

    ## setting the Y axis scales

    bcg_infil$hm_plots[[i]] <-
      list(x = bcg_infil$hm_plots[[i]],
           y = bcg_infil$test[[i]]) %>%
      pmap(function(x, y) x +
             scale_y_discrete(labels = set_names(y$axis_lab,
                                                 y$variable)))

  }

# Ribbon plots for the common significant effects ---------

  insert_msg('Ribbon plots for the common significant effects')

  for(i in names(bcg_infil$data)) {

    bcg_infil$ribbon_data[[i]] <- bcg_infil$data[[i]]

    for(j in names(bcg_infil$ribbon_data[[i]])) {

      bcg_infil$ribbon_data[[i]][[j]][bcg_infil$common_significant[[i]]] <-
        bcg_infil$ribbon_data[[i]][[j]][bcg_infil$common_significant[[i]]] %>%
        center_data

    }

    bcg_infil$ribbon_plots[[i]] <-
      list(data = bcg_infil$ribbon_data[[i]],
           plot_title = bcg_infil$algo_labs[i] %>%
             paste(globals$cohort_labs[names(bcg_infil$data[[i]])],
                   sep = ', ')) %>%
      pmap(draw_stat_panel,
           variables = bcg_infil$common_significant[[i]],
           split_factor = 'clust_id',
           stat = 'mean',
           err_stat = '2se',
           form = 'line',
           x_lab = 'Z-score, mean \u00B1 2\u00D7SEM',
           cust_theme = globals$common_theme +
             theme(axis.title.y = element_blank(),
                   axis.text.y = element_markdown()))

    ## setting the Y axis feature order: specificity for the cluster
    ## significant effects are labeled in bold

    bcg_infil$ribbon_plots[[i]] <-
      list(x = bcg_infil$ribbon_plots[[i]],
           y = bcg_infil$test[[i]]) %>%
      pmap(function(x, y) x +
             scale_fill_manual(values = globals$cluster_colors,
                               name = 'cluster') +
             scale_color_manual(values = globals$cluster_colors,
                                name = 'cluster') +
             scale_y_discrete(limits = bcg_infil$classification_tcga[[i]]$variable,
                              labels = set_names(y$axis_lab,
                                                 y$variable)))

  }

# END ------

  bcg_infil$data <- NULL
  bcg_infil$lexicons <- NULL
  bcg_infil$n_numbers <- NULL
  bcg_infil$ribbon_data <- NULL

  bcg_infil <- compact(bcg_infil)

  rm(i, j)

  insert_tail()
