# Response to approved drugs in the hormonal clusters

  insert_head()

# container -------

  bcg_approved <- list()

# drugs of interest, their predictions and analysis results ----

  insert_msg('Drugs of interest')

  ## drug search regular expressions

  bcg_approved$regex$plain <-
    c('platin', 'bleomycin', 'etoposide',
      'vinblastine', 'vincristine', 'ifosfamide',
      'methotrexate', 'actinomycin',
      'cyclophosphamide', 'gemcitabine')

  bcg_approved$regex$capital <- bcg_approved$regex$plain %>%
    stri_capitalize

  bcg_approved$regex$upper <- bcg_approved$regex$plain %>%
    toupper

  bcg_approved$regex <- bcg_approved$regex %>%
    map_chr(paste, collapse = '|') %>%
    paste(collapse = '|')

  ## ANOVA results

  bcg_approved$anova <- bcg_drugs$anova %>%
    map(map, reglook, bcg_approved$regex) %>%
    map(map, re_adjust, method = 'none')

  bcg_approved$anova <- bcg_approved$anova %>%
    map(map,
        mutate,
        eff_size = paste('\u03B7\u00B2 =', signif(effect_size, 2)),
        plot_cap = paste(eff_size, significance, sep = ', '),
        axis_lab = paste(drug_name, plot_cap, sep = '<br>'),
        axis_lab = ifelse(p_adjusted < 0.05 & effect_size >= 0.14,
                          html_bold(axis_lab), axis_lab)) %>%
    map(map,
        select,
        variable,
        drug_name,
        plot_cap,
        axis_lab)

  ## variable names

  bcg_approved$variables <- bcg_approved$anova %>%
    map(~.x[[1]]$variable)

  ## sensitivity estimates

  bcg_approved$data <-
    map2(drugs$predictions,
         bcg_approved$variables,
         function(data, variable) data %>%
           map(select, sample_id, all_of(variable))) %>%
    map(function(data) map2(bcg_globals$assignment,
                            data,
                            inner_join,
                            by = 'sample_id'))

# Box plots for single drugs -----

  insert_msg('Box plots')

  for(i in names(bcg_approved$data)) {

    for(j in names(bcg_approved$data[[i]])) {

      bcg_approved$plots[[i]][[j]] <-
        list(variable = bcg_approved$anova[[i]][[j]]$variable,
             plot_title = bcg_approved$anova[[i]][[j]]$drug_name %>%
               paste(globals$cohort_labs[j],
                     sep = ', '),
             plot_subtitle = bcg_approved$anova[[i]][[j]]$plot_cap) %>%
        pmap(plot_variable,
             bcg_approved$data[[i]][[j]],
             split_factor = 'clust_id',
             type = 'box',
             cust_theme = globals$common_theme,
             x_n_labs = TRUE,
             y_lab = paste(globals$drug_exp_labs[i],
                           globals$drug_unit_labs[i],
                           sep = '-trained predictions, '),
             x_lab = 'cluster') %>%
        map(~.x +
              scale_fill_manual(values =  globals$cluster_colors)) %>%
        set_names(bcg_approved$anova[[i]][[j]]$variable)

    }

  }

# Classification of the variables by cluster specificity, TCGA ------

  insert_msg('Classification of the variables')

  ## classification in the TCGA cohort applied to all other cohorts

  bcg_approved$classification_tcga <-
    list(data = map(bcg_approved$data, ~.x$tcga),
         variables = bcg_approved$variables) %>%
    pmap(classify,
         split_fct = 'clust_id') %>%
    map(~.x$classification)

# Heat map plots --------

  insert_msg('Heat map plots')

  for(i in names(bcg_approved$data)) {

    bcg_approved$heat_maps[[i]] <-
      list(data = bcg_approved$data[[i]],
           plot_title = globals$drug_exp_labs[i] %>%
             paste(globals$cohort_labs[names(bcg_approved$data[[i]])],
                   sep = '-trained predictions, ')) %>%
      pmap(heat_map,
           split_fct = 'clust_id',
           variables = bcg_approved$variables[[i]],
           normalize = TRUE,
           variable_classification = bcg_approved$classification_tcga[[i]],
           cust_theme = globals$common_theme +
             theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.y = element_blank(),
                   axis.text.y = element_markdown()),
           midpoint = 0,
           limits = c(-3, 3),
           oob = scales::squish,
           x_lab = 'cancer sample')

    ## additional styling

    bcg_approved$heat_maps[[i]] <-
      list(x = bcg_approved$heat_maps[[i]],
           y = bcg_approved$anova[[i]]) %>%
      pmap(function(x, y) x +
             scale_y_discrete(labels = set_names(y$axis_lab,
                                                 y$variable)))


  }

# END -----

  rm(i)

  bcg_approved$regex <- NULL
  bcg_approved$anova <- NULL
  bcg_approved$data <- NULL

  bcg_approved <- compact(bcg_approved)

  insert_tail()
