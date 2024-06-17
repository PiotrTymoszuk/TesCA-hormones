# Plots for the drug sensitivity analyses

  insert_head()

# container -----

  bcg_drugplots <- list()

# analysis globals ------

  insert_msg('Analysis globals')

  ## total numbers of analyzed drugs

  bcg_drugplots$n_total <- bcg_drugs$n_numbers %>%
    blast(experiment) %>%
    map_dbl(~.x$n_total[1])

  ## cluster levels

  bcg_drugplots$clust_levels <- levels(bcg_globals$assignment[[1]]$clust_id)

  ## common differentially regulated drugs and their predictions:
  ## to be displayed in heat maps

  bcg_drugplots$plot_variables <- bcg_drugs$common_significant %>%
    map(unlist) %>%
    map(unname) %>%
    map(unique)

  bcg_drugplots$predictions <-
    map2(drugs$predictions,
         bcg_drugplots$plot_variables,
         function(data, variable) data %>%
           map(select, sample_id, all_of(variable))) %>%
    map(function(data) map2(bcg_globals$assignment,
                            data,
                            inner_join, by = 'sample_id'))

# Numbers of significant drugs -------

  insert_msg('Numbers of significant drugs')

  bcg_drugplots$n_numbers <-
    list(x = blast(bcg_drugs$n_numbers,
                   experiment),
         y = paste0(globals$drug_exp_labs,
                    '-trained predictions'),
         z = bcg_drugplots$n_total) %>%
    pmap(function(x, y, z) x %>%
           mutate(plot_percent = ifelse(status == 'sensitive',
                                        -percent, percent)) %>%
           ggplot(aes(x = plot_percent,
                      y = cohort,
                      fill = status)) +
           facet_grid(. ~ clust_id) +
           geom_vline(xintercept = 0,
                      linetype = 'dashed') +
           geom_bar(stat = 'identity',
                    color = 'black') +
           scale_fill_manual(values = c(resistant = 'firebrick',
                                        sensitive = 'steelblue'),
                             name = 'Predicted response\nvs cohort mean') +
           scale_y_discrete(labels = globals$cohort_labs) +
           scale_x_continuous(labels = function(x) abs(x)) +
           globals$common_theme +
           theme(axis.title.y = element_blank()) +
           labs(title = y,
                subtitle = paste('total drugs: n =', z),
                x = '% of investigated drugs'))

# Volcano plots -------

  insert_msg('Volcano plots')

  for(i in names(bcg_drugs$dev_test)) {

    for(j in names(bcg_drugs$dev_test[[i]])) {

      bcg_drugplots$volcano_plots[[i]][[j]]<-
        list(data = blast(bcg_drugs$dev_test[[i]][[j]],
                          clust_id),
             plot_title = paste0(globals$drug_exp_labs[i],
                                 ' predictions, cluster ',
                                 bcg_drugplots$clust_levels,
                                 ', ', globals$cohort_labs[j])) %>%
        pmap(plot_volcano,
             regulation_variable = 'deviation_center',
             p_variable = 'p_adjusted',
             regulation_level = 0,
             x_lab = paste('\u0394',
                           globals$drug_unit_labs[i],
                           'vs cohort mean'),
             y_lab = expression('-log'[10] * ' pFDR'),
             top_regulated = 10,
             label_variable = 'drug_name',
             label_type = 'text',
             txt_size = 2.3,
             cust_theme = globals$common_theme) %>%
        map(~.x +
              labs(subtitle = .x$labels$tag) +
              theme(plot.tag = element_blank()) +
              scale_fill_manual(values = c(upregulated = 'firebrick',
                                           downregulated = 'steelblue',
                                           ns = 'gray70'),
                                labels = c(upregulated = 'resistant',
                                           downregulated = 'sensitive',
                                           ns = 'ns'),
                                name = 'Predicted response\nvs cohort mean'))

    }

  }

# Forest plots for the top regulated drugs -------

  insert_msg('Forest plots for the top regulated drugs')

  for(i in names(bcg_drugs$dev_test)) {

    for(j in names(bcg_drugs$dev_test[[i]])) {

      bcg_drugplots$top_plots[[i]][[j]]<-
        list(data = bcg_drugs$dev_test[[i]][[j]] %>%
               filter(regulation %in% c('resistant', 'sensitive')) %>%
               blast(clust_id),
             plot_title = paste0(globals$drug_exp_labs[i],
                                 ' predictions, cluster ',
                                 bcg_drugplots$clust_levels,
                                 ', ', globals$cohort_labs[j])) %>%
        pmap(plot_top,
             regulation_variable = 'deviation_center',
             p_variable = 'p_adjusted',
             label_variable = 'drug_name',
             regulation_level = 0,
             top_regulated = 20,
             lower_ci_variable = 'lower_ci',
             upper_ci_variable = 'upper_ci',
             x_lab = paste('\u0394',
                           globals$drug_unit_labs[i],
                           'vs cohort mean'),
             plot_subtitle = 'Top differences in predicted drug response',
             cust_theme = globals$common_theme) %>%
        map(~.x +
              scale_fill_manual(values = c(upregulated = 'firebrick',
                                           downregulated = 'steelblue',
                                           ns = 'gray70'),
                                labels = c(upregulated = 'resistant',
                                           downregulated = 'sensitive',
                                           ns = 'ns'),
                                name = 'Predicted response\nvs cohort mean') +
              scale_color_manual(values = c(upregulated = 'firebrick',
                                           downregulated = 'steelblue',
                                           ns = 'gray70'),
                                labels = c(upregulated = 'resistant',
                                           downregulated = 'sensitive',
                                           ns = 'ns'),
                                name = 'Predicted response\nvs cohort mean'))

    }

  }

# Classification of the common regulated drugs ------

  insert_msg('Classification of the common regulated drugs')

  ## done for the TCGA cohort and applied to the remaining collectives

  bcg_drugplots$classification_tcga <-
    list(data = map(bcg_drugplots$predictions, ~.x$tcga),
         variables = bcg_drugplots$plot_variables) %>%
    pmap(classify,
         split_fct = 'clust_id') %>%
    map(~.x$classification)

# Heat maps --------

  insert_msg('Heat maps')

  for(i in names(bcg_drugplots$predictions)) {

    bcg_drugplots$hm_plots[[i]] <-
      list(data = bcg_drugplots$predictions[[i]],
           plot_title = paste(globals$drug_exp_labs[i],
                               globals$cohort_labs[names(bcg_drugplots$predictions[[i]])],
                              sep = ', ')) %>%
      pmap(heat_map,
           variables = bcg_drugplots$plot_variables[[i]],
           split_fct = 'clust_id',
           normalize = TRUE,
           variable_classification = bcg_drugplots$classification_tcga[[i]],
           cust_theme = globals$common_theme +
             theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.y = element_blank()),
           x_lab = 'cancer sample',
           midpoint = 0,
           limits = c(-3, 3),
           oob = scales::squish)

  }

  ## Y axis styling

  bcg_drugplots$hm_plots$ctrp2 <-
    bcg_drugplots$hm_plots$ctrp2 %>%
    map(~.x +
          guides(y = guide_axis(n.dodge = 3)) +
          scale_y_discrete(labels = function(x) exchange(x,
                                                         drugs$lexicons$ctrp2,
                                                         value = 'drug_name')))

  bcg_drugplots$hm_plots$gdsc <-
    bcg_drugplots$hm_plots$gdsc %>%
    map(~.x +
          guides(y = guide_axis(n.dodge = 3)) +
          scale_y_discrete(labels = function(x) exchange(x,
                                                         drugs$lexicons$gdsc,
                                                         value = 'drug_name')))

# END ------

  bcg_drugplots <-
    bcg_drugplots[c("n_numbers",
                    "volcano_plots", "top_plots",
                    "classification_tcga", "hm_plots")]

  rm(i, j)

  insert_tail()
