# Plots for numbers of significantly regulated metabolic reactions,
# and regulation of selected metabolic subsystems

  insert_head()

# container ------

  bcg_metaplots <- list()

# Analysis globals -------

  insert_msg('Analysis globals')

  ## subsystems of interest and log2-fold estimates
  ## of their reaction regulation.
  ## we're filtering the pyrophosphatase reaction out, which is not
  ## essential a mamber of OxPhos subsystem

  bcg_metaplots$subsystems <-
    c('Citric acid cycle',
      'Oxidative phosphorylation',
      'Fatty acid oxidation',
      'Steroid metabolism',
      'Androgen and estrogen synthesis and metabolism')

  bcg_metaplots$estimates <- bcg_meta$estimates %>%
    transpose %>%
    map(map,
        filter,
        subsystem %in% bcg_metaplots$subsystems,
        !react_id %in% c('R_r0205', 'R_PPA')) %>%
    map(compress, names_to = 'clust_id') %>%
    map(transmute,
        clust_id = factor(clust_id,
                          levels(bcg_globals$assignment[[1]]$clust_id)),
        subsystem = factor(subsystem,
                           bcg_metaplots$subsystems),
        react_id = react_id,
        react_name = annotate_bigg(react_id),
        react_name = ifelse(is.na(react_name),
                            react_id,
                            stri_replace(react_name,
                                         regex = '^R_',
                                         replacement = '')),
        regulation = regulation,
        fold_reg = log2(fold_reg)) %>%
    map(arrange, subsystem, regulation, fold_reg)

  bcg_metaplots$reactions <- bcg_metaplots$estimates %>%
    map(~.x$react_id) %>%
    reduce(union)

  ## reactions significantly regulated both in the TCGA and GSE99420 cohorts

  bcg_metaplots$common_significant <- bcg_meta$common_significant %>%
    map(reduce, union) %>%
    map(~intersect(.x, bcg_metaplots$reactions))

# Numbers of significantly regulated reactions ------

  insert_msg('Reaction numbers')

  bcg_metaplots$n_numbers <- bcg_meta$react_numbers %>%
    mutate(percent = n/n_total * 100,
           plot_percent = ifelse(status == 'inhibited',
                                 -percent, percent),
           cohort = factor(cohort,
                           rev(globals$analysis_cohorts)),
           clust_id = factor(clust_id,
                             levels(bcg_globals$assignment[[1]]$clust_id))) %>%
    ggplot(aes(x = plot_percent,
               y = cohort,
               fill = status)) +
    geom_vline(xintercept = 0,
               linetype = 'dashed') +
    geom_bar(stat = 'identity',
             color = 'black') +
    facet_grid(. ~ clust_id) +
    scale_fill_manual(values = c(activated = 'firebrick',
                                 inhibited = 'steelblue'),
                      name = 'regulation') +
    scale_y_discrete(labels = globals$cohort_labs) +
    scale_x_continuous(labels = function(x) abs(x)) +
    globals$common_theme +
    theme(axis.title.y = element_blank()) +
    labs(title = 'Differentially regulated metabolic reactions',
         subtitle = paste('Cluster versus cohort mean, total reactions: n =',
                          bcg_meta$react_numbers$n_total[1]),
         x = '% of reactions')

# Bar plots with log2 fold-regulation estimates for selected subsystems --------

  insert_msg('Bar plots of fold-regulation estimates')

  for(i in names(bcg_metaplots$estimates)) {

    bcg_metaplots$bar_plots[[i]] <-
      list(x = blast(bcg_metaplots$estimates[[i]], subsystem),
           y = levels(bcg_metaplots$estimates[[i]]$subsystem)) %>%
      pmap(function(x, y) x %>%
             ggplot(aes(x = fold_reg,
                        y = reorder(react_id, fold_reg),
                        fill = regulation)) +
             geom_bar(stat = 'identity') +
             facet_grid(. ~ clust_id) +
             scale_fill_manual(values = c(activated = 'firebrick',
                                          inhibited = 'steelblue',
                                          ns = 'gray70'),
                               name = 'regulation\nvs hohort mean') +
             globals$common_theme +
             theme(axis.title.y = element_blank()) +
             labs(title = paste(y, globals$cohort_labs[i], sep = ', '),
                  x = expression('log'[2] * ' fold-regulation vs cohort mean')))

  }

# Bar plots for reactions regulated in both the TCGA and GSE99420 cohort-------

  insert_msg('Bar plots for the TCGA and GSE99420 cohorts')

  ## bar plot data, merging the steroid and sex hormone metabolism together

  bcg_metaplots$common_data <- bcg_metaplots$estimates[c("tcga", "gse99420")] %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = factor(cohort, c('tcga', 'gse99420'))) %>%
    blast(clust_id) %>%
    map2_dfr(bcg_metaplots$common_significant,
             ~filter(.x, react_id %in% .y)) %>%
    mutate(subsystem = ifelse(subsystem == 'Androgen and estrogen synthesis and metabolism',
                              'Steroid metabolism',
                              as.character(subsystem)),
           subsystem = factor(subsystem, bcg_metaplots$subsystems),
           subsystem = droplevels(subsystem),
           class = steroid_classifier(react_id))

  ## plots for particular subsystems

  bcg_metaplots$common_bar_plots <-
    list(x = blast(bcg_metaplots$common_data, subsystem),
         y = levels(bcg_metaplots$common_data$subsystem)) %>%
    pmap(function(x, y) x %>%
           ggplot(aes(x = fold_reg,
                      y = reorder(react_id, fold_reg),
                      fill = regulation)) +
           facet_grid(clust_id ~ cohort,
                      space = 'free_y',
                      scale = 'free_y',
                      labeller = labeller(.cols = globals$cohort_labs)) +
           geom_vline(xintercept = 0,
                      linetype = 'dashed') +
           geom_bar(stat = 'identity') +
           scale_fill_manual(values = c(activated = 'firebrick',
                                        inhibited = 'steelblue'),
                             name = 'Activity in cluster\nvs cohort mean') +
           scale_y_discrete(labels = metab_labeller) +
           globals$common_theme +
           theme(axis.title.y = element_blank()) +
           labs(title = y,
                x = expression('log'[2] * ' regulation vs cohort mean')))

  ## additional styling

  bcg_metaplots$common_bar_plots$`Oxidative phosphorylation` <-
    bcg_metaplots$common_bar_plots$`Oxidative phosphorylation` +
    theme(axis.title.y = element_text(size = 8, angle = 90)) +
    labs(y = 'mitochondrial complex')

  bcg_metaplots$common_bar_plots$`Fatty acid oxidation` <-
    bcg_metaplots$common_bar_plots$`Fatty acid oxidation` +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank())

# Rug plot with classification of steroid metabolism enzymes -------

  insert_msg('Rug plot for steroid metabolism enzymes')

  bcg_metaplots$rug_plot <- bcg_metaplots$common_data %>%
    filter(subsystem == 'Steroid metabolism') %>%
    ggplot(aes(x = 'reaction',
               y = reorder(react_id, fold_reg),
               fill = class)) +
    facet_grid(clust_id ~ .,
               space = 'free_y',
               scale = 'free_y') +
    geom_tile(color = 'black') +
    scale_fill_manual(values = c('androgens/testosterone' = 'plum4',
                                 'estrogens/progesteron' = 'aquamarine3',
                                 'steroids' = 'orangered3',
                                 'corticosteroids' = 'steelblue'),
                      name = 'Reaction\nclassification') +
    globals$common_theme +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank())

# END ------

  bcg_metaplots$subsystems <- NULL
  bcg_metaplots$estimates <- NULL
  bcg_metaplots$common_data <- NULL

  bcg_metaplots <- compact(bcg_metaplots)

  insert_tail()
