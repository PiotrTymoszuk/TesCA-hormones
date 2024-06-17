# Plots of the top enriched GO terms in the hormonal clusters

  insert_head()

# container --------

  bcg_goplots <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  ## top enriched GO terms: enrichment OR and pFDR

  bcg_goplots$top_gos <-
    map2(bcg_go$test,
         bcg_go$significant,
         function(x, y) map2(x, y,
                             ~filter(.x, term %in% .y))) %>%
    map(map,
        select,
        go_id, term,
        or, p_adjusted) %>%
    map(map,
        top_n,
        n = 20,
        or)

  ## common range of p values and OR

  bcg_goplots$p_val_range <- bcg_goplots$top_gos %>%
    map(map, ~.x$p_adjusted) %>%
    unlist %>%
    range

  bcg_goplots$or_range <- bcg_goplots$top_gos %>%
    map(map, ~.x$or) %>%
    unlist %>%
    range

# numbers of significantly enriched GO terms -------

  insert_msg('Numbers of significantly enricged GO terms')

  ## ready-to-use plot caption/labeller with the total numbers of GOs

  bcg_goplots$n_caption <- bcg_go$go_numbers %>%
    filter(cohort != 'common',
           !duplicated(cohort))

  bcg_goplots$n_caption <-
    map2_chr(bcg_goplots$n_caption$cohort,
             bcg_goplots$n_caption$n_total,
             ~paste(globals$cohort_labs[.x], .y, sep = '\ntotal GO: n = ')) %>%
    set_names(bcg_goplots$n_caption$cohort)

  ## bar plots

  bcg_goplots$n_numbers <- bcg_go$go_numbers %>%
    filter(cohort != 'common') %>%
    mutate(percent = n/n_total * 100,
           cohort = factor(cohort, globals$analysis_cohorts)) %>%
    ggplot(aes(x = percent,
               y = factor(clust_id, rev(levels(clust_id))),
               fill = clust_id)) +
    geom_bar(stat = 'identity',
             color = 'black') +
    facet_grid(. ~ cohort,
               labeller = as_labeller(bcg_goplots$n_caption)) +
    scale_fill_manual(values = globals$cluster_colors) +
    guides(fill = 'none') +
    globals$common_theme +
    labs(title = 'Numbers of enriched GO terms',
         x = '% of all terms',
         y = 'cluster')

# Top enriched GOs in the hormonal clusters --------

  insert_msg('Top enriched GO terms in the clusters')

  for(i in names(bcg_goplots$top_gos)) {

    bcg_goplots$top_plots[[i]] <-
      list(x = bcg_goplots$top_gos[[i]],
           y = paste('Top GO terms',
                     i,
                     globals$cohort_labs[names(bcg_goplots$top_gos[[i]])],
                     sep = ', ')) %>%
      pmap(function(x, y) x %>%
             ggplot(aes(x = or,
                        y = reorder(term, or),
                        size = -log10(p_adjusted),
                        fill = -log10(p_adjusted))) +
             geom_point(shape = 21) +
             globals$common_theme +
             scale_size_area(max_size = 6,
                             limits = -log10(rev(bcg_goplots$p_val_range)),
                             name = expression('-log'[10] * 'p FDR')) +
             scale_fill_gradient(low = 'white',
                                 high = 'firebrick',
                                 limits = -log10(rev(bcg_goplots$p_val_range)),
                                 name = expression('-log'[10] * 'p FDR')) +
             scale_x_continuous(limits = bcg_goplots$or_range) +
             guides(size = 'legend',
                    fill = 'legend') +
             theme(axis.title.y = element_blank()) +
             labs(title = y,
                  x = 'OR, enrichment over genome'))

  }

# END -----

  bcg_goplots$top_gos <- NULL
  bcg_goplots <- compact(bcg_goplots)

  rm(i)

  insert_tail()
