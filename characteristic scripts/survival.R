# Relapse- and progression-free survival in the hormonal clusters.
# The analysis is done for the TCGA cohort.
# Survival is compared between the clusters by Pato-Peto tests: global and
# for the cluster pairs

  insert_head()

# container ------

  bcg_surv <- list()

# analysis data -------

  insert_msg('Analysis data')

  bcg_surv$data <- tcga$clinic[c("sample_id",
                                 "rfs_days", "relapse",
                                 "pfs_days", "progression")]

  bcg_surv$data <-
    inner_join(bcg_surv$data,
               bcg_globals$assignment$tcga,
               by = 'sample_id') %>%
    filter(complete.cases(.))

# Numbers of cases and events --------

  insert_msg('Numbers of cases and events')

  ## global numbers of patients, relapses and progression cases
  ## ready-to-use subtitles for Kaplan-Meier plots

  for(i in list(c('rfs', 'relapse'),
                c('pfs', 'progression'))) {

    bcg_surv$n_numbers[[i[[1]]]] <-
      tibble(n_total = nrow(bcg_surv$data),
             n_events = sum(bcg_surv$data[[i[[2]]]]))

  }

  bcg_surv$n_captions <- bcg_surv$n_numbers %>%
    map(~paste0('total: n = ', .x[[1]], ', events: n = ', .x[[2]]))

  ## numbers of patients and events in the hormonal clusters
  ## ready-to-use legend labels for KM plots

  for(i in list(c('rfs', 'relapse'),
                c('pfs', 'progression'))) {

    bcg_surv$n_cluster[[i[[1]]]] <- bcg_surv$data %>%
      blast(clust_id) %>%
      map(~tibble(n_total = nrow(.x),
                  n_events = sum(.x[[i[[2]]]]))) %>%
      compress(names_to = 'clust_id') %>%
      mutate(legend_lab = paste0(clust_id, '\ntotal: n = ',
                                 n_total, '\nevents: n = ',
                                 n_events))

  }

  bcg_surv$legend_labs <- bcg_surv$n_cluster %>%
    map(~set_names(.x$legend_lab, .x$clust_id))

# surv-fit objects: globals and for pairwise comparisions -------

  insert_msg('Surv-fit objects')

  bcg_surv$survfit_objects <-
    list(rfs = Surv(rfs_days, relapse) ~ clust_id,
         pfs = Surv(pfs_days, progression) ~ clust_id) %>%
    map(make_survfit,
        data = bcg_surv$data)

# Median survival times in the clusters ---------

  insert_msg('Median survival times')

  bcg_surv$stats <- bcg_surv$survfit_objects %>%
    map(~.x$global) %>%
    map(surv_median) %>%
    map(mutate,
        clust_id = stri_replace(strata,
                                regex = '^clust_id=',
                                replacement = '')) %>%
    map(as_tibble)

# Tests -------

  insert_msg('Tests')

  ## adjustment of the post-hoc tests with the FDR method

  bcg_surv$test <- bcg_surv$survfit_objects %>%
    map(surv_pvalue, method = 'S1') %>%
    map(compress, names_to = 'comparison') %>%
    map(mutate,
        p_adjusted = NA,
        significance = NA) %>%
    map(as_tibble)

  for(i in names(bcg_surv$test)) {

    bcg_surv$test[[i]][-1, ] <- bcg_surv$test[[i]][-1, ] %>%
      re_adjust('pval', method = 'BH')

    ## the test for the global effects remain unchanged,
    ## we just need a nicely formatted significance label

    bcg_surv$test[[i]][1, ] <- bcg_surv$test[[i]][1, ] %>%
      re_adjust('pval', method = 'none')

  }

  bcg_surv$test <- bcg_surv$test %>%
    map(mutate,
        cluster1 = stri_split_fixed(comparison,
                                    pattern = '|',
                                    simplify = TRUE)[, 1],
        cluster2= stri_split_fixed(comparison,
                                    pattern = '|',
                                    simplify = TRUE)[, 2],
        cluster1 = factor(cluster1, levels(bcg_surv$data$clust_id)),
        cluster2 = factor(cluster2, levels(bcg_surv$data$clust_id)))

# Kaplan-Meier plots ------

  insert_msg('KM plots')

  bcg_surv$plots <-
    list(fit = map(bcg_surv$survfit_objects, ~.x$global),
         title = paste(globals$cohort_labs["tcga"],
                       c('relapse-free survival',
                         'progression-free survival')),
         xlab = c('relapse-free survival, days',
                  'progression-free survival, days'),
         legend.labs = map(bcg_surv$legend_labs, unname),
         pval = bcg_surv$test %>%
           map(filter,
               comparison == 'global') %>%
           map(~.x$significance)) %>%
    pmap(ggsurvplot,
         palette = unname(globals$cluster_colors),
         pval.size = 2.75,
         legend.title = 'Cluster') %>%
    map2(., bcg_surv$n_captions,
         ~.x$plot +
           labs(subtitle = .y) +
           globals$common_theme)

# Heat maps for pairwise comparisons between the clusters -------

  insert_msg('Heat maps for pairwise comparisons between the clusters')

  bcg_surv$post_hoc_plots <- bcg_surv$test %>%
    map(filter, comparison != 'global') %>%
    map(mutate,
        significant = ifelse(p_adjusted < 0.05, 'significant', 'ns'),
        significant = factor(significant, c('significant', 'ns'))) %>%
    map2(.,
         paste(globals$cohort_labs["tcga"],
               c('relapse-free survival',
                 'progression-free survival')),
         ~ggplot(.x,
                aes(x = cluster1,
                    y = cluster2,
                    fill = significant)) +
          geom_tile(color = 'black') +
          geom_text(aes(label = significance),
                    size = 2.75) +
          scale_fill_manual(values = c(significant = 'indianred',
                                       ns = 'gray80')) +
          globals$common_theme +
           labs(title = .y,
                subtitle = 'Peto-Peto post-hoc tests',
                x = 'cluster',
                y = 'cluster'))

# Result summary tables -------

  insert_msg('Result summary tables')

  bcg_surv$result_tbl <- bcg_surv$n_cluster %>%
    map(select, clust_id, n_total, n_events) %>%
    map2(.,
         map(bcg_surv$stats,
             select, median, lower, upper),
         cbind) %>%
    map2(.,
         map(bcg_surv$test,
             filter, comparison == 'global'),
         ~mutate(.x,
                 significance = .y$significance))

  bcg_surv$result_tbl <- bcg_surv$result_tbl %>%
    compress(names_to = 'response') %>%
    mutate(response = car::recode(response,
                                  "'rfs' = 'relapse-free survival';
                                  'pfs' = 'progression-free survival'"),
           median = paste0(signif(median, 2), ' [95% CI: ',
                           signif(lower, 2), ' to ',
                           signif(upper, 2), ']')) %>%
    select(response, clust_id,
           n_total, n_events,
           median, significance) %>%
    set_names(c('Survival response',
                'Hormonal cluster',
                'Total patients, N',
                'Events, N',
                'Median survival, days',
                'Significance'))

# END -------

  bcg_surv$data <- NULL
  bcg_surv$survfit_objects <- NULL

  bcg_surv <- compact(bcg_surv)

  rm(i)

  insert_tail()
