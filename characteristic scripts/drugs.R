# Comparison of drug sensitivity (log IC50 and AUC) between the hormonal
# clusters.
#
# Identification of differentially 'regulated' drugs follows a similar strategy
# to genes:
#
# 1) differences in drug sensitivity between the clusters are investigated first
# by one-way ANOVA with eta-square effect size statistic.
#
# 2) Difference in drug sensitivity as compared with the cohort mean are
# assessed by one-sample T test.
#
# 3) Differentially regulated drugs are defined by pFDR(ANOVA) < 0.05,
# eta-square >= 0.14, and pFDR(T-test) < 0.05.
#
# The GDSC1 and GDSC2 experiment-trained drug sensitivity estimates are
# analyzed jointly.

  insert_head()

# container -------

  bcg_drugs <- list()

# parallel backend --------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis globals -------

  insert_msg('Analysis globals')

  ## drugs: drugs with zero variance are excluded

  bcg_drugs$variables <- expl_drugs$top_variables

  ## drug sensitivity data for single cancer samples:
  ## AUC for the CTRP2 experiment-trained estimates
  ## log IC50 for the GDSC experiment-trained estimates

  bcg_drugs$data <-
    map2(drugs$predictions[names(bcg_drugs$variables)],
         bcg_drugs$variables,
         function(data, variable) data %>%
           map(select, sample_id, all_of(variable))) %>%
    map(function(data) map2(bcg_globals$assignment[names(data)],
                            data,
                            inner_join, by = 'sample_id')) %>%
    map(map, column_to_rownames, 'sample_id')

# Descriptive stats -------

  insert_msg('Descriptive stats')

  bcg_drugs$stats <- bcg_drugs$data %>%
    map(map,
        fast_num_stats,
        split_factor = 'clust_id')

# one-way ANOVA -------

  insert_msg('One-way ANOVA')

  for(i in names(bcg_drugs$data)) {

    bcg_drugs$anova[[i]] <- bcg_drugs$data[[i]] %>%
      future_map(test_anova,
                 split_fct = 'clust_id',
                 variables = bcg_drugs$variables[[i]],
                 adj_method = 'BH',
                 .options = furrr_options(seed = TRUE)) %>%
      map(~.x$anova) %>%
      map(mutate,
          variable = response,
          drug_name = exchange(variable,
                               drugs$lexicons[[i]],
                               value = 'drug_name'))

    ## identification of the significant differences between
    ## the clusters

    bcg_drugs$anova_significant[[i]] <- bcg_drugs$anova[[i]] %>%
      map(filter,
          p_adjusted < 0.05,
          effect_size >= 0.14) %>%
      map(~.x$variable)

  }

# Differences between the cluster means and the grand averages -------

  insert_msg('Differences between the clusters and the cohort means')

  for(i in names(bcg_drugs$data)) {

    bcg_drugs$dev_test[[i]] <- bcg_drugs$data[[i]] %>%
      map(avg_deviation,
          split_fct = 'clust_id',
          variables = bcg_drugs$variables[[i]])

    ## formatting the results and regulation signs

    bcg_drugs$dev_test[[i]] <-
      map2(bcg_drugs$dev_test[[i]],
           bcg_drugs$anova_significant[[i]],
           ~mutate(.x,
                   anova_significant = ifelse(variable %in% .y,
                                              'yes', 'no'))) %>%
      map(mutate,
          regulation = ifelse(p_adjusted >= 0.05 | anova_significant == 'no',
                              'ns',
                              ifelse(deviation_center > 0, 'resistant',
                                     ifelse(deviation_center < 0,
                                            'sensitive', 'ns'))),
          regulation = factor(regulation,
                              c('resistant', 'sensitive', 'ns')),
          drug_name = exchange(variable,
                               drugs$lexicons[[i]],
                               value = 'drug_name'))


  }

# Significant differences in drug sensitivity --------

  insert_msg('Significant effects')

  for(i in names(bcg_drugs$dev_test)) {

    ## in single cohorts

    bcg_drugs$significant[[i]] <- bcg_drugs$dev_test[[i]] %>%
      map(filter, regulation %in% c('resistant', 'sensitive')) %>%
      map(blast, clust_id) %>%
      transpose %>%
      map(map, blast, regulation) %>%
      map(transpose) %>%
      map(map, map, ~.x$variable)

    ## common effects: shared by the TCGA and GSE999420 cohorts

    bcg_drugs$common_significant[[i]] <- bcg_drugs$significant[[i]] %>%
      map(map, ~.x[c('tcga', 'gse99420')]) %>%
      map(map, reduce, intersect)

  }

# Numbers of differentially regulated drugs -------

  insert_msg('Numbers of differentially regulated drugs')

  bcg_drugs$n_numbers$totals <- bcg_drugs$variables %>%
    map_dbl(length) %>%
    compress(names_to = 'experiment',
             values_to = 'n_total')

  bcg_drugs$n_numbers$drugs <- bcg_drugs$significant %>%
    map(map, map, map_dbl, length) %>%
    map(map, map, compress,
        names_to = 'cohort',
        values_to = 'n')%>%
    map(map, compress, names_to = 'status') %>%
    map(compress, names_to = 'clust_id') %>%
    compress(names_to = 'experiment')

  bcg_drugs$n_numbers <-
    left_join(bcg_drugs$n_numbers$drugs,
              bcg_drugs$n_numbers$totals,
              by = 'experiment') %>%
    mutate(experiment = factor(experiment, names(bcg_drugs$data)),
           percent = n/n_total * 100)

# END ----

  rm(i)

  bcg_drugs <- bcg_drugs[c("stats", "anova", "anova_significant",
                           "dev_test", "significant", "common_significant",
                           "n_numbers")]

  plan('sequential')

  insert_tail()
