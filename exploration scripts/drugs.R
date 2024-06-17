# Exploratory analysis for drug sensitivity estimates, i.e. AUC for the CTRP2
# and log IC50 for the GDSC studies.

  insert_head()

# container -------

  expl_drugs <- list()

# analysis globals --------

  insert_msg('Analysis globals')

  expl_drugs$data <- drugs$predictions %>%
    map(map, column_to_rownames, 'sample_id') %>%
    map2(., map(drugs$lexicons, ~.x$variable),
         function(data, variable) data %>%
           map(select, all_of(variable)))

  ## log transformation of the CTRP2 data

  #expl_drugs$data$ctrp2 <- expl_drugs$data$ctrp2 %>%
   # map(map_dfc, log)

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# Distribution stats -------

  insert_msg('Distribution stats')

  ## mean, SD, median, interquartile ranges, ranges

  expl_drugs$stats <- expl_drugs$data %>%
    map(map,
        ~tibble(variable = names(.x),
                mean = colMeans(.x),
                sd = colSDs(.x),
                min = colMins(.x),
                max = colMax(.x),
                median = colMedians(.x),
                q25 = colQuantiles(.x, 0.25)[, 1],
                q75 = colQuantiles(.x, 0.75)[, 1])) %>%
    map(map,
        mutate,
        diff_mean_median = mean - median)

# Looking for non-variant drugs -------

  insert_msg('Looking for non-variant drugs')

  expl_drugs$nzv <- expl_drugs$data %>%
    map(map, distr_stats)

  for(i in names(expl_drugs$nzv)) {

    expl_drugs$top_variables[[i]] <- expl_drugs$nzv[[i]] %>%
      map(filter, !nzv) %>%
      map(~.x$variable) %>%
      reduce(intersect)

  }

# Normality of distribution -------

  insert_msg('Normality of distribution')

  ## tested for variant drugs only

  for(i in names(expl_drugs$data)) {

    expl_drugs$normality[[i]] <- expl_drugs$data[[i]] %>%
      future_map(~explore(.x,
                          variables = expl_drugs$top_variables[[i]],
                          what = 'normality',
                          pub_styled = FALSE),
                 .options = furrr_options(seed = TRUE))

  }

  ## counting normality faults: W < 0.95
  ##
  ## the GDSC data (log IC50) are more or less nicely normally distributed
  ## for up to 52 drugs in the CTRP2 experiment, lacking normality may be
  ## a problem. This can not be mended by log or sqrt transformation

  expl_drugs$normality <- expl_drugs$normality %>%
    map(map,
        select,
        variable, stat, p_value)

  expl_drugs$non_normal <- expl_drugs$normality %>%
    map(map, filter, stat <= 0.95) %>%
    map(map, ~.x$variable)

# END ------

  rm(i)

  expl_drugs$data <- NULL
  expl_drugs <- compact(expl_drugs)

  plan('sequential')

  insert_tail()

