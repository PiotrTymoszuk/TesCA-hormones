# Modeling of progression-free survival with RIDGE Cox regression in the
# TCGA cohort.
# Candidate explanatory factors are: age, detailed histology, marker status,
# and, in one of the models, hormonal cluster assignment.
# pT stage was not included in the explanatory variable set because of
# insufficient numbers of progressions in late stage cancers,
# which leads to unexpected HR values
# (stages II and III as favorable explanatory factors).
# Histology: yolk sac tumors are collapsed with embryonic carcinomas.

  insert_head()

# container ------

  bcg_glmcox <- list()

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis data and globals -------

  insert_msg('Analysis data and globals')

  ## analysis data, collapsing stages, histologies and marker status
  ## setting the intercept for the hormonal clusters to cluster #2,
  ## which is expected to have favorable prognosis. Clusters #1 and #3
  ## are collapsed

  bcg_glmcox$clinical_variables <- c("age",
                                     #"histology",
                                     "histology_icd",
                                     #"pt_stage",
                                     "marker_status")

  bcg_glmcox$data <-
    inner_join(bcg_globals$assignment$tcga,
               tcga$clinic[, c("sample_id",
                               "pfs_days",
                               "progression",
                               bcg_glmcox$clinical_variables)],
               by = 'sample_id') %>%
    filter(complete.cases(.)) %>%
    mutate(histology_icd = fct_collapse(histology_icd,
                                        TT = c('benign teratoma',
                                               'malignant teratoma',
                                               'teratocarcinoma'),
                                        `EMBCA/TYST` = c('yolk sac cancer',
                                                       'embryonal carcinoma'),
                                        MGCT = c('germinal mixed histology'),
                                        SEM = c('seminoma')),
           histology_icd = fct_relevel(histology_icd, 'SEM'),
           marker_status = fct_collapse(marker_status,
                                        `S2+` = c('S2', 'S3')),
           clust_id = fct_collapse(clust_id,
                                   `#1/#3` = c('#1', '#3')),
           clust_id = fct_relevel(clust_id, '#2'))

  ## Y objects: survival objects
  ##
  ## X matrices: matrices with explanatory factors in three flavors:
  ## for the full model with the cluster assignment and clinical predictors,
  ## for the model with clinical explanatory variables only,
  ## for the model with clustter assignment as the sole explanatory factor

  bcg_glmcox$y <-
    Surv(bcg_glmcox$data$pfs_days, bcg_glmcox$data$progression)

  bcg_glmcox$x$full <- bcg_glmcox$data %>%
    select(-pfs_days, -progression) %>%
    column_to_rownames('sample_id')

  bcg_glmcox$x$clinical <- bcg_glmcox$x$full %>%
    select(-clust_id)

  bcg_glmcox$x$clusters <- bcg_glmcox$x$full %>%
    select(clust_id)

  bcg_glmcox$x <- bcg_glmcox$x[c('clusters',
                                 'clinical',
                                 'full')] %>%
    map(~model.matrix(~ .,  data = .x))

  ## CV folds

  bcg_glmcox$n_rep <- 1:200 %>%
    set_names(paste0('rep_', 1:200))

  set.seed(12345)

  bcg_glmcox$fold_id <- bcg_glmcox$n_rep %>%
    map(function(x) createFolds(bcg_glmcox$data$progression,
                                k = 5,
                                list = FALSE,
                                returnTrain = TRUE))

  ## model labels and colors

  bcg_glmcox$model_labs <- c(clusters = 'Cluster-only PFS model',
                             clinical = 'Clinical PFS model',
                             full = 'Cluster/clinical PFS model')

  bcg_glmcox$model_colors <- c(clusters = 'darkorange2',
                               clinical = 'darkseagreen3',
                               full = 'steelblue4')

  bcg_glmcox$model_legends <- c(clusters = 'clusters',
                                clinical = 'clinical factors',
                                full = 'clusters\n+clinical factors')

# Numbers of patients and progression cases -------

  insert_msg('Numbers of patients and progression cases')

  ## and ready-to-use plot captions

  bcg_glmcox$n_caption <-
    paste0('total: n = ', nrow(bcg_glmcox$data),
           ', events: n = ', sum(bcg_glmcox$data$progression))

# Tuning of the GLMNET cox model ------

  insert_msg('Model tuning')

  for(i in names(bcg_glmcox$x)) {

    bcg_glmcox$tune_models[[i]] <- bcg_glmcox$fold_id %>%
      future_map(~cv.glmnet(x = bcg_glmcox$x[[i]],
                            y = bcg_glmcox$y,
                            type.measure = 'default',
                            family = 'cox',
                            foldid = .x,
                            alpha = 0),
                 .options = furrr_options(seed = TRUE))

  }

  ## lambda corresponding to the lowest deviance

  bcg_glmcox$tune_stats <- bcg_glmcox$tune_models %>%
    map(map, ~.x[c('lambda', 'cvm', 'nzero')]) %>%
    map(map, as_tibble) %>%
    map(map, filter, cvm == min(cvm)) %>%
    map(map_dfr, ~.x[1, ])

  ## the best tune

  bcg_glmcox$best_tune <- bcg_glmcox$tune_stats %>%
    map(filter, cvm == min(cvm))

  bcg_glmcox$best_tune <- bcg_glmcox$best_tune %>%
    map(~.x[1, ])

  bcg_glmcox$tune_stats <- bcg_glmcox$tune_stats %>%
    map(mutate,
        best = ifelse(cvm == min(cvm),
                      'yes', 'no'))

  ## tuning process plot

  bcg_glmcox$tune_plots <-
    list(x = bcg_glmcox$tune_stats,
         y = paste(bcg_glmcox$model_labs,
                   globals$cohort_labs["tcga"],
                   sep = ', '),
         z = bcg_glmcox$best_tune) %>%
    pmap(function(x, y, z) x %>%
           ggplot(aes(x = lambda,
                      y = cvm,
                      fill = best)) +
           geom_point(shape = 21,
                      size = 2) +
           scale_fill_manual(values = c('no' = 'steelblue',
                                        'yes' = 'orangered2')) +
           globals$common_theme +
           labs(title = y,
                subtitle = paste0('minimal error: dev = ', signif(z$cvm, 3),
                                  ', \u03BB = ', signif(z$lambda, 3)),
                x = expression(lambda),
                y = 'deviance, 200-repeats 10-fold CV'))

# Training a GLMNET model for the optimal lambda -------

  insert_msg('Training the GMLNET models')

  bcg_glmcox$glmnet_models <-
    list(x = bcg_glmcox$x,
         lambda = map(bcg_glmcox$best_tune, ~.x$lambda)) %>%
    pmap(glmnet,
         y = bcg_glmcox$y,
         alpha = 0,
         family = 'cox')

# Linear predictor scores and an univariable Cox models for evaluation -------

  insert_msg('Linear predictor scores, univariable cox models')

  bcg_glmcox$lp_scores <-
    list(object = bcg_glmcox$glmnet_models,
         newx = bcg_glmcox$x) %>%
    pmap(predict) %>%
    map(set_colnames, 'lp_score') %>%
    map(as.data.frame) %>%
    map(rownames_to_column, 'sample_id') %>%
    map(left_join,
        bcg_glmcox$data[c('sample_id', 'pfs_days', 'progression')],
        by = 'sample_id') %>%
    map(as_tibble)

  bcg_glmcox$cox_models <- bcg_glmcox$lp_scores %>%
    map(~call2('coxph',
               formula = Surv(pfs_days, progression) ~ lp_score,
               data = .x,
               x = TRUE,
               y = TRUE)) %>%
    map(eval) %>%
    map2(bcg_glmcox$lp_scores,
         as_coxex)

# Model assumptions and fit stats --------

  insert_msg('Model assumptions and fit stats')

  ## the proportional hazard assumption is met

  bcg_glmcox$assumptions <- bcg_glmcox$cox_models %>%
    map(summary, 'assumptions')

  bcg_glmcox$stats <- bcg_glmcox$cox_model %>%
    map(summary, 'fit') %>%
    compress(names_to = 'model_type') %>%
    relocate(model_type) %>%
    mutate(model_type = factor(model_type, names(bcg_glmcox$cox_models)))

# Bubble plot of the numeric stats of model performance -------

  insert_msg('Bubble plots of model stats')

  bcg_glmcox$stat_plot <- bcg_glmcox$stats %>%
    ggplot(aes(x = c_index,
               y = 1 - ibs_model,
               size = raw_rsq,
               fill = model_type,
               color = model_type)) +
    geom_vline(xintercept = 0.5,
               linetype = 'dashed') +
    geom_hline(yintercept = 0.75,
               linetype = 'dashed') +
    geom_point(shape = 21,
               color = 'black') +
    geom_errorbarh(aes(xmin  = lower_ci,
                       xmax = upper_ci),
                   height = 0,
                   linewidth = 0.75,
                   show.legend = FALSE) +
    scale_fill_manual(labels = bcg_glmcox$model_legends,
                      values = bcg_glmcox$model_colors,
                      name = '') +
    scale_color_manual(labels = c(clusters = 'clusters',
                                  clinical = 'clinical factors',
                                  full = 'clinical factors + clusters'),
                       values = bcg_glmcox$model_colors,
                       name = '') +
    scale_size_area(max_size = 4,
                    limits = c(0, 0.2),
                    name = expression('R'^2)) +
    globals$common_theme +
    labs(title = 'PFS model performance',
         subtitle = bcg_glmcox$n_caption,
         x = 'C-index, 95% CI',
         y = '1 - IBS')

# Brier scores for the unique time points ------

  insert_msg('Brier scores for the unique time points')

  ## Brier score frame

  bcg_glmcox$brier_scores <- bcg_glmcox$cox_models %>%
    map(surv_brier)

  bcg_glmcox$brier_scores$clusters <- bcg_glmcox$brier_scores$clusters %>%
    select(time, reference, training) %>%
    set_names(c('time', 'reference', 'clusters'))

  bcg_glmcox$brier_scores[c("clinical", "full")] <-
    bcg_glmcox$brier_scores[c("clinical", "full")] %>%
    map(select, time, training) %>%
    map2(., names(.),
         ~set_names(.x, c('time', .y)))

  bcg_glmcox$brier_scores <- bcg_glmcox$brier_scores %>%
    reduce(left_join, by = 'time')

  ## Brier score trajectories

  bcg_glmcox$brier_plot <- bcg_glmcox$brier_scores %>%
    pivot_longer(cols = c(reference,
                          clusters,
                          clinical,
                          full),
                 names_to = 'model_type',
                 values_to = 'brier_score') %>%
    mutate(model_type =  factor(model_type, names(bcg_glmcox$cox_models))) %>%
    ggplot(aes(x = time,
               y = brier_score,
               color = model_type)) +
    geom_line() +
    scale_color_manual(labels = c(bcg_glmcox$model_legends,
                                  reference = 'NULL model'),
                       values = c(bcg_glmcox$model_colors,
                                  reference = 'gray70'),
                       name = 'Explanatory factors') +
    globals$common_theme +
    labs(title = 'RIDGE Cox model calibration',
         subtitle = bcg_glmcox$n_caption,
         x = 'progression-free survival, days',
         y = 'Brier score')

# Survival in tertiles of linear predictor scores -------

  insert_msg('Survival in tertiles of linear predictor scores')

  ## done for the clinical and full models

  bcg_glmcox$calibrator_objects <-
    bcg_glmcox$cox_models[c("clinical",
                            "full")] %>%
    future_map(calibrate.coxex,
               n = 3,
               labels = paste0('Q', 1:3),
               .options = furrr_options(seed = TRUE))

  ## Kaplan-Meier plots

  bcg_glmcox$km_plots <- bcg_glmcox$calibrator_objects %>%
    map(plot,
        cust_theme = globals$common_theme,
        palette = c('firebrick1',
                    'firebrick4',
                    'black'))

  bcg_glmcox$km_plots <-
    list(x = bcg_glmcox$km_plots,
         y = paste(bcg_glmcox$model_labs[c("clinical",
                                           "full")],
                   globals$cohort_labs["tcga"],
                   sep = ', ')) %>%
    pmap(function(x, y) x +
           labs(title = y,
                subtitle = paste(bcg_glmcox$n_caption,
                                 x$labels$subtitle,
                                 sep = '\n'),
                x = 'progression-free survival, days') +
           theme(plot.tag = element_blank()))

# Model coefficient estimates -------

  insert_msg('Model coefficient estimates')

  bcg_glmcox$coefs <- bcg_glmcox$glmnet_models %>%
    map(coef) %>%
    map(as.matrix) %>%
    map(as.data.frame) %>%
    map(set_names, 'beta') %>%
    map(rownames_to_column, 'parameter') %>%
    map(filter, beta != 0) %>%
    map(as_tibble)

  ## extraction of the variables and levels,
  ## ready-to-use plot axis labels

  bcg_glmcox$coefs <- bcg_glmcox$coefs %>%
    map(mutate,
        hr = exp(beta),
        variable = stri_extract(parameter,
                                regex = paste(sort(c(bcg_glmcox$clinical_variables,
                                                     'clust_id'),
                                                   decreasing = TRUE),
                                              collapse = '|')),
        level = stri_replace(parameter,
                             regex = paste(sort(c(bcg_glmcox$clinical_variables,
                                                  'clust_id'),
                                                decreasing = TRUE),
                                           collapse = '|'),
                             replacement = ''),
        axis_lab = car::recode(variable,
                               "'clust_id' = 'cluster';
                               'histology_icd' = '';
                               'pt_stage' = 'stage';
                               'age' = 'age, years';
                               'marker_status' = 'markers'"),
        axis_lab = ifelse(axis_lab == '',
                          level,
                          ifelse(level == '',
                                 axis_lab,
                                 paste(axis_lab, level, sep = ': '))))

# Bar plots of the HRs ------

  insert_msg('Bar plots of HRs')

  bcg_glmcox$coef_plots <-
    list(x = bcg_glmcox$coefs,
         y = paste(bcg_glmcox$model_labs,
                   globals$cohort_labs["tcga"],
                   sep = ', ')) %>%
    pmap(function(x, y) x %>%
           ggplot(aes(x = beta,
                      y = reorder(axis_lab, beta),
                      fill = factor(sign(beta)),
                      color = factor(sign(beta)))) +
           geom_vline(xintercept = 0,
                      linetype = 'dashed') +
           geom_bar(stat = 'identity',
                    color = 'black') +
           geom_text(aes(label = signif(hr, 3),
                         #x = beta * 1.05,
                         hjust = ifelse(beta < 0, 1.3, -1.3)),
                     size = 2.5,
                     show.legend = FALSE) +
           scale_fill_manual(values = c('-1' = 'steelblue',
                                        '0' = 'gray60',
                                        '1' = 'firebrick'),
                             labels = c('-1' = 'favorable',
                                        '0' = 'ns',
                                        '1' = 'unfavorable'),
                             name = 'Association\nwith risk') +
           scale_color_manual(values = c('-1' = 'steelblue',
                                         '0' = 'gray60',
                                         '1' = 'firebrick'),
                              labels = c('-1' = 'favorable',
                                         '0' = 'ns',
                                         '1' = 'unfavorable'),
                              name = 'Association\nwith progression risk') +
           globals$common_theme +
           theme(axis.title.y = element_blank()) +
           labs(title = y,
                x = expression('log HR'[RIDGE])))

# END -----

  bcg_glmcox$data <- NULL
  bcg_glmcox$y <- NULL
  bcg_glmcox$x <- NULL
  bcg_glmcox$fold_id <- NULL

  bcg_glmcox <- compact(bcg_glmcox)

  rm(i)

  plan('sequential')

  insert_tail()
