# Prediction of drug sensitivity based on cell line drug screening data obtained
# by the  CTRP and GDSC projects.
#
# Prediction of individual IC50 (GDSC) or AUC values for cancer samples with
# RIDGE regression
#
# The strategy:
#
# 1) selection of compounds and cell lines for the training data set.
# The cell lines are expected to be of epithelial origin; bone, nervous system,
# and blood cancers are excluded. In case of compounds shared by both the GDSC1
# and GDSC2 experiments, the GDSC2 data are preferred, as recommended by the
# project authors.
#
# 2) selection of explanatory factors: those are genes shared by the training
# data set and the bulk cancer data sets. Genes with 10% lowest median
# expression are removed. Next, numeric matrices
# with expression of the selected modeling variables are subjected to batch
# adjustment with ComBat.
#
# 3) Training of RIDGE models. Lambdas are selected with the minimal
# mean squared errors in cross-validation. Successful models are defined
# by cross-validated R^2 >= 0.13 and correlation of the predicted and observed
# outcome with r >= 0.5 (Cohen 1988 and Cohen 2013).
#
# 4) Subsequently, log IC50 [µM] and AUC are predicted for the cancer samples
# based on expression of the selected genes (ComBat-adjusted).


  insert_head()

# container -------

  drugs <- list()

# test data: gene expression ------

  insert_msg('Test data: gene expression')

  ## conversion to numeric matrices with genes in rows
  ## and observations in columns, data are already log2-transformed
  ## genes detected in all cohorts are used for modeling

  drugs$cancer_genes <- globals$cohort_expr %>%
    eval %>%
    map(~.x$annotation) %>%
    map(filter, !duplicated(gene_symbol)) %>%
    map(~.x$gene_symbol) %>%
    map(unname) %>%
    reduce(intersect)

  drugs$test_x <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    map(filter, tissue == 'tumor') %>%
    map(column_to_rownames, 'sample_id') %>%
    map(select, all_of(drugs$cancer_genes)) %>%
    map(t)

# Training data: reading the cached data -------

  insert_msg('Reading the cached training data')

  list(cache_path = c('./data/ctrp2.RData',
                      './data/gdsc.RData'),
       script_path = c('./import scripts/ctrp2.R',
                       './import scripts/gdsc.R'),
       message = c('Cached data for the CTRP2 experiment',
                   'Cached data for the GDSC experiments')) %>%
    pwalk(access_cache)

# Selection of the cell lines for training of the GLMNET models ------

  insert_msg('Cell line selection')

  drugs$sample_ids <- list(ctrp2 = ctrp2$cell_lexicon,
                           gdsc1 = gdsc$cell_lexicon$gdsc1,
                           gdsc2 = gdsc$cell_lexicon$gdsc2) %>%
    map(mutate, tissue = tolower(tissue)) %>%
    map(filter,
        !is.na(tissue),
        tissue %in% c('bowel',
                      'skin',
                      'bladder/urinary tract',
                      'lung',
                      'ovary/fallopian tube',
                      'breast',
                      'kidney',
                      'soft tissue',
                      'thyroid',
                      'prostate',
                      'esophagus/stomach',
                      'biliary tract',
                      'head and neck',
                      'uterus',
                      'liver',
                      'soft_tissue',
                      'aero_digestive_tract',
                      'urogenital_system',
                      'digestive_system')) %>%
    map(~.x$sample_id)

# Selection of drugs -------

  insert_msg('Selection of drugs')

  drugs$variables <- list(ctrp2 = ctrp2$drug_lexicon,
                          gdsc1 = gdsc$drug_lexicon$gdsc1,
                          gdsc2 = gdsc$drug_lexicon$gdsc2) %>%
    map(~.x$variable)

# Training data: drug sensitivity -------

  insert_msg('Training data: drug sensitivity')

  drugs$train_y <- list(ctrp2 = ctrp2$auc,
                        gdsc1 = gdsc$log_ic50$gdsc1,
                        gdsc2 = gdsc$log_ic50$gdsc2) %>%
    map(column_to_rownames, 'sample_id') %>%
    map(as.matrix)

  ## GDSC data sets: unit conversion:
  ## pM to avoid negative values upon log transformation

  drugs$train_y[c('gdsc1', 'gdsc2')] <-
    drugs$train_y[c('gdsc1', 'gdsc2')] %>%
    map(~log(exp(.x) * 1e6))

  ## samples and drugs of interest

  drugs$train_y <-
    list(x = drugs$train_y,
         y = drugs$sample_ids,
         z = drugs$variables) %>%
    pmap(function(x, y, z) x[y, z])

# Training data: gene expression --------

  insert_msg('Training data: gene expression')

  ## expression data are already log2 transformed.
  ## matrix format: genes in rows and samples in columns

  drugs$train_x <- list(ctrp2 = ctrp2$expression,
                        gdsc1 = gdsc$expression,
                        gdsc2 = gdsc$expression) %>%
    map(column_to_rownames, 'sample_id') %>%
    map(t)

# Selection of modeling genes and pre-processing -------

  insert_msg('Selection of explanatory factors')

  ## defining cutoffs for the median expression
  ## (genes with 10% lowest values are removed)

  drugs$data_objects <- drugs$train_x %>%
    map(multi_process,
        test = drugs$test_x,
        median_quantile = 0.1)

# cleaning the environment prior to modeling ------

  insert_msg('Cleaning the environment')

  drugs <- drugs[c("test_x", "train_x", "train_y", "data_objects")]

# Training and evaluation of the RIDGE models ---------

  insert_msg('Model selection, training and evaluation')

  for(i in names(drugs$train_y)) {

    plan('multisession')

    drugs$train_objects[[i]] <-
      train(x = drugs$data_objects[[i]],
            y = drugs$train_y[[i]],
            standardize = TRUE,
            alpha = 0,
            family = 'gaussian',
            type.measure = 'mse')

    plan('sequential')

  }

# Model evaluation stats -------

  drugs$stats <- drugs$train_objects %>%
    map(summary) %>%
    map(mutate,
        qc_passed = ifelse(rsq_oof >= 0.13 & pearson >= 0.5,
                           'yes', 'no'))

  ## plotting of the cross-validated R^2 and Pearson's correlation

  drugs$plots <- drugs$train_objects %>%
    map(plot) %>%
    map(~.x$rsq_spearman) %>%
    map(~.x +
          geom_vline(xintercept = 0.13,
                     linetype = 'dashed') +
          geom_hline(yintercept = 0.5,
                     linetype = 'dashed') +
          globals$common_theme)

  ## additional styling

  drugs$plots <-
    list(x = drugs$plots,
         y = paste('Model performance,',
                   c('CTRP2', 'GDSC1', 'GDSC2')),
         z = drugs$stats) %>%
    pmap(function(x, y, z) x  +
           labs(title = y,
                subtitle = paste0('total: n = ', nrow(z),
                                  ', QC passed: n = ', table(z$qc_passed)['yes'])))

# Compound lexicons ---------

  insert_msg('Compound lexicon for the QC-passing models')

  drugs$lexicons <- drugs$stats %>%
    map(filter, qc_passed == 'yes') %>%
    map(transmute,
        variable = response) %>%
    map2(.,
         list(ctrp2 = ctrp2$drug_lexicon,
              gdsc1 = gdsc$drug_lexicon$gdsc1,
              gdsc2 = gdsc$drug_lexicon$gdsc2),
         left_join, by = 'variable')

# Predictions -------

  insert_msg('Predictions of IC50')

  for(i in names(drugs$train_objects)) {

    drugs$predictions[[i]] <-
      predict(drugs$train_objects[[i]],
              newdata = drugs$data_objects[[i]],
              type = 'response')

  }

  ## bringing the GDSC predictions back to µM

  drugs$predictions[c("gdsc1", "gdsc2")] <-
    drugs$predictions[c("gdsc1", "gdsc2")] %>%
    map(map, ~exp(.x) * 1e-6) %>%
    map(map, log)

  ## as data frames

  drugs$predictions <- drugs$predictions %>%
    map(map, as.data.frame) %>%
    map(map, rownames_to_column, 'sample_id') %>%
    map(map, as_tibble)

# Caching -------

  insert_msg('Caching')

  rm(ctrp2, gdsc, i)

  drugs <-
    drugs[c("stats", "plots", "lexicons", "predictions")]

  save(drugs, file = './data/drugs.RData')

# END ------

  insert_tail()
