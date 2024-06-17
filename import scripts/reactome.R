# ssGSEA scores of the Reactome pathway gene signatures

  insert_head()

# container --------

  reactome <- list()

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# signature database --------

  insert_msg('Signature database')

  db <- load_dbsig('./data/signatures/msigdb.v7.5.1.symbols.gmt')

  reactome$lexicon <- db %>%
    filter(stri_detect(sign_name, regex = '^REACTOME'))

# input expression data --------

  insert_msg('Input expression data')

  exp_data <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    map(filter, tissue == 'tumor') %>%
    map(select, -tissue) %>%
    map(column_to_rownames, 'sample_id')

# ssGSEA scores ---------

  insert_msg('ssGSEA scores')

  reactome$scores <- exp_data %>%
    future_map(calculate.dbsig,
               x = reactome$lexicon,
               .options = furrr_options(seed = TRUE))

  ## appending the score frames with sample IDs

  reactome$scores <-
    map2(reactome$scores,
         exp_data,
         ~mutate(.x, sample_id = rownames(.y))) %>%
    map(relocate, sample_id)

  ## updating the signature lexicon

  reactome$lexicon <- reactome$lexicon %>%
    mutate(variable = sign_name,
           label = stri_replace(sign_name,
                                regex = '^REACTOME_',
                                replacement = ''),
           label = stri_replace_all(label,
                                    fixed = '_',
                                    replacement = ' ')) %>%
    filter(variable %in% reduce(map(reactome$scores, names), intersect))

# Caching --------

  insert_msg('Caching')

  save(reactome, file = './data/reactome.RData')

# END -------

  plan('sequential')

  rm(exp_data)

  insert_tail()


