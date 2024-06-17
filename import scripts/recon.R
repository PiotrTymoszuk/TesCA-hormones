# ssGSEA scores of the Recon metabolic subsystem gene signatures

  insert_head()

# container -----

  recon <- list()

# parallel backend --------

  insert_msg('Parallel backend')

  plan('multisession')

# gene signature list ---------

  insert_msg('Gene signature list')

  subs <- extract_subsystems(Recon2D)
  genes <- extract_genes(Recon2D)

  recon$lexicon <-
    inner_join(subs,
               genes[c('react_id', 'entrez_id')],
               by = 'react_id') %>%
    blast(subsystem) %>%
    map(~.x$entrez_id) %>%
    map(reduce, union) %>%
    compress(names_to = 'label',
             values_to = 'entrez_id') %>%
    mutate(variable = make.names(label),
           gene_symbol = map(entrez_id,
                             ~mapIds(org.Hs.eg.db,
                                     keys = .x,
                                     keytype = 'ENTREZID',
                                     column = 'SYMBOL')),
           gene_symbol = unname(gene_symbol))

# input expression data --------

  insert_msg('Input expression data')

  exp_data <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    map(filter, tissue == 'tumor') %>%
    map(select, -tissue) %>%
    map(column_to_rownames, 'sample_id')

# ssGSEA score calculation --------

  insert_msg('ssGSEA scores')

  recon$scores <- exp_data %>%
    future_map(calculate.default,
               x = set_names(recon$lexicon$gene_symbol,
                             recon$lexicon$variable),
               .options = furrr_options(seed = TRUE))

  ## appending with the sample IDs

  recon$scores <-
    map2(recon$scores,
         exp_data,
         ~mutate(.x,
                 sample_id = rownames(.y))) %>%
    map(relocate, sample_id)

  ## updating the lexicon

  recon$lexicon <- recon$lexicon %>%
    filter(variable %in% reduce(map(recon$scores, names), intersect))

# Caching --------

  insert_msg('Caching')

  save(recon, file = './data/recon.RData')

# END -------

  rm(exp_data)

  plan('sequential')

  insert_tail()
