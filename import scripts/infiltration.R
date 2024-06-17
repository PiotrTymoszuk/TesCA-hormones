# Infiltration estimates obtained by immunedeconvolution with the
# QuanTIseq, xCell, and MCPcounter algorithms

  insert_head()

# containers -------

  quantiseq <- list()
  xcell <- list()
  mcp <- list()

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# Input expression data -------

  insert_msg('Input data')

  exp_data <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    map(filter, tissue == 'tumor') %>%
    map(select, -tissue) %>%
    map(column_to_rownames, 'sample_id') %>%
    map(t)

  exp_data <- exp_data %>%
    map(~2^.x)

# Immunedeconvolution ----------

  insert_msg('Immunedeconvolution')

  quantiseq$infiltration <-
    list(gene_expression = exp_data,
         arrays = c(FALSE, TRUE, TRUE)) %>%
    future_pmap(deconvolute,
                method = 'quantiseq',
                .options = furrr_options(seed = TRUE))

  xcell$infiltration <-
    list(gene_expression = exp_data,
         arrays = c(FALSE, TRUE, TRUE)) %>%
    future_pmap(deconvolute,
                method = 'xcell',
                .options = furrr_options(seed = TRUE))

  mcp$infiltration <-
    list(gene_expression = exp_data,
         arrays = c(FALSE, TRUE, TRUE)) %>%
    future_pmap(deconvolute,
                method = 'mcp_counter',
                .options = furrr_options(seed = TRUE))

# Cell type lexicons --------

  insert_msg('Cell type lexicons')

  quantiseq$lexicon <-
    c('B', 'TAM M1', 'TAM2',
      'Mono', 'Neutro', 'NK',
      'T CD4+', 'T CD8+',
      'Treg', 'mDC', 'other') %>%
    set_names(quantiseq$infiltration[[1]]$cell_type) %>%
    compress(names_to = 'variable',
             values_to = 'label')

  xcell$lexicon <-
    c('act mDC',
      'B',
      'T CD4+ memory',
      'T CD4+ naive',
      'T CD4+',
      'T CD4+ cm',
      'T CD4+ em',
      'T CD8+ naive',
      'T CD8+',
      'T CD8+ cm',
      'T CD8+ em',
      'B switch',
      'CLP', 'CMP', 'mDC',
      'EC', 'Eosino', 'CAF',
      'GMP', 'HSC', 'TAM',
      'TAM M1', 'TAM M2',
      'Mast', 'B memory',
      'Mono', 'B naive',
      'Neutro', 'NK', 'NKT',
      'pDC', 'B plasma', 'T gamma/delta',
      'T CD4+ Th1', 'T CD4+ Th2',
      'Treg', 'immune score', 'stroma score', 'TME score') %>%
    set_names(xcell$infiltration[[1]]$cell_type) %>%
    compress(names_to = 'variable',
             values_to = 'label')

  mcp$lexicon <-
    c('T', 'T CD8+', 'CTX score',
      'NK', 'B', 'Mono',
      'TAM', 'mDC', 'Neutro',
      'EC', 'CAF') %>%
    set_names(mcp$infiltration[[1]]$cell_type) %>%
    compress(names_to = 'variable',
             values_to = 'label')

# Formatting the infiltration estimates -------

  insert_msg('Formatting the infiltration estimates')

  quantiseq$infiltration <- quantiseq$infiltration %>%
    map(column_to_rownames, 'cell_type') %>%
    map(t) %>%
    map(as.data.frame) %>%
    map(rownames_to_column, 'sample_id') %>%
    map(as_tibble)

  xcell$infiltration <- xcell$infiltration %>%
    map(column_to_rownames, 'cell_type') %>%
    map(t) %>%
    map(as.data.frame) %>%
    map(rownames_to_column, 'sample_id') %>%
    map(as_tibble)

  mcp$infiltration <- mcp$infiltration %>%
    map(column_to_rownames, 'cell_type') %>%
    map(t) %>%
    map(as.data.frame) %>%
    map(rownames_to_column, 'sample_id') %>%
    map(as_tibble)

# Caching the results ------

  insert_msg('Caching')

  save(quantiseq, file = './data/quantiseq.RData')
  save(xcell, file = './data/xcell.RData')
  save(mcp, file = './data/mcp.RData')

# END -----

  plan('sequential')

  rm(exp_data)

  insert_tail()
