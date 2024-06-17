# Import of drug sensitivity and expression information for the GDSC data set
# of cell lines subjected to treatment with various anti-cancer drugs.
#
# This data will be used for prediction of drug sensitivity for testicular
# cancer samples.
#
# Data from the GDSC1 and GDSC2 experiments are formatted jointly. In case of
# drugs investigated in both experiments, the GDSC2 results are used for
# predictions of drug response.

  insert_head()

# container -------

  gdsc <- list()

# Drug lexicon ------

  insert_msg('Drug lexicon')

  gdsc$drug_lexicon <-
    c(gdsc1 = './data/GDSC/compounds_gdsc1.tab',
      gdsc2 = './data/GDSC/compounds_gdsc2.tab') %>%
    map(read_tsv) %>%
    map(mutate,
        variable = paste(drug_name, drug_id, sep = '_'),
        targets = ifelse(targets == '-', NA, targets),
        pathway_name = ifelse(pathway_name == '-', NA, pathway_name),
        targets = stri_split_fixed(targets,
                                   pattern = ', ',
                                   simplify = FALSE))

  gdsc$drug_lexicon$gdsc1 <- gdsc$drug_lexicon$gdsc1 %>%
    filter(!drug_name %in% gdsc$drug_lexicon$gdsc2$drug_name)

# Experiment design data --------

  insert_msg('Experiment design and IC50')

  gdsc$design <- c('gdsc1' = './data/GDSC/GDSC1_IC50.csv',
                   'gdsc2' = './data/GDSC/GDSC2_IC50.csv') %>%
    map(read_csv)

  ## clearing

  gdsc$design <- gdsc$design %>%
    map(transmute,
        variable = paste(`Drug Name`, `Drug ID`, sep = '_'),
        drug_id = `Drug ID`,
        drug_name = `Drug Name`,
        sample_id = paste0('cosmic_', `Cosmic ID`),
        cosmic_id = `Cosmic ID`,
        cell_line = `Cell Line Name`,
        tcga_entity = `TCGA Classification`,
        tissue = Tissue,
        tissue_subtype = `Tissue Sub-type`,
        log_ic50 = IC50,
        auc = AUC,
        rmse = RMSE)

  ## restricting to the drugs specified by the lexicon

  gdsc$design <- gdsc$design %>%
    map2(., gdsc$drug_lexicon,
         ~filter(.x, variable %in% .y$variable))

# logIC50 in wide format ---------

  insert_msg('logIC50 in wide format')

  gdsc$log_ic50 <- gdsc$design %>%
    map(select, sample_id, variable, log_ic50) %>%
    map(pivot_wider,
        id_cols = 'sample_id',
        names_from = 'variable',
        values_from = 'log_ic50')

# Cell line lexicon -------

  insert_msg('Cell line lexicon')

  gdsc$cell_lexicon <- gdsc$design %>%
    map(select,
        sample_id,
        cosmic_id,
        cell_line,
        tcga_entity,
        tissue,
        tissue_subtype) %>%
    map(filter,
        !duplicated(sample_id))

# Reading the expression data -------

  insert_msg('Reading the gene expression data')

  gdsc$expression <-
    read_tsv('./data/GDSC/Cell_line_RMA_proc_basalExp.txt') %>%
    mutate(probe_id = paste0('probe_', 1:nrow(.)))

# Probe/gene annotation ------

  insert_msg('Probe/gene annotation')

  gdsc$annotation <- gdsc$expression %>%
    transmute(probe_id = probe_id,
              gene_symbol = GENE_SYMBOLS) %>%
    filter(complete.cases(.))

  ## Entrez ID: some of 'gene symbols' are in fact aliases

  gdsc$annotation <- gdsc$annotation %>%
    mutate(entrez_id = mapIds(org.Hs.eg.db,
                              keys = gene_symbol,
                              keytype = 'SYMBOL',
                              column = 'ENTREZID')) %>%
    mutate(entrez_id = ifelse(!is.na(entrez_id), entrez_id,
                              mapIds(org.Hs.eg.db,
                                     keys = gene_symbol,
                                     keytype = 'ALIAS',
                                     column = 'ENTREZID'))) %>%
    filter(complete.cases(.))

  ## official gene symbols

  gdsc$annotation <- gdsc$annotation %>%
    mutate(gene_symbol = mapIds(org.Hs.eg.db,
                                keys = entrez_id,
                                keytype = 'ENTREZID',
                                column = 'SYMBOL')) %>%
    filter(complete.cases(.))

# Clearing of the expression data set -----

  insert_msg('Clearing of the expression data set')

  gdsc$expression<- gdsc$expression %>%
    select(-GENE_SYMBOLS, -GENE_title) %>%
    column_to_rownames('probe_id') %>%
    integrate_expression(annotation = gdsc$annotation)

  ## sample ID: official COSMIC identifier

  gdsc$expression <- gdsc$expression %>%
    mutate(sample_id = stri_replace(sample_id,
                                    regex = '^DATA.',
                                    replacement = 'cosmic_'))

# Caching the formatted data sets -------

  insert_msg('Caching')

  gdsc <-
    gdsc[c("log_ic50",
           "drug_lexicon", "cell_lexicon",
           "expression", "annotation")]

  save(gdsc, file = './data/gdsc.RData')

# END -----

  insert_tail()
