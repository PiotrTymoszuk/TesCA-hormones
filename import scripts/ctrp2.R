# Import of drug sensitivity and expression information for the CTRP2 data set
# of cell lines subjected to treatment with various anti-cancer drugs.
#
# This data will be used for prediction of drug sensitivity for testicular
# cancer samples.

  insert_head()

# container ------

  ctrp2 <- list()

# sensitivity data: AUC ------

  insert_msg('Sensitivity data')

  ctrp2$auc <-
    read_csv('./data/CTRP2/Drug_sensitivity_AUC_(CTD^2).csv')

  names(ctrp2$auc)[1] <- 'sample_id'

# drug lexicon -------

  insert_msg('Drug lexicon')

  ctrp2$drug_lexicon <-
    tibble(variable = names(ctrp2$auc)[-1]) %>%
    mutate(drug_id = stri_extract(variable,
                                  regex = '\\(.*\\)$'),
           drug_id = stri_replace_all(drug_id,
                                      regex = '\\(|\\)',
                                      replacement = ''),
           drug_id = stri_replace(drug_id,
                                  regex = '^CTRP:',
                                  replacement = ''))

  ## appending the lexicon with drug targets and MOA

  ctrp2$appendix <- read_tsv('./data/CTRP2/v20.meta.per_compound.txt')

  ctrp2$appendix <- ctrp2$appendix %>%
    transmute(drug_id = as.character(master_cpd_id),
              drug_name = cpd_name,
              label = drug_name,
              targets = stri_split_fixed(gene_symbol_of_protein_target,
                                         pattern = ';'),
              moa = target_or_activity_of_compound,
              broad_cpd_id = broad_cpd_id)

  ctrp2$drug_lexicon <-
    left_join(ctrp2$drug_lexicon,
              ctrp2$appendix,
              by = 'drug_id')

# Cell line lexicon ------

  insert_msg('Cell line details')

  ctrp2$cell_lexicon <-
    read_csv('./data/CTRP2/Model.csv')

  ctrp2$cell_lexicon <- ctrp2$cell_lexicon %>%
    transmute(sample_id = ModelID,
              cell_line = StrippedCellLineName,
              tcga_entity = OncotreeCode,
              tissue = OncotreeLineage,
              tissue_subtype = OncotreeSubtype) %>%
    filter(sample_id %in% ctrp2$auc$sample_id)

# reading the baseline expression data -------

  insert_msg('Reading the baseline expression data')

  ctrp2$expression <-
    read_csv('./data/CTRP2/Expression_Public_23Q4.csv')

  names(ctrp2$expression)[1] <- 'sample_id'

# Annotation -------

  insert_msg('Annotation')

  ctrp2$annotation <-
    tibble(probe_id = names(ctrp2$expression)[-1]) %>%
    mutate(entrez_id = mapIds(org.Hs.eg.db,
                              keys = probe_id,
                              keytype = 'SYMBOL',
                              column = 'ENTREZID'),
           entrez_id = ifelse(!is.na(entrez_id),
                              entrez_id,
                              mapIds(org.Hs.eg.db,
                                     keys = probe_id,
                                     keytype = 'ALIAS',
                                     column = 'ENTREZID'))) %>%
    filter(complete.cases(.)) %>%
    mutate(gene_symbol = mapIds(org.Hs.eg.db,
                                keys = entrez_id,
                                keytype = 'ENTREZID',
                                column = 'SYMBOL')) %>%
    filter(complete.cases(.))

  ## re-naming the expression data set: official symbols

  ctrp2$expression <- ctrp2$expression %>%
    select(sample_id, all_of(ctrp2$annotation$probe_id)) %>%
    set_names(c('sample_id', ctrp2$annotation$gene_symbol))

# Removal of duplicated features --------

  insert_msg('Unique genes only')

  ctrp2$unique_genes <- names(ctrp2$expression)[-1] %>%
    unique

  ctrp2$expression <- ctrp2$expression %>%
    select(sample_id, all_of(ctrp2$unique_genes))

# Caching the formatted data -------

  insert_msg('Caching')

  ctrp2 <-
    ctrp2[c("auc",
            "drug_lexicon", "cell_lexicon",
            "expression", "annotation")]

  save(ctrp2, file = './data/ctrp2.RData')

# END ------

  insert_tail()
