# Tables for the manuscript

  insert_head()

# container ------

  tabs <- list()

# cohort characteristic --------

  insert_msg('Cohort characteristic')

  ## only for the analyzed cohorts

  tabs$cohort <- expl_cohorts$result_tbl %>%
    filter(!Variable %in% c('TSS',
                            'RFS',
                            'Progression-free survival',
                            'Tumor-specific death')) %>%
    mutate(Variable = stri_replace(Variable,
                                   fixed = 'OS',
                                   replacement = 'Follow-up, months')) %>%
    mdtable(label = 'cohort_characteristic',
            ref_name = 'cohorts',
            caption = paste('Characteristic of the investigated cohorts of',
                            'testicular cancer patients.',
                            'Numeric variables are presented as medians with',
                            'interquartile ranges and ranges.',
                            'Qualitative variables are shown as percentages',
                            'and counts of the categories within the complete',
                            'observation set.'))

# Genes of interest -------

  insert_msg('Genes of interest')

  tabs$genes <- globals$gene_lexicon %>%
    filter(gene_symbol %in% expl_dist$top_variables) %>%
    mutate(entrez_id = mapIds(org.Hs.eg.db,
                              keys = gene_symbol,
                              keytype = 'SYMBOL',
                              column = 'ENTREZID')) %>%
    select(class, gene_symbol, entrez_id) %>%
    set_names(c('Gene classification',
                'Gene symbol',
                'Entrez ID')) %>%
    mdtable(label = 'genes_of_interest',
            ref_name = 'genes',
            caption = 'Sex hormone-related genes of interest.')

# Saving on the disc ------

  insert_msg('Saving on the disc')

  tabs %>%
    save_excel(path = './paper/tables.xlsx')

# END -----

  insert_tail()
