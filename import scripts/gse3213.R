# Import of data sets for the GSE3218 study

  insert_head()

# container -------

  gse3218 <- list()

# raw data sets -------

  insert_msg('Raw data sets')

  gse3218$raw <- getGEO('GSE3218', destdir = './data/GSE3218')

# Clinical information -------

  insert_msg('Clinical information')

  ## both clinical data frames contain the sample clinical information

  gse3218$clinic <- gse3218$raw %>%
    map(pData) %>%
    map(as_tibble)

  gse3218$clinic[[1]] <- gse3218$clinic[[1]] %>%
    transmute(sample_id = stri_extract(title, regex = '\\d+.*$'),
              sample_id = stri_replace(sample_id, regex = '\\s{1}.*',
                                       replacement = ''),
              patient_id = stri_extract(sample_id, regex = '\\d+'),
              geo_accession1 = geo_accession,
              histology_details = characteristics_ch1,
              tissue = ifelse(stri_detect(histology_details, fixed = 'normal'),
                              'normal', 'tumor'),
              tissue = factor(tissue, c('normal', 'tumor')),
              histology = ifelse(stri_detect(histology_details,
                                             regex = '\\((s|S)eminoma\\)'),
                                 'seminoma', 'NSGCT'),
              histology = factor(histology, c('seminoma', 'NSGCT')),
              histology_icd = ifelse(stri_detect(histology_details,
                                                 regex = '^mixed'),
                                     'germinal mixed histology',
                                     stri_extract(histology_details,
                                                  regex = '\\(.*\\)')),
              histology_icd = stri_replace_all(histology_icd,
                                               regex = '\\(|\\)',
                                               replacement = ''),
              histology_icd = car::recode(tolower(histology_icd),
                                          "'anaplastic seminoma' = 'seminoma';
                                          'teratoma' = 'benign teratoma';
                                          'teratoma with secondary somatic malignancy' = 'malignant teratoma';
                                          'yolk sac tumor' = 'yolk sac cancer'"),
              histology_icd = factor(histology_icd),
              histology_icd = fct_relevel(histology_icd,
                                          'seminoma',
                                          'germinal mixed histology'))

  gse3218$clinic[[2]] <- gse3218$clinic[[2]] %>%
    transmute(sample_id = stri_extract(title, regex = '\\d+.*$'),
              sample_id = stri_replace(sample_id, regex = '\\s{1}.*',
                                       replacement = ''),
              patient_id = stri_extract(sample_id, regex = '\\d+'),
              geo_accession2 = geo_accession)

  ## merging, setting histology levels for the normal tissue

  gse3218$clinic <- gse3218$clinic %>%
    reduce(left_join, by = c('sample_id', 'patient_id')) %>%
    relocate(sample_id, patient_id, geo_accession1, geo_accession2) %>%
    mutate(histology = ifelse(tissue == 'normal',
                              'normal', as.character(histology)),
           histology = factor(histology, c('normal', 'seminoma', 'NSGCT')),
           histology_icd = ifelse(tissue == 'normal',
                                  'normal', as.character(histology_icd)),
           histology_icd = factor(histology_icd),
           histology_icd = fct_relevel(histology_icd,
                                       'normal',
                                       'seminoma',
                                       'germinal mixed histology'))

# Annotation --------

  insert_msg('Annotation')

  gse3218$annotation <- gse3218$raw %>%
    map(fData) %>%
    map(as_tibble) %>%
    map(transmute,
        probe_id = ID,
        entrez_id = stri_extract(ENTREZ_GENE_ID, regex = '\\d+')) %>%
    map(~filter(.x, complete.cases(.x))) %>%
    map(mutate,
        gene_symbol = mapIds(org.Hs.eg.db,
                             keys = entrez_id,
                             keytype = 'ENTREZID',
                             column = 'SYMBOL')) %>%
    map(~filter(.x, complete.cases(.x)))

  ## a common annotation frame

  gse3218$annotation <- gse3218$annotation %>%
    reduce(rbind) %>%
    filter(!duplicated(probe_id))

# Expression -------

  insert_msg('Expression')

  gse3218$expression <- gse3218$raw %>%
    map(exprs) %>%
    map(~.x[rownames(.x) %in% gse3218$annotation$probe_id, ])

  ## setting the study-provided columns names

  for(i in seq_along(gse3218$expression)) {

    colnames(gse3218$expression[[i]]) <-
      colnames(gse3218$expression[[i]]) %>%
      exchange(gse3218$clinic,
               key = paste0('geo_accession', i),
               value = 'sample_id') %>%
      unname

  }

  ## merging and collapsing expression of genes targeted by multiple samples
  ## by arithmetic mean

  gse3218$expression <- gse3218$expression %>%
    reduce(rbind) %>%
    integrate_expression(annotation = gse3218$annotation)

  ## appending the expression data frame with tissue type

  gse3218$expression <-
    inner_join(gse3218$clinic[c('sample_id', 'tissue')],
               gse3218$expression,
               by = 'sample_id')

  gse3218$annotation <- gse3218$annotation %>%
    filter(!duplicated(gene_symbol))

# Caching the cleared data sets -------

  insert_msg('Caching the cleared data sets')

  gse3218$raw <- NULL

  gse3218 <- compact(gse3218)

  save(gse3218, file = './data/gse3218.RData')

# END -----

  insert_tail()
