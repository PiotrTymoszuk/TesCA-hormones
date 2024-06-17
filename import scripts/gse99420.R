# Import of data sets for the GSE99420 study

  insert_head()

# container -------

  gse99420 <- list()

# raw data sets -------

  insert_msg('Raw data sets')

  gse99420$raw <- getGEO('GSE99420', destdir = './data/GSE99420')

# Clinical information ------

  insert_msg('Clinical information')

  gse99420$clinic <- gse99420$raw[[1]] %>%
    pData %>%
    as_tibble

  gse99420$clinic <- gse99420$clinic %>%
    transmute(sample_id = geo_accession,
              tissue = factor('tumor', c('normal', 'tumor')),
              histology = ifelse(`non-seminoma (ns) and seminoma (s):ch1` == 'S',
                                 'seminoma', 'NSGCT'),
              histology = factor(histology, c('seminoma', 'NSGCT')),
              relapse = ifelse(`relapsed (r) vs. non-relapsed (nr):ch1` == 'R',
                               1, 0),
              relapse_factor = ifelse(`relapsed (r) vs. non-relapsed (nr):ch1` == 'R',
                                      'yes', 'no'),
              relapse_factor = factor(relapse_factor, c('no', 'yes')))

# Annotation ------

  insert_msg('Annotation')

  gse99420$annotation <- gse99420$raw[[1]] %>%
    fData %>%
    as_tibble

  gse99420$annotation <- gse99420$annotation %>%
    transmute(probe_id = ID,
              entrez_id = as.character(Entrez_Gene_ID)) %>%
    filter(complete.cases(.)) %>%
    mutate(gene_symbol = mapIds(org.Hs.eg.db,
                                keys = entrez_id,
                                keytype = 'ENTREZID',
                                column = 'SYMBOL')) %>%
    filter(complete.cases(.))

# Expression -----

  insert_msg('Expression')

  gse99420$expression <- gse99420$raw[[1]] %>%
    exprs %>%
    integrate_expression(annotation = gse99420$annotation) %>%
    mutate(tissue = factor('tumor', c('normal', 'tumor')))

  gse99420$annotation <- gse99420$annotation %>%
    filter(!duplicated(gene_symbol))

# Caching the cleared data sets ------

  insert_msg('Caching the cleared data sets')

  gse99420$raw <- NULL

  gse99420 <- compact(gse99420)

  save(gse99420, file = './data/gse99420.RData')

# END -----

  insert_tail()



