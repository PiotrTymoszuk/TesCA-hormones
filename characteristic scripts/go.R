# GO enrichment analysis for genes found to be differentially regulated in the
# hormonal clusters as compared with the cohort average expression.
#
# The enrichment in biological process terms is investigated by GOANA with OR
# as an effect size statistic for enrichment in the gene set as compared with
# the genomic occurrence.
#
# Significant GOs are identified by pFDR < 0.05.

  insert_head()

# container ------

  bcg_go <- list()

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis data --------

  insert_msg('Analysis data')

  ## universe genes: all genes detected by the respective RNAseq or microarray
  ## platform

  bcg_go$universe <- bcg_dge$dev_test %>%
    map(~.x$entrez_id) %>%
    map(unique)

  ## differentially regulated genes: collapsing up- and downregulated
  ## features in a particular cluster

  bcg_go$dge <-
    list(`#1` = bcg_dge$significant[c("#1.upregulated",
                                      "#1.downregulated")],
         `#2` = bcg_dge$significant[c("#2.upregulated",
                                      "#2.downregulated")],
         `#3` = bcg_dge$significant[c("#3.upregulated",
                                      "#3.downregulated")],
         `#4` = bcg_dge$significant[c("#4.upregulated",
                                      "#4.downregulated")]) %>%
    map(transpose) %>%
    map(map, reduce, c)


  bcg_go$dge <- bcg_go$dge %>%
    map(map, ~mapIds(org.Hs.eg.db,
                keys = .x,
                keytype = 'SYMBOL',
                column = 'ENTREZID')) %>%
    map(map, ~.x[!is.na(.x)]) %>%
    map(map, unname)

# Enrichment analysis ---------

  insert_msg('Enrichment analysis')

  for(i in names(bcg_go$dge)) {

    bcg_go$test[[i]] <-
      list(de = bcg_go$dge[[i]],
           universe = bcg_go$universe) %>%
      future_pmap(GOana,
                  ontology = 'BP',
                  adj_method = 'BH')

  }

# significantly regulated GO terms --------

  insert_msg('Significant entichment')

  ## in single cohorts

  bcg_go$significant <- bcg_go$test %>%
    map(map,
        filter,
        p_adjusted < 0.05) %>%
    map(map, ~.x$term)

  ## common significant: shared by the TCGA and GSE99420 cohorts

  bcg_go$common_significant <- bcg_go$significant %>%
    map(~.x[c('tcga', 'gse99420')]) %>%
    map(reduce, intersect)

# Numbers of significantly regulated GOs -------

  insert_msg('GO numbers')

  ## total numbers of analyzed GOs

  bcg_go$go_numbers$total <- bcg_go$test %>%
    map(map_dbl, nrow) %>%
    map(compress,
        names_to = 'cohort',
        values_to = 'n_total') %>%
    compress(names_to = 'clust_id')

  ## numbers of significantly enriched GOs

  bcg_go$go_numbers$clusters <- bcg_go$significant %>%
    map(map_dbl, length) %>%
    map(compress,
        names_to = 'cohort',
        values_to = 'n') %>%
    compress(names_to = 'clust_id')

  ## common

  bcg_go$go_numbers$common <- bcg_go$common_significant %>%
    map_dbl(length) %>%
    compress(names_to = 'clust_id',
             values_to = 'n') %>%
    mutate(cohort = 'common')

  ## merging into one table

  bcg_go$go_numbers <-
    rbind(bcg_go$go_numbers$clusters,
          bcg_go$go_numbers$common) %>%
    left_join(bcg_go$go_numbers$total,
              by = c('clust_id', 'cohort')) %>%
    mutate(clust_id = factor(clust_id,
                             levels(bcg_globals$assignment[[1]]$clust_id)))

# Caching the results --------

  insert_msg('Caching the results')

  bcg_go <-
    bcg_go[c("test", "significant", "common_significant", "go_numbers")]

  save(bcg_go, file = './cache/bcg_go.RData')

# END ------

  rm(i)

  plan('sequential')

  insert_tail()
