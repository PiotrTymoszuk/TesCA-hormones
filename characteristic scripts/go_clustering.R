# Semantic clustering of the common regulated GOs

  insert_head()

# container ----

  bcg_goclust <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  ## common enriched GOs and their identifiers

  bcg_goclust$go_lexicon <- bcg_go$test[[1]] %>%
    map_dfr(select, go_id, term) %>%
    filter(!duplicated(go_id))

  bcg_goclust$go_ids <- bcg_go$common_significant %>%
    exchange(bcg_goclust$go_lexicon,
             key = 'term',
             value = 'go_id') %>%
    compact %>%
    map(unname)

  ## GOdata object

  bcg_goclust$godata <- godata(OrgDb = org.Hs.eg.db,
                          ont = 'BP')

# Wang distance calculation and MDS objects -------

  insert_msg('Wang distances')

  bcg_goclust$distances <- bcg_goclust$go_ids %>%
    map(go_sem,
        semData = bcg_goclust$godata,
        .parallel = TRUE)

  ## two-dimensional MDS

  bcg_goclust$component_tbl <- bcg_goclust$distances %>%
    map(cmdscale, k = 2) %>%
    map(as.data.frame) %>%
    map(set_names, c('comp_1', 'comp_2')) %>%
    map(rownames_to_column, 'observation')

  bcg_goclust$mds_objects <- bcg_goclust$component_tbl %>%
    map(~list(red_obj = .x,
              red_fun = 'mds',
              dist_method = 'custom',
              component_tbl = .x,
              loadings = NULL,
              data = quo(!!.x))) %>%
    map(red_analysis)

# PAM clustering of the MDS components -------

  insert_msg('PAM clustering of the MDS components')

  bcg_goclust$clust_objects <-
    list(data = bcg_goclust$mds_objects,
         k = c(3, 4, 3)) %>%
    pmap(kcluster,
         clust_fun = 'pam',
         distance_method = 'euclidean')

  ## cluster assignment and GO term names

  bcg_goclust$assignment <- bcg_goclust$clust_objects %>%
    map(extract, 'assignment') %>%
    map(mutate,
        term = exchange(observation,
                        bcg_goclust$go_lexicon,
                        key = 'go_id',
                        value = 'term'))

# Caching the results -------

  insert_msg('Caching the results')

  bcg_goclust <-
    bcg_goclust[c('distances', 'mds_objects',
                  'clust_objects', 'assignment')]

  save(bcg_goclust, file = './cache/bcg_goclust.RData')

# END -----

  insert_tail()
