# Analysis of targets of drugs differing in sensitivity between the hormonal
# clusters. Targets of the significant drugs in the TCGA and GSE99420 cohorts
# are analyzed.

  insert_head()

# container ------

  bcg_targets <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  ## significantly regulated drugs and drug lexicons

  bcg_targets$drugs <- bcg_drugs$significant %>%
    map(map, map, ~.x[c('tcga', 'gse99420')])

  ## lexicons: extracting additional targets from the
  ## MOA (mechanism of action) and pathway columns

  bcg_targets$regex_lexicon <-
    c('DNA synthesis' = 'DNA\\s{1}',
      'DNA synthesis' = 'DNA$',
      'DNA synthesis' = '(p|P)yrimidine',
      'DNA synthesis' = '(p|P)urine',
      'DNA synthesis' = '(t|T)opoisomerase',
      'DNA synthesis' = 'Antimetabolite',
      'DNA synthesis' = 'DHFR|TYMS|WRN|FEN',
      'DNA synthesis' = 'TOP\\d{1}',
      'DNA synthesis' = '(T|t)elomerase',
      'ROS' = 'ROS',
      'Apoptosis' = 'MCL1|CASP|BIRC|DIABLO|XIAP|SMAC',
      'Apoptosis' = 'Bcl|BCL|bcl|Bax|BAX|PARP|IAP',
      'Apoptosis' = '(A|a)poptosis',
      'Cell cycle' = '(A|a)urora',
      'Cell cycle' = 'PLK1|AURK|CDK',
      'Cell cycle' = '(C|c)ell\\s+cycle',
      'Cell cycle' = '(C|c)yclin',
      'Cell cycle' = '(M|m)itosis',
      'DNA repair' = 'p53|TP53|MDM|WRN|CHEK|ATM|ATR|GADD34|MRE11|PPM1D',
      'Microtubule' = 'tubule|centrosome',
      'Cytoskeleton' = '(c|C)ytoskeleton',
      'Cytoskeleton' = '(m|M)icrotubule',
      'Translation' = '(t|T)ranslation',
      'Translation' = 'EF\\s+',
      'Translation' = '(e|E)IF\\s+',
      'Respiration' = '(r|R)espiration',
      'Respiration' = 'ATP\\s{1}synthase',
      'Chromatin' = '(H|h)istone',
      'Chromatin' = '(C|c)hromatin',
      'Chromatin' = 'HDAC',
      'WNT' = 'WNT',
      'Angiogenesis' = 'KDR|VEGF|PDGFR|Eph|EPH|Ephrin|ephrin|FLT',
      'RTK signaling' = 'VEGFR|KDR|FLT|KIT|PDGFR|MET|RET|FGFR|TIE',
      'ERBB signaling' = 'EGFR|ERBB\\d{1}|HER\\d{1}')

  bcg_targets$lexicons <-
    drugs$lexicons[names(bcg_targets$drugs)] %>%
    map(target_extractor,
        target_column = 'targets',
        pathway_column = 'targets',
        regex_lexicon = bcg_targets$regex_lexicon) %>%
    map2(c('moa', 'pathway_name'),
         ~target_extractor(.x,
                           pathway_column = .y,
                           regex_lexicon = bcg_targets$regex_lexicon))

# Total numbers of drugs per target ------

  insert_msg('Total number of drugs per molecular targets')

  bcg_targets$n_totals <- bcg_targets$lexicons %>%
    map(~.x$targets) %>%
    map(unlist) %>%
    map(table) %>%
    map(compress,
        names_to = 'target',
        values_to = 'n_total') %>%
    map(mutate,
        n_total = as.numeric(n_total))

# Identification of the drug targets ------

  insert_msg('Identification of the drug targets')

  for(i in names(bcg_targets$drugs)) {

    for(j in names(bcg_targets$drugs[[i]])) {

      bcg_targets$stats[[i]][[j]] <-
        bcg_targets$drugs[[i]][[j]] %>%
        map(map,
            ~filter(bcg_targets$lexicons[[i]],
                    variable %in% .x)) %>%
        map(map,  ~compact(.x$targets)) %>%
        map(map, unlist) %>%
        map(map, table) %>%
        map(map,
            compress,
            names_to = 'target',
            values_to = 'n') %>%
        map(compress,
            names_to = 'cohort') %>%
        compress(names_to = 'status')

    }

    bcg_targets$stats[[i]] <- bcg_targets$stats[[i]] %>%
      compress(names_to = 'clust_id') %>%
      mutate(clust_id = factor(clust_id,
                               levels(bcg_globals$assignment[[1]]$clust_id)),
             n = as.numeric(n),
             status = factor(status, c('resistant', 'sensitive')),
             cohort = factor(cohort, c('tcga', 'gse99420')))

    ## appending with the total drug numbers per target
    ## Fisher's exact tests for enrichment

    bcg_targets$stats[[i]] <-
      left_join(bcg_targets$stats[[i]],
                bcg_targets$n_totals[[i]],
                by = 'target') %>%
      test_targets

    ## top most frequent targets shared by both cohorts and
    ## hit by more than two drugs

    bcg_targets$common_targets[[i]] <- bcg_targets$stats[[i]] %>%
      blast(clust_id) %>%
      map(blast, status) %>%
      map(map, blast, cohort) %>%
      map(map, map, filter, n > 2) %>%
      map(map, map, ~.x$target) %>%
      map(map, reduce, intersect)

    bcg_targets$top_targets[[i]] <-
      bcg_targets$stats[[i]] %>%
      blast(clust_id) %>%
      map(blast, status)

    for(j in names(bcg_targets$top_targets[[i]])) {

      bcg_targets$top_targets[[i]][[j]] <-
        map2_dfr(bcg_targets$top_targets[[i]][[j]],
                 bcg_targets$common_targets[[i]][[j]],
                 ~filter(.x, target %in% .y))

    }

    bcg_targets$top_targets[[i]] <-
      bcg_targets$top_targets[[i]] %>%
      reduce(rbind)

  }


# Bar plots for the top most frequent targets ------

  insert_msg('Bar plots for the top most frequent targets')

  ## faceted bar plots

  bcg_targets$top_plots <-
    list(x = bcg_targets$top_targets %>%
           map(mutate,
               plot_n = ifelse(status == 'sensitive',
                               n, -n)),
         y = paste0('Top drug targets, ',
                   globals$drug_exp_labs,
                   '-trained drug sensitivity predictions')) %>%
    pmap(function(x, y) x %>%
           ggplot(aes(x = n,
                      y = reorder(target, n),
                      fill = status)) +
           geom_bar(stat = 'identity',
                    color = 'black') +
           facet_grid(clust_id + status ~ cohort,
                      scales = 'free',
                      space = 'free',
                      labeller = labeller(.cols = globals$cohort_labs)) +
           scale_fill_manual(values = c(resistant = 'firebrick',
                                        sensitive = 'steelblue'),
                             name = 'Predicted response\nvs cohort mean') +
           globals$common_theme +
           theme(axis.title.y = element_blank(),
                 strip.text.y = element_text(angle = 0,
                                             hjust = 0)) +
           labs(title = y,
                subtitle = 'Molecular targets of common significant drugs',
                x = '# of drugs'))

# END ------

  insert_tail()
