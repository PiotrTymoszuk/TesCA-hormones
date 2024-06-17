# Differential gene expression for the hormonal clusters.
#
# Differentially regulated genes are identified by the following procedure:
#
# 1) one-way ANOVA with eta-square effect size metric
#
# 2) differences mean cluster expression and the cohort average investigated by
# one-sample T test
#
# 3) Differentially regulated genes are identified by pFDR(ANOVA) < 0.05,
# eta-square >= 0.14, pFDR(T-test) < 0.05 significant regulation of expression
# versus cohort mean.

  insert_head()

# container -----

  bcg_dge <- list()

# analysis variables and data ------

  insert_msg('Analysis variables and data')

  bcg_dge$annotation <- globals$cohort_expr %>%
    eval %>%
    map(~.x$annotation) %>%
    map(~filter(.x,
                complete.cases(.x),
                !duplicated(gene_symbol)))

  bcg_dge$data <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression)

  bcg_dge$data <-
    map2(bcg_globals$assignment[names(bcg_dge$data)],
         bcg_dge$data,
         left_join, by = 'sample_id')

# Descriptive stats -------

  insert_msg('Descriptive stats')

  bcg_dge$stats <-
    map2(bcg_dge$data,
         bcg_dge$annotation,
         ~select(.x, clust_id, all_of(.y$gene_symbol))) %>%
    map(fast_num_stats,
        split_factor = 'clust_id')

# ANOVA and significant effects ------

  insert_msg('ANOVA')

  bcg_dge$anova <-
    list(data = bcg_dge$data,
         variables = map(bcg_dge$annotation, ~.x$gene_symbol)) %>%
    pmap(test_anova,
         split_fct = 'clust_id',
         adj_method = 'BH',
         .parallel = TRUE) %>%
    map(~.x$anova)

  ## annotation of the ANOVA results with Entrez ID, identification
  ## of the significant effects

  bcg_dge$anova <-
    map2(bcg_dge$anova,
         bcg_dge$annotation,
         ~mutate(.x,
                 gene_symbol = response,
                 entrez_id = exchange(gene_symbol,
                                      .y,
                                      key = 'gene_symbol',
                                      value = 'entrez_id'),
                 eff_size = paste('\u03B7\u00B2 =', signif(effect_size)))) %>%
    ## for nicely formatted significance text
    map(re_adjust, method = 'none')

  bcg_dge$anova_significant <- bcg_dge$anova %>%
    map(filter,
        p_adjusted < 0.05,
        effect_size >= 0.14) %>%
    map(~.x$gene_symbol)

# Mean deviations of expression in the hormonal clusters -------

  insert_msg('Expression regulation as compared with grand means')

  bcg_dge$dev_test <-
    list(data = bcg_dge$data,
         variables = map(bcg_dge$annotation, ~.x$gene_symbol)) %>%
    pmap(avg_deviation,
         split_fct = 'clust_id',
         adj_method = 'BH')

  ## annotation of the results with Entrez ID

  bcg_dge$dev_test <-
    map2(bcg_dge$dev_test,
         bcg_dge$anova_significant,
         ~mutate(.x,
                 gene_symbol = variable,
                 anova_significant = ifelse(gene_symbol %in% .y,
                                            'yes', 'no'))) %>%
    map2(.,
         bcg_dge$annotation,
         ~mutate(.x,
                 entrez_id = exchange(gene_symbol,
                                      .y,
                                      key = 'gene_symbol',
                                      value = 'entrez_id')))

  ## regulation sign

  bcg_dge$dev_test <- bcg_dge$dev_test %>%
    map(mutate,
        ## effect size of regulation: delta of means divided by SD,
        ## i.e. Cohen's d
        effect_size = deviation_center/deviation_sd,
        ## regulation sign
        regulation = ifelse(p_adjusted >= 0.05 | anova_significant == 'no',
                            'ns',
                            ifelse(deviation_center > 0,
                                   'upregulated',
                                   ifelse(deviation_center < 0,
                                          'downregulated', 'ns'))),
        regulation = factor(regulation,
                            c('upregulated', 'downregulated', 'ns')))

# Identification of differentially regulated genes -------

  insert_msg('Differentially regulated genes')

  ## in single cohorts

  bcg_dge$significant <- bcg_dge$dev_test %>%
    map(filter, regulation %in% c('upregulated', 'downregulated')) %>%
    map(blast, clust_id, regulation) %>%
    transpose %>%
    map(map, ~.x$gene_symbol)

  ## common significant, shared by the TCGA and GSE99420 cohorts

  bcg_dge$common_significant <- bcg_dge$significant %>%
    map(~.x[c('tcga', 'gse99420')]) %>%
    map(reduce, intersect)

# Numbers of regulated genes -------

  insert_msg('Numbers of regulated genes')

  bcg_dge$dge_numbers$total <- bcg_dge$annotation %>%
    map_dbl(nrow) %>%
    compress(names_to = 'cohort',
             values_to = 'n_total')

  ## genes regulated in particular clusters

  bcg_dge$dge_numbers$clusters <- bcg_dge$significant %>%
    map(map_dbl, length) %>%
    map(compress,
        names_to = 'cohort',
        values_to = 'n') %>%
    compress(names_to = 'status')

  ## common regulated genes

  bcg_dge$dge_numbers$common <- bcg_dge$common_significant %>%
    map_dbl(length) %>%
    compress(names_to = 'status',
             values_to = 'n') %>%
    mutate(cohort = 'common')

  ## merging into one data frame

  bcg_dge$dge_numbers <-
    bcg_dge$dge_numbers[c("clusters", "common")] %>%
    reduce(rbind) %>%
    left_join(bcg_dge$dge_numbers$total,
              by = 'cohort') %>%
    mutate(clust_id = stri_split_fixed(status,
                                       pattern = '.',
                                       simplify = TRUE)[, 1],
           regulation = stri_split_fixed(status,
                                         pattern = '.',
                                         simplify = TRUE)[, 2],
           clust_id = factor(clust_id, levels(bcg_dge$data[[1]]$clust_id)),
           regulation = factor(regulation,
                               c('upregulated', 'downregulated')))

# Caching the results ------

  insert_msg('Caching the results')

  bcg_dge <- bcg_dge[c("stats",
                       "anova", "anova_significant", "dev_test",
                       "significant", "common_significant", "dge_numbers")]

  save(bcg_dge, file = './cache/bcg_dge.RData')

# END ----

  insert_tail()
