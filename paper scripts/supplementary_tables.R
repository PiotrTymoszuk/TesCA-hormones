# Supplementary tables for the manuscript

  insert_head()

# container --------

  suppl_tabs <- list()

# Differential gene expression: main histologies -------

  insert_msg('Differential gene expression, SEM vs NGSCT')

  suppl_tabs$histology <- expl_histo$result_tbl[c("tcga", "gse99420")] %>%
    map(~set_colnames(.x, stri_capitalize(names(.x)))) %>%
    compress(names_to = 'Cohort') %>%
    mutate(Cohort = globals$cohort_labs[Cohort]) %>%
    relocate(Cohort) %>%
    mdtable(label = 'hormone_genes_histology',
            ref_name = 'histology',
            caption = paste('Differential expression of hormone-related genes',
                            'in seminoma and NGSCT in the TCGA and GSE99420',
                            'cohorts.',
                            'Log2-transformed expression levels are presented',
                            'as medians with interquartile ranges and ranges.',
                            'Significant effects are shown.',
                            'The full table is available as a supplementary',
                            'Excel table.'))

# Differential gene expression, histology details, TCGA -------

  insert_msg('Differential gene expression, histology details')

  suppl_tabs$histology_icd <- expl_icd$result_tbl$tcga %>%
    mdtable(label = 'hormone_genes_histology_icd',
            ref_name = 'histology_icd',
            caption = paste('Expression of hormone-related genes in testicular',
                            'cancers of the TCGA cohort stratified by ICD-O',
                            'histological subtypes.',
                            'Log2-transformed expression levels are presented',
                            'as medians with interquartile ranges and ranges.'))

# Network stats for the hormone-related genes -------

  insert_msg('Network stats for the hormone-related genes')

  suppl_tabs$hormone_nets <- expl_net$stats[c("tcga", "gse99420")] %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = globals$cohort_labs[cohort]) %>%
    select(cohort, variable, class, degree, betweenness, hub_score) %>%
    map_dfc(function(x) if(is.numeric(x)) signif(x, 2) else x) %>%
    set_names(c('Cohort', 'Gene symbol', 'Gen classification',
                'Degree', 'Betweenness', 'Hub score')) %>%
    mdtable(label = 'hormone_gene_network_stats',
            ref_name = 'hormone_nets',
            caption = paste('Statistics for co-expression networks of',
                            'hormone-related genes in',
                            'in the TCGA and GSE99420 cohorts.',
                            'Top genes with the largest degree, betweenness,',
                            'and hub scores in each cohort are presented.',
                            'The complete table is available as a',
                            'supplementary Excel file.'))

# Cluster performance stats -------

  insert_msg('Cluster performance stats')

  suppl_tabs$cluster_stats <- clust_eval$stats %>%
    filter(cohort %in% c('tcga', 'gse99420')) %>%
    mutate(cohort = globals$cohort_labs[cohort]) %>%
    select(cohort,
           sil_width,
           frac_misclassified,
           frac_var,
           neighborhood_error) %>%
    map_dfc(function(x) if(is.numeric(x)) signif(x, 2) else x) %>%
    set_names('Cohort',
              'Silhouette width',
              'Misclassification rate',
              'Explained clustering variance',
              'Neighborhood misclassification') %>%
    mdtable(label = 'clustering_stats',
            ref_name = 'cluster_stats',
            caption = paste('Metrics of cluster separation, potential',
                            'misclassification, explained variance, and',
                            'neighborhood misclassification for the hormonal',
                            'clusters in the TCGA training cohort and the',
                            'GSE99420 test collective.'))

# Expression of the cluster-defining hormone-related genes ------

  insert_msg('Expression of hormone-related cluster-defining genes')

  suppl_tabs$cluster_factors <- clust_ft$result_tbl %>%
    filter(Cohort %in% c('TCGA', 'GSE99420')) %>%
    mdtable(label = 'clustering_features',
            ref_name = 'cluster_factors',
            caption = paste('Expression of the hormone-related cluster-defining',
                            'genes in the hormonal clusters of testicular cancer',
                            'samples in the training TCGA cohort and the GSE99420',
                            'test collective.',
                            'Statistical significance of differences between the',
                            'clusters was determined by Kruskal-Wallis test with',
                            'with eta-square effect size statistic.',
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'Log2-transformed expression levels are shown as',
                            'medians with interquartile ranges and ranges.',
                            'The table is available as a supplementary',
                            'Excel file.'))

# Clinical characteristic of the clusters ------

  insert_msg('Clinical characteristic of the clusters')

  suppl_tabs$cluster_clinic <- bcg_clinic$result_tbl %>%
    filter(Cohort %in% c('TCGA', 'GSE99420')) %>%
    mdtable(label = 'clusters_clinics',
            ref_name = 'cluster_clinics',
            caption = paste('Clinical characteristic of the hormonal clusters',
                            'in the TCGA and GSE99420 cohorts.',
                            'Quantitative variables are presented as medians',
                            'with interquartile ranges and ranges.',
                            'Qualitative variables are presented as percentages',
                            'and counts of the categories within the clusters.'))

# Performance metrics for the RIDGE PFS Cox models ---------

  insert_msg('Cox model performance')

  suppl_tabs$pfs_cox <- bcg_glmcox$stats %>%
    transmute(`Model type` = bcg_glmcox$model_labs[model_type],
              `Progression cases, N` = n_events,
              `Observations, N` = n_complete,
              `C-index` = paste0(signif(c_index, 2),
                               ' [95% CI: ', signif(lower_ci, 2),
                               ' to ', signif(upper_ci, 2),
                               ']'),
              `R-square` = signif(raw_rsq, 2),
              `IBS` = signif(ibs_model, 2)) %>%
    mdtable(label = 'pfs_cox',
            ref_name = 'pfs_cox',
            caption = paste('Numeric statistics of performance of RIDGE',
                            'Cox proportional hazard models of progression-free',
                            'survival in the TCGA cohort.'))

# Clusters, infiltration -------

  insert_msg('Clusters, infiltration')

  suppl_tabs$infiltration <- bcg_infil$result_tbl %>%
    map(filter, Cohort %in% c('TCGA', 'GSE99420')) %>%
    compress(names_to = 'Algorithm') %>%
    relocate(Algorithm) %>%
    mdtable(label = 'clusters_infiltration',
            ref_name = 'infiltration',
            caption = paste('Non-malignant cell content in the hormonal',
                            'clusters of the TCGA and GSE99420 cohorts was',
                            'estimated by the QuanTIseq, xCell, and MCP Counter',
                            'algorithms.',
                            'Differences in the predicted cell levels between',
                            'the clusters were investigated',
                            'by Kruskal-Wallis test with eta-square effect size',
                            'statistic.',
                            'P values were corrected for multiple tesing with',
                            'the false discovery rate method.',
                            'Median infiltration levels with interquartile',
                            'ranges and ranges are shown.',
                            'The table is available as a supplementary Excel',
                            'file.'))

# Reactome pathways -------

  insert_msg('Reactome pathways')

  suppl_tabs$reactome <- bcg_reactome$result_tbl %>%
    filter(Cohort %in% c('TCGA', 'GSE99420')) %>%
    mdtable(label = 'clusters_reactome_pathways',
            ref_name = 'reactome',
            caption = paste('Differences in single sample gene set enrichment',
                            'analysis scores (ssGSEA scores) of the Reactome',
                            'pathway gene signatures were compared between the',
                            'hormonal clusters by one-way ANOVA with eta-square',
                            'effect size statistic.',
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'Median ssGSEA scores with interquartile ranges',
                            'and ranges are shown for gene signatures found',
                            'to differ significantly between the clusters',
                            'in both the TCGA and GSE99420 cohorts.',
                            'The table is available as a supplementary Excel',
                            'file.'))

# Differential gene expression -------

  insert_msg('Differential gene expression')

  suppl_tabs$dge <- bcg_dge$dev_test %>%
    map(filter,
        regulation %in% c('upregulated', 'downregulated')) %>%
    map(select,
        clust_id, regulation,
        gene_symbol, entrez_id,
        deviation_center, deviation_sem,
        lower_ci, upper_ci,
        p_adjusted) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = globals$cohort_labs[cohort]) %>%
    relocate(cohort) %>%
    set_names(c('Cohort',
                'Cluster', 'Regulation vs cohort mean',
                'Gene symbol', 'Entrez ID',
                'log2 fold-regulation vs cohort mean',
                'standard error',
                '95% CI, lower bound', '95% CI, upper bound',
                'pFDR')) %>%
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x)

  suppl_tabs$dge <- suppl_tabs$dge %>%
    mdtable(label = 'dge',
            ref_name = 'dge',
            caption = paste('Genes differentially regulated in the hormonal',
                            'clusters as compared with the cohort means.',
                            'Statistical significance of differences in',
                            'log2-transformed expression levels between the',
                            'clusters was assessed by one-way ANOVA with',
                            'eta-square effect size statistic.',
                            'Statistical significance of differences between',
                            'log2-transformed expression in the cluster and',
                            'the cohort mean was determined by one-sample',
                            "two-tailed T test",
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'Differential regulation in a particular cluster',
                            'was considered for',
                            'pFDR(ANOVA) < 0.05, eta-squared of at least 0.14,',
                            "and pFDR(T test) < 0.05.",
                            'The table is available as a supplementary',
                            'Excel file.'))

# Biological process GO enrichment -------

  insert_msg('GO enrichment in the clusters')

  suppl_tabs$go <- bcg_go$test %>%
    map(~.x[c('tcga', 'gse99420')]) %>%
    map(compress, names_to = 'cohort') %>%
    map2(.,
         bcg_go$common_significant,
         ~filter(.x, term %in% .y)) %>%
    compress(names_to = 'clust_id') %>%
    mutate(cohort = globals$cohort_labs[cohort]) %>%
    select(clust_id, cohort, term, go_id, or, p_adjusted) %>%
    map_dfc(function(x) if(is.numeric(x)) signif(x, 2) else x) %>%
    set_names(c('Cluster',
                'Cohort',
                'GO term name',
                'GO term ID',
                'Enrichment OR',
                'pFDR')) %>%
    mdtable(label = 'clusters_go_enrichment',
            ref_name = 'go',
            caption = paste('Biological process gene ontology (GO) term',
                            'enrichment analysis for genes found differentially',
                            'regulated in the hormonal clusters was performed',
                            'with the goana algorithm.',
                            'Enrichment p values were corrected for multiple',
                            'testing with the false discovery rate method',
                            '(FDR).',
                            'Odds ratio of enrichment in the differentially',
                            'regulated gene set as compared with the entire',
                            'genome served as enrichment effect size metric.',
                            'GO terms found to be significantly enriched in',
                            'the corresponding clusters of both the TCGA',
                            'and GSE99420 cohort are shown.',
                            'The table is available as a supplementary Excel',
                            'file.'))

# collecTRI regulons ------

  insert_msg('CollecTRI regulons')

  suppl_tabs$collectri <- bcg_collectri$test %>%
    map(~.x[c('tcga', 'gse99420')]) %>%
    map(compress, names_to = 'cohort') %>%
    map2(.,
         map(bcg_collectri$common_significant,
             reduce, union),
         ~filter(.x, source %in% .y)) %>%
    compress(names_to = 'clust_id') %>%
    mutate(cohort = globals$cohort_labs[cohort],
           score = signif(score, 2)) %>%
    select(clust_id,
           source,
           cohort,
           regulation,
           score,
           significance) %>%
    set_names(c('Cluster',
                'Regulon name',
                'Cohort',
                'Status vs cohort mean',
                'LM score',
                'Significance')) %>%
    mdtable(label = 'collectri',
            ref_name = 'collectri',
            caption = paste('Differential modulation of transcriptional',
                            'collecTRI regulons in the hormonal clusters as',
                            'compared with the cohort mean.',
                            'Analysis of regulon activity was investigated by',
                            'univariable linear modeling algorithm from the',
                            'decoupleR package fed with T statistic values of',
                            'differential gene expression of all available',
                            'genes.',
                            'Magnitude of differential regulon activity was',
                            'measured by linear modeling score (LM score).',
                            'P values of non-zero LM score were corrected for',
                            'multiple testing with the false discovery rate',
                            'method.',
                            'Regulons found to be differentially modulated',
                            'in both the TCGA and GSE99420 cohorts are shown.',
                            'The table is available as a supplementary',
                            'Excel file.'))

# PROGENy signaling pathways -------

  insert_msg('Progeny signaling pathways')

  suppl_tabs$progeny <- bcg_progeny$test %>%
    map(~.x[c('tcga', 'gse99420')]) %>%
    map(compress, names_to = 'cohort') %>%
    map2(.,
         map(bcg_progeny$common_significant,
             reduce, union),
         ~filter(.x, source %in% .y)) %>%
    compress(names_to = 'clust_id') %>%
    mutate(cohort = globals$cohort_labs[cohort],
           score = signif(score, 2)) %>%
    select(clust_id,
           source,
           cohort,
           regulation,
           score,
           significance) %>%
    arrange(clust_id, source, cohort) %>%
    set_names(c('Cluster',
                'Signaling pathway',
                'Cohort',
                'Status vs cohort mean',
                'LM score',
                'Significance')) %>%
    mdtable(label = 'progeny',
            ref_name = 'progeny',
            caption = paste('Differential regulation of PROGENy signaling',
                            'pathways in the hormonal clusters as',
                            'compared with the cohort mean.',
                            'Analysis of signaling pathway activity was',
                            'investigated by',
                            'multivariable linear modeling algorithm from the',
                            'decoupleR package fed with T statistic values of',
                            'differential gene expression of all available',
                            'genes.',
                            'Magnitude of differential pathway activity was',
                            'measured by linear modeling score (LM score).',
                            'P values of non-zero LM score were corrected for',
                            'multiple testing with the false discovery rate',
                            'method.',
                            'Signaling pathways found to be differentially',
                            'regulated in both the TCGA and GSE99420 cohorts',
                            'are shown.'))

# Differentially regulated metabolic reactions -------

  insert_msg('Metabolic reactions')

  suppl_tabs$reactions <- bcg_meta$estimates %>%
    map(~.x[c('tcga', 'gse99420')]) %>%
    map(compress, names_to = 'cohort') %>%
    map2(.,
         map(bcg_meta$common_significant,
             reduce, union),
         ~filter(.x, react_id %in% .y)) %>%
    compress(names_to = 'clust_id') %>%
    mutate(cohort = globals$cohort_labs[cohort],
           react_name = annotate_bigg(react_id),
           react_name = ifelse(is.na(react_name),
                               stri_replace(react_id,
                                            regex = '^R_',
                                            replacement = ''),
                               react_name),
           fold_reg = signif(log2(fold_reg), 2),
           lower_ci = signif(log2(lower_ci), 2),
           upper_ci = signif(log2(upper_ci), 2)) %>%
    select(clust_id,
           subsystem,
           react_id,
           react_name,
           regulation,
           cohort,
           fold_reg,
           upper_ci,
           lower_ci,
           p_adjusted) %>%
    arrange(clust_id,
            subsystem,
            react_id,
            cohort) %>%
    set_names(c('Cluster',
                'RECON subsystem',
                'Reaction ID',
                'Reaction name',
                'Status vs cohort mean',
                'Cohort',
                'log2 fold-regulation vs cohort mean',
                '95% CI, lower bound',
                '95% CI, upper bound',
                'pFDR')) %>%
    mdtable(label = 'metabolic_reactions',
            ref_name = 'reactions',
            caption = paste('Differential regulation of RECON2 model reactions',
                            'in the hormonal clusters as compared with the',
                            'cohort mean was investigated by Monte Carlo',
                            'simulation fed with log2 fold-regulation',
                            'estimates of differential gene expression and',
                            'their standard errors for all available genes.',
                            'Log2 fold-regulation estimates of reaction activity',
                            'with 95% confidence intervals (95% CI) and',
                            'false discovery rate (FDR) corrected p values are',
                            'presented for metabolic reactions found to be',
                            'significantly regulated both in the TCGA and',
                            'GSE99420 cohort.',
                            'The table is available as a supplementary',
                            'Excel file.'))

# Metabolic subsystem enrichment analysis -------

  insert_msg('Metabolic subsystem enrichment analysis')

  suppl_tabs$subsystems <- bcg_meta$enrichment %>%
    map(~.x[c('tcga', 'gse99420')]) %>%
    map(compress, names_to = 'cohort') %>%
    map2(.,
         map(bcg_meta$enrichment_common_significant,
             reduce, union),
         ~filter(.x,
                 subsystem %in% .y,
                 p_value < 0.05,
                 status %in% c('activated', 'inhibited'))) %>%
    compress(names_to = 'clust_id') %>%
    mutate(cohort = globals$cohort_labs[cohort]) %>%
    select(clust_id,
           subsystem,
           status,
           cohort,
           OR,
           p_value) %>%
    arrange(clust_id, subsystem, status, cohort) %>%
    set_names(c('Cluster',
                'RECON subsystem',
                'Reaction regulation mode',
                'Cohort',
                'Enrichment OR',
                'p value')) %>%
    map_dfc(function(x) if(is.numeric(x)) signif(x, 2) else x) %>%
    mdtable(label = 'subsystems',
            ref_name = 'subsystems',
            caption = paste('Enrichment of RECON metabolic subsystems with',
                            'significantly activated and inhibited reactions',
                            'in the hormonal subsets of testicular cancer.',
                            'Statistical significance was investigated by',
                            'comparing frequency of the subsystem reaction',
                            'within the activated or inhibited reaction set',
                            'with 100000 random draws from the total reaction',
                            'pool.',
                            'Odds ratio (OR) of enrichment in the regulated',
                            'reaction set over the entire reaction pool',
                            'served as an effect size metric.',
                            'Because p values decrease with increasing numbers',
                            'or random draws, no multiple testing correction',
                            'was applied.',
                            'Metabolic subsystems found significantly enriched',
                            'with activated or inhibited reactions in both the',
                            'TCGA and GSE99420 cohorts are presented.',
                            'The table is available as a supplementary Excel',
                            'file.'))

# Differential protein expression ------

  insert_msg('Differential protein expression')

  suppl_tabs$protein <- bcg_protein$dev_test %>%
    filter(regulation %in% c('upregulated', 'downregulated')) %>%
    select(clust_id, regulation,
           variable,
           deviation_center, deviation_sem,
           lower_ci, upper_ci,
           p_adjusted) %>%
    set_names(c('Cluster', 'Regulation vs cohort mean',
                'Protein name',
                'log2 fold-regulation vs cohort mean',
                'Standard error',
                '95% CI, lower bound', '95% CI, upper bound',
                'pFDR')) %>%
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x)

  suppl_tabs$protein <- suppl_tabs$protein %>%
    mdtable(label = 'protein',
            ref_name = 'protein',
            caption = paste('Expression of proteins in the hormonal clusters of',
                            'the TCGA cohorts was compared with the cohort means.',
                            'Statistical significance of differences in',
                            'log2-transformed expression levels between the',
                            'clusters was assessed by one-way ANOVA with',
                            'eta-square effect size statistic.',
                            'Statistical significance of differences between',
                            'log2-transformed expression in the cluster and',
                            'the cohort mean was determined by one-sample',
                            'two-tailed T test.',
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'Differential regulation in a particular cluster',
                            'was considered for',
                            'pFDR(ANOVA) < 0.05, eta-squared of at least 0.14,',
                            "and pFDR(T test) < 0.05.",
                            'The table is available as a supplementary',
                            'Excel file.'))

# Genetic alteration numbers --------

  insert_msg('Genetic alteration number')

  suppl_tabs$genetic_numbers <- bcg_burdens$result_tbl %>%
    mdtable(label = 'genetic_numbers',
            ref_name = 'genetic_numbers',
            caption = paste('Total mutation numbers, scores of microsatellite',
                            'instability, and numbers of gene mutations',
                            'deletions, and amplifications in the hormonal',
                            'clusters of the TCGA cohort.',
                            'Medians with interquartile ranges and ranges are',
                            'shown.'))

# Top genetic features of the clusters -------

  insert_msg('Top genetic features of the clusters')

  suppl_tabs$genetics <- bcg_genet$result_tbl %>%
    mdtable(label = 'genetics',
            ref_name = 'genetics',
            caption = paste('Frequencies of gene mutations, deletions, and',
                            'amplifications in the hormonal clusters of the',
                            'TCGA cohort.',
                            'Statistical significance of differences between',
                            'the clusters was determined by chi-square test with',
                            "Cramer's V effect size statistic.",
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'The table is available as a supplementary Excel',
                            'file.'))

# Predicted drug response ------

  insert_msg('Predicted drug response')

  suppl_tabs$drugs <- bcg_drugs$dev_test %>%
    map(map,
        filter,
        regulation %in% c('sensitive', 'resistant')) %>%
    map(map,
        select,
        clust_id, regulation,
        variable, drug_name,
        deviation_center, deviation_sem,
        lower_ci, upper_ci,
        p_adjusted) %>%
    map(compress,
        names_to = 'cohort') %>%
    compress(names_to = 'training_data') %>%
    mutate(cohort = globals$cohort_labs[cohort],
           training_data = globals$drug_exp_labs[training_data]) %>%
    relocate(training_data, cohort) %>%
    set_names(c('Training data set',
                'Cohort',
                'Cluster',
                'Drug response vs cohort mean',
                'Drug variable',
                'Drug name',
                'log2 fold-regulation vs cohort mean',
                'Standard error',
                '95% CI, lower bound',
                '95% CI, upper bound',
                'pFDR')) %>%
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x)

  suppl_tabs$drugs <- suppl_tabs$drugs %>%
    mdtable(label = 'drugs',
            ref_name = 'drugs',
            caption = paste('Drug response in form of log IC50 and area under',
                            'response curve (AUC) was predicted for cancer',
                            'samples by whole-transcriptome RIDGE linear models',
                            'trained with the CTRP2 and GDSC drug screening',
                            'data sets.',
                            'Statistical significance of differences in',
                            'the drug response metrics between the',
                            'clusters was assessed by one-way ANOVA with',
                            'eta-square effect size statistic.',
                            'Statistical significance of differences between',
                            'the drug response metrics in the cluster and',
                            'the cohort mean was determined by one-sample',
                            "two-tailed T test",
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'Differential drug response in a particular cluster',
                            'was considered for',
                            'pFDR(ANOVA) < 0.05, eta-squared of at least 0.14,',
                            "and pFDR(T test) < 0.05.",
                            'The table is available as a supplementary',
                            'Excel file.'))

# Saving tables on the disc ------

  insert_msg('Saving tables on the disc')

  suppl_tabs %>%
    save_excel(path = './paper/supplementary_tables.xlsx',
               prefix = 'Table S')

# END -------

  insert_tail()
