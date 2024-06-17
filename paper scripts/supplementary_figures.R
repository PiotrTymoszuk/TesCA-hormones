# Supplementary figures

  insert_head()

# container --------

  suppl_figs <- list()

# differential gene expression and histology details --------

  insert_msg('Histology details and gene expression')

  ## upper panel: heat map of Z-scores and rug heat map of gene classification

  suppl_figs$histology$upper <-
    plot_grid(expl_icd$hm_plots$tcga +
                theme(legend.position = 'none',
                      strip.text.y = element_blank(),
                      strip.background.y = element_blank(),
                      strip.text.x = element_text(angle = 90),
                      axis.text.y = element_text(size = 7)),
              expl_icd$rug_hm_plot +
                theme(legend.position = 'none'),
              ncol = 2,
              rel_widths = c(0.97, 0.03),
              align = 'h',
              axis = 'tblr') %>%
    plot_grid(plot_grid(get_legend(expl_icd$hm_plots$tcga +
                                     theme(legend.position = 'bottom',
                                           legend.text = element_text(size = 7),
                                           legend.title = element_text(size = 7))),
                        get_legend(expl_icd$rug_hm_plot +
                                     theme(legend.position = 'bottom',
                                           legend.text = element_text(size = 7),
                                           legend.title = element_text(size = 7))),
                        ncol = 2),
              nrow = 2,
              rel_heights = c(0.9, 0.1))

  ## bottom panel: selected genes specific for particular histological subtypes

  suppl_figs$histology$bottom <-
    expl_icd$plots$tcga[c("LHB", "CGA", "STAR",
                          "PRL", "SRD5A1", "SHBG")] %>%
    map(~.x +
          guides(x = guide_axis(angle = 45)) +
          theme(legend.position = 'none',
                axis.text.x = element_text(size = 7))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr')

  ## the entire figure

  suppl_figs$histology <-
    plot_grid(suppl_figs$histology$upper,
              suppl_figs$histology$bottom,
              nrow = 2,
              rel_heights = c(1.5, 1.5),
              labels = LETTERS,
              label_size = 10) %>%
    as_figure(label = 'differential_gene_expression_histology_icd',
              ref_name = 'histology',
              caption = paste('Differential expression of hormone-related',
                              'genes in histological subtypes of testicular',
                              'cancer in the TCGA cohort.'),
              w = 180,
              h = 230)

# Co-expression of the hormone-related genes --------

  insert_msg('Co-expression of the hormone-related genes')

  ## upper panel: graph plots

  suppl_figs$co_expression$upper <-
    expl_net$graph_plots[c("tcga", "gse99420")] %>%
    map(~.x +
          theme(legend.position = 'none',
                plot.subtitle = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(expl_net$graph_plots[[1]] +
                           theme(legend.position = 'bottom',
                                 #legend.box = 'vertical',
                                 legend.title = element_text(size = 7),
                                 legend.text = element_text(size = 7))),
              nrow = 2,
              rel_heights = c(0.85, 0.15))

  ## bottom panel: network stats

  suppl_figs$co_expression$bottom <-
    expl_net$stat_plots[c("tcga", "gse99420")] %>%
    map(~.x +
          theme(legend.position = 'none',
                plot.subtitle = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(expl_net$stat_plots[[1]] +
                           theme(legend.position = 'bottom',
                                 #legend.box = 'vertical',
                                 legend.title = element_text(size = 7),
                                 legend.text = element_text(size = 7))),
              nrow = 2,
              rel_heights = c(0.85, 0.15))

  ## the entire figure

  suppl_figs$co_expression <-
    plot_grid(suppl_figs$co_expression$upper,
              suppl_figs$co_expression$bottom,
              nrow = 2,
              rel_heights = c(1.1, 0.9),
              labels = LETTERS,
              label_size = 10) %>%
    as_figure(label = 'hormone_gene_co_expression',
              ref_name = 'co_expression',
              caption = paste('Networks of co-expressed hormone-related genes',
                              'in testicular carcinoma.'),
              w = 180,
              h = 230)

# Cluster development --------

  insert_msg('Cluster development')

  ## we're displaying only the TCGA and GSE99420 cohorts
  ## the GSE3218 results are available in the repository's R pipeline

  ## upper panel:
  ## statistics for different k number scenarios
  ## cluster sizes

  suppl_figs$clusters$upper <-
    list(clust_dev$stat_plots$train +
           guides(color = guide_legend(nrow = 2)) +
           labs(title = paste('Cluster number choice,',
                              globals$cohort_labs["tcga"]),
                subtitle = paste('samples: n =',
                                 nobs(clust_dev$algos$k4)$observations)),
         clust_eval$clust_size_plot)

  suppl_figs$clusters$upper[[2]]$data <-
    suppl_figs$clusters$upper[[2]]$data %>%
    filter(cohort %in% c('tcga', 'gse99420'))

  suppl_figs$clusters$upper <-
    suppl_figs$clusters$upper %>%
    map(~.x +
          theme(legend.position = 'bottom',
                legend.title = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr',
              labels = LETTERS,
              label_size = 10)

  ## bottom panel: pairwise distances between the observations
  ## and mean distances between the training and test cohort clusters

  suppl_figs$clusters$bottom <-
    c(clust_eval$homo_heat_maps[c("tcga", "gse99420")],
      clust_eval$hetero_heat_maps["gse99420"])

  suppl_figs$clusters$bottom <- suppl_figs$clusters$bottom %>%
    map(~.x +
          scale_fill_gradient2(low = 'firebrick',
                               mid = 'white',
                               high = 'steelblue',
                               midpoint = suppl_figs$clusters$bottom %>%
                                 map(~.x$data$mean) %>%
                                 reduce(c) %>%
                                 range %>%
                                 mean,
                               limits = suppl_figs$clusters$bottom %>%
                                 map(~.x$data$mean) %>%
                                 reduce(c) %>%
                                 range,
                               name = 'Mean\nsum-of-square\ndistance')) %>%
    map2(.,
         c(paste('Pairwise distances,',
                 globals$cohort_labs[c("tcga", "gse99420")]),
           paste('Cross-distances,',
                 globals$cohort_labs["tcga"],
                 'vs',
                 globals$cohort_labs["gse99420"])),
         ~.x +
           labs(title = .y) +
           theme(plot.subtitle = element_blank()))

  suppl_figs$clusters$bottom <- suppl_figs$clusters$bottom %>%
    map(~.x + theme(legend.position = 'none')) %>%
    c(list(get_legend(suppl_figs$clusters$bottom[[1]] +
                        theme(legend.position = 'bottom')))) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr',
              labels = 'C',
              label_size = 10)

  ## the entire figure

  suppl_figs$clusters <-
    plot_grid(suppl_figs$clusters$upper,
              suppl_figs$clusters$bottom,
              nrow = 2,
              rel_heights = c(1, 1.4)) %>%
    as_figure(label = 'cluster_development_evaluation',
              ref_name = 'clusters',
              caption = paste('Development of the hormonal clusters with the',
                              'TCGA training data set and evaluation of',
                              'the clustering structures in the TCGA',
                              'and GSE99420 cohorts.'),
              w = 180,
              h = 230)

# Key hormone-related genes for cluster discrimination ------

  insert_msg('Key hormone-related genes for the clusters')

  ## cluster #1: PRL, GNRH1
  ## cluster #2: SRD5A3
  ## cluster #3: CYP19A1, CGA
  ## cluster #4: STAR, HSD17B3

  suppl_figs$cluster_key <- clust_ft$plots[c("tcga", "gse99420")] %>%
    map(~.x[c('PRL', 'GNRH1',
              'SRD5A3',
              'CYP19A1', 'CGA',
              'STAR', 'HSD17B3')]) %>%
    transpose %>%
    unlist(recursive = FALSE) %>%
    map(~.x +
          theme(legend.position = 'none',
                axis.title.x = element_blank()))

  ## N numbers places in the common plot legends

  suppl_figs$cluster_key <-
    c(suppl_figs$cluster_key[1:6],
      list(get_legend(clust_ft$plots$tcga[[1]] +
                        labs(fill = globals$cohort_labs["tcga"])),
           get_legend(clust_ft$plots$gse99420[[1]] +
                        labs(fill = globals$cohort_labs["gse99420"]))),
      suppl_figs$cluster_key[7:length(suppl_figs$cluster_key)])

  suppl_figs$cluster_key <- suppl_figs$cluster_key %>%
    plot_grid(plotlist = .,
              ncol = 4,
              align = 'hv',
              axis = 'tblr') %>%
    as_figure(label = 'clusters_key_hormone_genes',
              ref_name = 'cluster_key',
              caption = paste('Expression of the key cluster-defining',
                              'hormone-related genes.'),
              w = 180,
              h = 210)

# Relapse-free survival and RIDGE modeling --------

  insert_msg('Relapse-free survival and RIDGE Cox modeling')

  ## top panel: RFS Kaplan-Meier analysis

  suppl_figs$surv$upper <- bcg_surv[c("plots", "post_hoc_plots")] %>%
    map(~.x$rfs +
          guides(color = guide_legend(nrow = 2)) +
          theme(legend.position = 'bottom',
                legend.title = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr')

  ## bottom panel: RIDGE Cox modeling of PFS

  suppl_figs$surv$bottom <-
    plot_grid(bcg_glmcox$stat_plot +
                theme(legend.position = 'bottom',
                      legend.box = 'vertical'),
              bcg_glmcox$coef_plots$full +
                scale_x_continuous(limits = c(-0.05, 0.15)) +
                theme(legend.position = 'bottom'),
              ncol = 2,
              align = 'h',
              axis = 'tblr')

  ## the entire figure

  suppl_figs$surv <- plot_grid(suppl_figs$surv$upper,
                               suppl_figs$surv$bottom,
                               nrow = 2,
                               labels = LETTERS,
                               label_size = 10) %>%
    as_figure(label = 'clusters_survival_rfs_pfs_modeling',
              ref_name = 'surv',
              caption = paste('Relapse-free survival in the hormonal clusters.',
                              'Multi-parameter modeling of progression-free',
                              'survival.'),
              w = 180,
              h = 210)

# Infiltration, MCP counter -------

  insert_msg('Infiltration, QuanTIseq and MCP counter')

  ## upper panel: quanTIseq
  ## bottom panel: MCP counter

  suppl_figs$infiltration$upper <-
    bcg_infil$plots$quantiseq[c("T cell CD8+", "B cell")]

  suppl_figs$infiltration$bottom <-
    bcg_infil$plots$mcp[c("Cancer associated fibroblast",
                          "T cell",
                          "T cell CD8+",
                          "B cell")]

  suppl_figs$infiltration <- suppl_figs$infiltration %>%
    map(map, ~.x[c('tcga', 'gse99420')]) %>%
    map(unlist, recursive = FALSE) %>%
    unlist(recursive = FALSE) %>%
    map(~.x +
          theme(legend.position = 'none',
                axis.title.x = element_blank(),
                plot.subtitle = element_text(size = 7)) +
          scale_x_discrete(labels = function(x) stri_replace(x,
                                                             regex = '\\n.*$',
                                                             replacement = '')))

  ## stitching the plot panel and appending it with legends
  ## containing the N number information

  suppl_figs$infiltration <-
    c(suppl_figs$infiltration,
      list(get_legend(bcg_infil$plots[[1]][[1]]$tcga +
                        labs(fill = globals$cohort_labs["tcga"])),
           get_legend(bcg_infil$plots[[1]][[1]]$gse99420 +
                        labs(fill = globals$cohort_labs["gse99420"])))) %>%
    plot_grid(plotlist = .,
              ncol = 4,
              align = 'hv',
              axis = 'tblr',
              labels = c('A', '', '', '',
                         'B', '', '', ''),
              label_size = 10) %>%
    as_figure(label = 'clusters_infiltration_quantiseq_mcp',
              ref_name = 'infiltration',
              caption = paste('Infiltration of cancer-associated fibroblasts,',
                              'T and B cells in the hormonal clusters predicted',
                              'by the QuanTIseq and MCP Counter algorithms.'),
              w = 180,
              h = 210)

# Differential gene expression: top genes per cluster -------

  insert_msg('Differential gene expression')

  suppl_figs$dge <- bcg_dgeplots$cmm_top_plots %>%
    map(~.x + theme(legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              ncol = 2) %>%
    as_figure(label = 'clusters_differential_gene_expression_common',
              ref_name = 'dge',
              caption = paste('Differential gene expression in the hormonal',
                              'clusters: common top regulated genes.'),
              w = 180,
              h = 230)

# GO enrichment --------

  insert_msg('GO enrichment')

  ## common significant GOs, semantic clustering results

  suppl_figs$go_clust <- bcg_gocplots$mds_plots %>%
    map(~.x +
          guides(fill = guide_legend(nrow = 2)) +
          theme(legend.position = 'bottom')) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    as_figure(label = 'clusters_go_enrichment',
              ref_name = 'go_clust',
              caption = paste('Biological process GO term enrichment in the',
                            'hormonal clusters.'),
              w = 180,
              h = 230)

# Differential expression of genes related to immunosuppression ------

  insert_msg('Genes related to immune checkpoint')

  ## common significant genes of relevance for immune
  ## therapy found significantly differentially regulated
  ## in both the TCGA and GSE99420 cohort

  suppl_figs$immun <- bcg_tex$box_plots[c("tcga", "gse99420")] %>%
    map(~.x[c('CTLA4', 'TIGIT', 'PDCD1')]) %>%
    unlist(recursive = FALSE) %>%
    map(~.x + theme(axis.title.x = element_blank(),
                    legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr') %>%
    as_figure(label = 'clusters_immune_checkpoint',
              ref_name = 'immun',
              caption = paste('Expression of genes related to immune checkpoint',
                              'in the hormonal clusters.'),
              w = 180,
              h = 140)

# Expression of cancer testis antigens and tumor markers ------

  insert_msg('Expression of cancer testis antigens')

  ## upper panel: heat maps

  suppl_figs$cta$upper <-
    plot_grid(bcg_cta$hm_plots$tcga +
                theme(strip.text.y = element_blank(),
                      strip.background.y = element_blank(),
                      legend.position = 'none',
                      axis.text.y = element_text(size = 7)),
              bcg_cta$hm_plots$gse99420 +
                theme(strip.text.y = element_blank(),
                      strip.background.y = element_blank(),
                      axis.text.y = element_blank(),
                      legend.position = 'none'),
              ncol = 2,
              align = 'h',
              axis = 'tblr',
              rel_widths = c(1.7, 1)) %>%
    plot_grid(get_legend(bcg_cta$hm_plots[[1]] +
                           theme(legend.position = 'bottom')),
              nrow = 2,
              rel_heights = c(0.87, 0.13))

  ## bottom panel: selected

  suppl_figs$cta$bottom <- bcg_cta$box_plots[c("tcga", "gse99420")] %>%
    map(~.x[c('AFP', 'CGA', 'LDHA', 'LDHC')]) %>%
    transpose %>%
    unlist(recursive = FALSE) %>%
    map(~.x +
          scale_x_discrete(labels = function(x) stri_replace(x,
                                                             regex = '\\n.*$',
                                                             replacement = '')) +
          theme(legend.position = 'none',
                axis.title.x = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 4,
              align = 'hv',
              axis = 'tblr')

  ## the entire figure

  suppl_figs$cta <-
    plot_grid(suppl_figs$cta$upper,
              suppl_figs$cta$bottom,
              nrow = 2,
              rel_heights = c(1.15, 0.85),
              labels = LETTERS,
              label_size = 10) %>%
    as_figure(label = 'clusters_cancer_antigens_markers',
              ref_name = 'cta',
              caption = paste('Expression of cancer testis antigens and',
                              'markers of testicular carcinoma in the hormonal',
                              'clusters.'),
              w = 180,
              h = 230)

# Expression of EGFR/ERBB/FGFR genes and the corresponding ligand genes -------

  insert_msg('Expression of ERBB/FGFR receptor and ligand genes')

  ## upper panel: heat maps for the ERBB2/FGFR system,
  ## common regulated genes

  suppl_figs$gfr$upper <-
    bcg_gfr$hm_plots[c("tcga", "gse99420")] %>%
    map(~.x +
          theme(legend.position = 'none',
                strip.background.y = element_blank(),
                strip.text.y = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(bcg_gfr$hm_plots[[1]] +
                           theme(legend.position = 'bottom')),
              nrow = 2,
              rel_heights = c(0.9, 0.1))

  ## bottom panel: box plots for the major sex hormone receptors
  ## found to be differentially regulated in both the TCGA and GSE99420 cohort

  suppl_figs$gfr$bottom <- bcg_hr$box_plots[c("tcga", "gse99420")] %>%
    map(~.x[c('FSHR', 'LHCGR', 'ESR2', 'AR')]) %>%
    transpose %>%
    unlist(recursive = FALSE) %>%
    map(~.x +
          scale_x_discrete(labels = function(x) stri_replace(x,
                                                             regex = '\\n.*',
                                                             replacement = '')) +
          theme(legend.position = 'none',
                    axis.title.x = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 4,
              align = 'hv',
              axis = 'tblr')

  ## the entire figure

  suppl_figs$gfr <- plot_grid(suppl_figs$gfr$upper,
                              suppl_figs$gfr$bottom,
                              nrow = 2,
                              rel_heights = c(1, 1),
                              labels = LETTERS,
                              label_size = 10) %>%
    as_figure(label = 'clusters_gf_hormone_receptors',
              ref_name = 'gfr',
              caption = paste('Expression of genes coding of ERBB and FGFR',
                              'family receptors, ERBB/FGFR ligands, and',
                              'receptors for gonadotropins and sex hormones in',
                              'the hormonal clusters.'),
              w = 180,
              h = 230)

# Differential expression of estrogen- and androgen-responsive genes -------

  insert_msg('Differential expression of androgen- and estrogen-responsive genes')

  suppl_figs$erar_resp <-
    c(bcg_estro$hm_plots[c("tcga", "gse99420")],
      bcg_andro$hm_plots[c("tcga", "gse99420")]) %>%
    map(~.x +
          theme(legend.position = 'none',
                axis.text.y = element_markdown(size = 7),
                strip.background.y = element_blank(),
                strip.text.y = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr',
              labels = c('A', '', 'B', ''),
              label_size = 10) %>%
    plot_grid(get_legend(bcg_andro$hm_plots[[1]] +
                           theme(legend.position = 'bottom')),
              nrow = 2,
              rel_heights = c(0.9, 0.1)) %>%
    as_figure(label = 'estrogen_androgen_responsive_genes',
              ref_name = 'erar_resp',
              caption = paste('Expression of estrogen- and',
                              'androgen-responsive genes in the hormonal',
                              'clusters.'),
              w = 180,
              h = 230)

# Matrisome genes -------

  insert_msg('Differential expression of matrisome genes')

  suppl_figs$ecm <-
    plot_grid(bcg_ecm$hm_plots$tcga +
                theme(legend.position = 'none',
                      strip.text.y = element_blank(),
                      strip.background.y = element_blank(),
                      axis.text.y = element_markdown(size = 7)),
              bcg_ecm$hm_plots$gse99420 +
                theme(legend.position = 'none',
                      strip.text.y = element_blank(),
                      strip.background.y = element_blank(),
                      axis.text.y = element_blank()),
              bcg_ecm$hm_rug +
                theme(legend.position = 'none'),
              ncol = 3,
              align = 'h',
              axis = 'tblr',
              rel_widths = c(1.77, 1, 0.12)) %>%
    plot_grid(plot_grid(get_legend(bcg_ecm$hm_plots[[1]] +
                                     theme(legend.position = 'bottom')),
                        get_legend(bcg_ecm$hm_rug +
                                     guides(fill = guide_legend(ncol = 2)) +
                                     theme(legend.position = 'bottom')),
                        ncol = 2),
              nrow = 2,
              rel_heights = c(0.9, 0.1)) %>%
    as_figure(label = 'clusters_matrisome',
              ref_name = 'ecm',
              caption = paste('Expression of matrisome genes in the hormonal',
                              'clusters of testicular carcinoma.'),
              w = 180,
              h = 230)

# Metabolic subsystem enrichment ------

  insert_msg('Metabolism subsystem enrichment')

  suppl_figs$metabolism <- bcg_subplots$bar_plots %>%
    map(~.x +
          theme(plot.subtitle = element_blank(),
                legend.position = 'none',
                strip.background.y = element_blank(),
                strip.text.y = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(bcg_subplots$bar_plots[[1]] +
                           theme(legend.position = 'bottom')),
              nrow = 2,
              rel_heights = c(0.92, 0.08)) %>%
    as_figure(label = 'clusters_metabolism',
              ref_name = 'metabolism',
              caption = paste('Metabolic subsystem enrichment analysis',
                              'for the hormonal clusters.'),
              w = 180,
              h = 230)

# Energy metabolism: TCA and Oxphos ------

  insert_msg('Energy metabolism reactions')

  suppl_figs$energy_reactions <-
    plot_grid(bcg_metaplots$common_bar_plots$`Citric acid cycle` +
                theme(legend.position = 'none'),
              bcg_metaplots$common_bar_plots$`Oxidative phosphorylation` +
                theme(legend.position = 'none'),
              nrow = 2,
              align = 'hv',
              axis = 'tblr',
              rel_heights = c(1, 0.7),
              labels = LETTERS,
              label_size = 10) %>%
    plot_grid(get_legend(bcg_metaplots$common_bar_plots[[1]]),
              ncol = 2,
              rel_widths = c(0.85, 0.15)) %>%
    as_figure(label = 'clusters_energy_metabolism_reactions',
              ref_name = 'energy_reations',
              caption = paste('Predicted regulation of reactions of citric acid',
                              'cycle and oxidative phosphorylation in the',
                              'hormonal clusters.'),
              w = 180,
              h = 230)

# Steroid metabolism reactions -------

  insert_msg('Steroid metabolism reactions')

  suppl_figs$steroid_reactions <-
    plot_grid(bcg_metaplots$common_bar_plots$`Steroid metabolism` +
                theme(legend.position = 'none',
                      strip.text.y = element_blank(),
                      strip.background.y = element_blank(),
                      axis.text.y = element_text(size = 7)),
              bcg_metaplots$rug_plot +
                theme(legend.position = 'none',
                      plot.margin = ggplot2::margin(l = 0,
                                                    r = 2,
                                                    unit = 'mm')),
              ncol = 2,
              align = 'h',
              axis = 'tblr',
              rel_widths = c(0.9, 0.1)) %>%
    plot_grid(plot_grid(get_legend(bcg_metaplots$common_bar_plots$`Steroid metabolism` +
                                     theme(legend.position = 'right')),
                        get_legend(bcg_metaplots$rug_plot +
                                     guides(fill = guide_legend(ncol = 3)) +
                                     theme(legend.position = 'right')),
                        ncol = 2),
              nrow = 2,
              rel_heights = c(0.89, 0.11)) %>%
    as_figure(label = 'clusters_steroid_metabolism_reactions',
              ref_name = 'steroid_reactions',
              caption = paste('Predicted activity of steroid metabolism',
                              'reactions in the hormonal clusters.'),
              w = 180,
              h = 230)

# Expression of proteins of the ERBB - AKT - mTOR pathway -------

  insert_msg('Expression of proteins of the ERBB - AKT - mTOR pathway')

  suppl_figs$gf_signaling <-
    plot_grid(bcg_protein$detail_plot +
                scale_x_continuous(limits = c(-1, 1)) +
                theme(legend.position = 'bottom')) %>%
    as_figure(label = 'clusters_growth_factor_signaling_proteins',
              ref_name = 'gf_signaling',
              caption = paste('Regulation of growth factor pathway signaling',
                              'proteins in the hormonal clusters.'),
              w = 180,
              h = 210)

# Protein networks --------

  insert_msg('Protein networks')

  suppl_figs[c('pro_net_clust1_2',
               'pro_net_clust3_4')] <-
    list(bcg_pronet$graph_plots$comm_id[c("#1", "#2")],
         bcg_pronet$graph_plots$comm_id[c("#3", "#4")]) %>%
    map(map,
        ~.x +
          theme(legend.position = 'none',
                plot.subtitle = element_blank())) %>%
    map(~plot_grid(plotlist = .,
                   nrow = 2)) %>%
    map(~plot_grid(.x,
                   get_legend(bcg_pronet$graph_plots$comm_id[[1]] +
                                guides(fill = 'none',
                                       color = 'none') +
                                theme(legend.position = 'right')),
                   ncol = 2,
                   rel_widths = c(0.85, 0.15)))

  suppl_figs[c('pro_net_clust1_2',
               'pro_net_clust3_4')] <-
    suppl_figs[c('pro_net_clust1_2',
                 'pro_net_clust3_4')] %>%
    list(x = .,
         label = c('cluster1_2_protein_networks',
                   'cluster3_4_protein_networks'),
         ref_name = names(.),
         caption = paste('Co-expression protein networks in the hormonal',
                         'subsets of the TCGA cohort:',
                         c('clusters #1 and #2.',
                           'clusters #3 and #4.'))) %>%
    pmap(as_figure,
         w = 180,
         h = 200)

# Burdens of genetic alterations ------

  insert_msg('Burdens of genetic alterations')

  suppl_figs$burdens <-
    map2(bcg_burdens$plots,
         c(rep('identity', 4),
           rep('sqrt', 2)),
         ~.x +
           scale_y_continuous(trans = .y) +
           theme(legend.position = 'none',
                 axis.title.x = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr') %>%
    as_figure(label = 'clusters_genetic_alteration_numbers',
              ref_name = 'burdens',
              caption = paste('Numbers of genetic alterations in the hormonal',
                              'subsets.'),
              w = 180,
              h = 140)

# Genetics -------

  insert_msg('Genetics of the clusters')

  ## bar plot panels with alterations differing between the clusters
  ## with effect size of at least 0.2

  suppl_figs$genetics <-
    bcg_genet$panels[c("mutations", "deletions", "amplifications")] %>%
    map(~.x +
          theme(strip.text.y = element_blank(),
                strip.background.y = element_blank()))

  suppl_figs$genetics <-
    plot_grid(suppl_figs$genetics$mutations,
              suppl_figs$genetics$deletions,
              nrow = 3,
              rel_heights = c(5, 4, 2)) %>%
    plot_grid(suppl_figs$genetics$amplifications,
              ncol = 2,
              rel_widths = c(1, 1.13)) %>%
    as_figure(label = 'clusters_top_alterations',
              ref_name = 'genetics',
              caption = 'Cluster-specific genetic alterations.',
              w = 180,
              h = 180)

# Drugs: CTRP2 -------

  insert_msg('Drugs CTRP2')

  ## upper panel: heat maps

  suppl_figs$drugs$upper <-
    list(bcg_drdetails$hm_plots$ctrp2$tcga +
           guides(y = guide_axis(n.dodge = 2,
                                 check.overlap = TRUE)) +
           theme(axis.text.y = element_text(size = 7),
                 plot.subtitle = element_blank()),
         bcg_drdetails$hm_plots$ctrp2$gse99420 +
           theme(axis.text.y = element_blank(),
                 plot.subtitle = element_blank()),
         bcg_drdetails$rug_hm_plots$ctrp2 +
           theme(plot.margin = ggplot2::margin(l = 0,
                                               r = 2,
                                               unit = 'mm'))) %>%
    map(~.x +
          theme(legend.position = 'none',
                strip.background.y = element_blank(),
                strip.text.y = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'h',
              axis = 'tblr',
              rel_widths = c(1.6, 1, 0.1))

  ## middle panel: legends

  suppl_figs$drugs$middle <-
    plot_grid(get_legend(bcg_drdetails$hm_plots$ctrp2[[1]] +
                           theme(legend.position = 'bottom',
                                 legend.title = element_text(size = 7),
                                 legend.text = element_text(size = 7))),
              get_legend(bcg_drdetails$rug_hm_plots$ctrp2 +
                           guides(fill = guide_legend(ncol = 4)) +
                           theme(legend.position = 'bottom',
                                 legend.title = element_text(size = 7),
                                 legend.text = element_text(size = 7))),
              ncol = 2,
              rel_widths = c(0.6, 1.4))

  ## bottom panel: selected drugs

  suppl_figs$drugs$bottom <-
    bcg_drdetails$plots$ctrp2[c("tcga", "gse99420")] %>%
    map(~.x[c('ABITREXATE (CTRP:30371)',
              'DASATINIB (CTRP:52882)',
              'GX15-070 (CTRP:606142)')]) %>%
    unlist(recursive = FALSE) %>%
    map(~.x +
          labs(y = 'predicted AUC') +
          theme(legend.position = 'none',
                axis.title.x = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr')

  ## the entire figure

  suppl_figs$drugs <- plot_grid(suppl_figs$drugs$upper,
                                suppl_figs$drugs$middle,
                                suppl_figs$drugs$bottom,
                                nrow = 3,
                                rel_heights = c(1.1, 0.21, 0.9),
                                labels = c('A', '', 'B'),
                                label_size = 10) %>%
    as_figure(label = 'cluster_predicted_drug_response_ctrp2',
              ref_name = 'drugs',
              caption = paste('Predicted anti-cancer drug response in the',
                              'hormonal clusters.'),
              w = 180,
              h = 230)

# Result summary -------

  insert_msg('Result summary')

  suppl_figs$summary <-
    plot_grid(ggdraw() +
                draw_image('./schemes/result_summary.png') +
                theme(plot.margin = ggplot2::margin(2, 2, 2, 2,
                                                    unit = 'mm'))) %>%
    as_figure(label = 'result_summary',
              ref_name = 'summary',
              caption = 'Summary of the analysis results.',
              w = 180,
              h = 2684/2724 * 180)

# Saving the supplements on the disc ------

  insert_msg('Saving on the disc')

  suppl_figs %>%
    number_figures(prefix = 'supplementary_figure_') %>%
    walk(pickle,
         path = './paper/supplementary figures',
         device = cairo_pdf)

# END ------

  insert_tail()
