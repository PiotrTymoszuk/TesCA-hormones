# Figures for the main manuscript

  insert_head()

# container -------

  figs <- list()

# Hormone-related genes and histology -------

  insert_msg('Hormone-related genes and histology')

  ## upper panel: Volcano plots

  figs$histology$upper <-
    expl_histo$volcano_plots[c("tcga", "gse99420")] %>%
    map(~.x +
          theme(legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(expl_histo$volcano_plots[[1]] +
                           guides(alpha = 'none') +
                           theme(legend.position = 'bottom',
                                 legend.title = element_text(size = 7),
                                 legend.text = element_text(size = 7))),
              nrow = 2,
              rel_heights = c(0.85, 0.15))

  ## bottom panel: selected common regulated genes

  figs$histology$bottom <- expl_histo$plots[c("tcga", "gse99420")] %>%
    map(~.x[c('SRD5A1', 'HSD17B2', 'SHBG', 'SRD5A3')]) %>%
    transpose %>%
    unlist(recursive = FALSE) %>%
    map(~.x + theme(legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              ncol = 4,
              align = 'hv',
              axis = 'tblr')

  ## the entire figure

  figs$histology <-
    plot_grid(figs$histology$upper,
              figs$histology$bottom,
              nrow = 2,
              rel_heights = c(1.3, 1.7),
              labels = LETTERS,
              label_size = 10) %>%
    as_figure(label = 'differential_gene_expression_histology',
              ref_name = 'histology',
              caption = paste('Differential expression of hormone-related',
                              'genes in NSGCT and seminoma.'),
              w = 180,
              h = 210)

# Cluster-defining factors in the hormonal clusters -------

  insert_msg('Cluster-defining factors in the hormonal clusters')

  ## upper panel: heat maps

  figs$clusters$upper <-
    list(ft_hm = clust_ft$hm_plots[c("tcga", "gse99420")],
         sample_rug = clust_ft$hm_sample[c("tcga", "gse99420")]) %>%
    pmap(make_ft_heat_map,
         gene_rug = clust_ft$hm_gene) %>%
    plot_grid(plotlist = .,
              ncol = 2)

  ## middle panel: legends

  figs$clusters$middle <-
    plot_grid(get_legend(clust_ft$hm_plots[[1]] +
                           theme(legend.position = 'bottom',
                                 legend.text = element_text(size = 7.5))),
              get_legend(clust_ft$hm_sample[[1]] +
                           theme(legend.position = 'bottom',
                                 legend.text = element_text(size = 7.5))),
              get_legend(clust_ft$hm_gene +
                           theme(legend.position = 'bottom',
                                 legend.text = element_text(size = 7.5))),
              nrow = 2)

  ## bottom panel: histology

  figs$clusters$bottom <- bcg_clinic$plots[c("tcga", "gse99420")] %>%
    map(~.x[c('histology', 'histology_icd')]) %>%
    unlist(recursive = FALSE) %>%
    compact %>%
    map2(., c(1, 2, 1),
         ~.x +
           guides(fill = guide_legend(nrow = .y)) +
           theme(legend.position = 'bottom',
                 axis.title.x = element_blank(),
                 legend.box = 'vertical',
                 legend.text = element_text(size = 7.5))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr')

  ## the entire figure

  figs$clusters <-
    plot_grid(figs$clusters$upper,
              figs$clusters$middle +
                theme(plot.margin = ggplot2::margin(l = 10, r = 2, unit = 'mm')),
              figs$clusters$bottom,
              nrow = 3,
              rel_heights = c(0.85, 0.2, 0.6),
              labels = c('A', '', 'B'),
              label_size = 10) %>%
    as_figure(label = 'hormonal_cluster_definition',
              ref_name = 'clusters',
              caption = 'Hormonal clusters of testicular cancers.',
              w = 180,
              h = 230)

# Clinical features and survival in the hormonal clusters -------

  insert_msg('Clinics and survival in the clusters')

  ## upper panel: age, tumor stage, marker status

  figs$clinic$upper <-
    map2(bcg_clinic$plots$tcga[c("age", "pt_stage", "marker_status")],
         c('none', 'bottom', 'bottom'),
         ~.x +
           theme(legend.position = .y,
                 axis.title.x = element_blank(),
                 legend.title = element_text(size = 7.2),
                 legend.text = element_text(size = 7.2))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr')

  ## bottom panel: survival

  figs$clinic$bottom <- bcg_surv[c("plots", "post_hoc_plots")] %>%
    map(~.x$pfs +
          guides(color = guide_legend(nrow = 2)) +
          theme(legend.position = 'bottom',
                legend.title = element_blank(),
                legend.text = element_text(size = 7.2))) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr')

  ## the entire figure

  figs$clinic <- plot_grid(figs$clinic$upper,
                           figs$clinic$bottom,
                           nrow = 2,
                           rel_heights = c(1, 1.3),
                           labels = LETTERS,
                           label_size = 10) %>%
    as_figure(label = 'clusters_clinic_survival',
              ref_name = 'clinic',
              caption = paste('Clinical and prognostic characteristic of the',
                              'hormonal clusters of testicular cancer.'),
              w = 180,
              h = 180)

# Infiltration and Reactome pathways -------

  insert_msg('Infiltration and Reactome pathways')

  ## upper panel: infiltration, xCell

  figs$biology$upper <-
    bcg_infil$plots$xcell[c("Cancer associated fibroblast",
                            "T cell CD8+",
                            "Class-switched memory B cell")] %>%
    map(~.x[c('tcga', 'gse99420')]) %>%
    transpose %>%
    unlist(recursive = FALSE) %>%
    map(~.x +
          theme(legend.position = 'none',
                axis.title.x = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr')

  ## bottom panel: Reactome pathways

  figs$biology$bottom <-
    bcg_reactome$hm_plots[c("tcga", "gse99420")] %>%
    map(~.x + theme(legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(bcg_reactome$hm_plots[[1]] +
                           theme(legend.position = 'bottom')),
              nrow = 2,
              rel_heights = c(0.9, 0.1))

  ## the entire figure

  figs$biology <- plot_grid(figs$biology$upper,
                            figs$biology$bottom,
                            nrow = 2,
                            rel_heights = c(1, 1.3),
                            labels = LETTERS,
                            label_size = 10) %>%
    as_figure(label = 'clusters_tme_reactome_pathways',
              ref_name = 'biology',
              caption = paste('Tumor microenvironment composition and',
                              'Reactome pathway gene signatures in the hormonal',
                              'clusters of testicular carcinoma.'),
              w = 180,
              h = 230)

# Regulons and signaling pathways -------

  insert_msg('Regulons and signaling')

  ## upper part: collecTRI regulons

  figs$signaling$upper <- bcg_collectri$top_plot +
    theme(plot.subtitle = element_blank())

  ## bottom part: PROGENy signaling pathways

  figs$signaling$bottom <-
    bcg_progeny$bubble_plots[c("tcga", "gse99420")] %>%
    map(~.x +
          theme(legend.position = 'none',
                plot.subtitle = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(bcg_progeny$bubble_plots[[1]]),
              ncol = 2,
              rel_widths = c(0.85, 0.15))

  ## the entire figure

  figs$signaling <- plot_grid(figs$signaling$upper,
                              figs$signaling$bottom,
                              nrow = 2,
                              rel_heights = c(1.4, 1),
                              labels = LETTERS,
                              label_size = 10) %>%
    as_figure(label = 'clusters_regulons_signaling',
              ref_name = 'signaling',
              caption = paste('Modulation of transcriptional regulons',
                              'and signaling pathways in the hormonal clusters',
                              'of testicular cancer.'),
              w = 180,
              h = 220)

# Protein markers of the clusters -------

  insert_msg('Protein markers of the clusters')

  ## upper panel: volcano plots

  figs$protein$upper <- bcg_protein$volcano_plots %>%
    map(~.x +
          theme(legend.position = 'none',
                plot.subtitle = element_text(size = 7))) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr')

  ## bottom panel, box plots for top markers
  ## YAP1|YAP_pS127
  ## KIT|c-Kit
  ## SRC|SRCpY527
  ## FN1|Fibrinectin
  ## ERBB2|HER2
  ## VHL|VHL

  figs$protein$bottom <-
    bcg_protein$box_plots[c("YAP1|YAP_pS127",
                            "SRC|Src_pY527",
                            "ERBB2|HER2",
                            "KIT|c-Kit",
                            "FN1|Fibronectin",
                            "VHL|VHL")] %>%
    map(~.x +
          theme(legend.position = 'none',
                axis.title.x = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr')

  ## the entire figure

  figs$protein <- plot_grid(figs$protein$upper,
                            figs$protein$bottom,
                            nrow = 2,
                            rel_heights = c(1.1, 0.9),
                            labels = LETTERS,
                            label_size = 10) %>%
    as_figure(label = 'clusters_protein_expression',
              ref_name = 'protein',
              caption = paste('Differential protein expression in the hormonal',
                               'clusters of testicular cancers in the TCGA',
                               'cohort.'),
              w = 180,
              h = 230)

# Predicted drug response -------

  insert_msg('Predicted drug response')

  ## GDSC-trained predictions are shown in heat maps

  ## upper panel: heat maps

  figs$drugs$upper <-
    list(bcg_drdetails$hm_plots$gdsc$tcga +
           guides(y = guide_axis(n.dodge = 2,
                                 check.overlap = TRUE)) +
           theme(axis.text.y = element_text(size = 6),
                 plot.subtitle = element_blank()),
         bcg_drdetails$hm_plots$gdsc$gse99420 +
           theme(axis.text.y = element_blank(),
                 plot.subtitle = element_blank()),
         bcg_drdetails$rug_hm_plots$gdsc +
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
              rel_widths = c(1.8, 1, 0.1))

  ## middle panel: legends

  figs$drugs$middle <-
    plot_grid(get_legend(bcg_drdetails$hm_plots$gdsc[[1]] +
                           theme(legend.position = 'bottom',
                                 legend.title = element_text(size = 7),
                                 legend.text = element_text(size = 7))),
              get_legend(bcg_drdetails$rug_hm_plots$gdsc +
                           guides(fill = guide_legend(ncol = 4)) +
                           theme(legend.position = 'bottom',
                                 legend.title = element_blank(),
                                 legend.text = element_text(size = 7))),
              ncol = 2,
              rel_widths = c(0.6, 1.4))

  ## bottom panel: selected drugs

  figs$drugs$bottom <-
    bcg_drdetails$plots$gdsc[c("tcga", "gse99420")] %>%
    map(~.x[c('Vinblastine_1004', 'Ponatinib_155', 'Venetoclax_1909')]) %>%
    unlist(recursive = FALSE) %>%
    map(~.x +
          labs(y = 'predicted log IC50 [µM]') +
          theme(legend.position = 'none',
                axis.title.x = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr')

  ## the entire figure

  figs$drugs <- plot_grid(figs$drugs$upper,
                          figs$drugs$middle,
                          figs$drugs$bottom,
                          nrow = 3,
                          rel_heights = c(1.1, 0.21, 0.9),
                          labels = c('A', '', 'B'),
                          label_size = 10) %>%
    as_figure(label = 'cluster_predicted_drug_response',
              ref_name = 'drugs',
              caption = paste('Predicted anti-cancer drug response in the',
                              'hormonal clusters.'),
              w = 180,
              h = 230)

# Saving figures on the disc -------

  insert_msg('Saving figures on the disc')

  figs %>%
    number_figures %>%
    walk(pickle,
         path = './paper/figures',
         device = cairo_pdf)

# END ------

  insert_tail()
