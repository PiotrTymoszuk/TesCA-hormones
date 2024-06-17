# Differential expression of proteins between the hormonal clusters in the
# TCGA cohort.
#
# Statistical significance is determined by one-way ANOVA with eta-square
# effect size

  insert_head()

# container ------

  bcg_protein <- list()

# analysis data -------

  insert_msg('Analysis data')

  bcg_protein$lexicon <- tcga$protein_annotation

  bcg_protein$data <-
    inner_join(bcg_globals$assignment$tcga,
               tcga$protein,
               by = 'sample_id') %>%
    filter(tissue == 'tumor') %>%
    select(-tissue)

  bcg_protein$clust_levels <- levels(bcg_protein$data$clust_id)

# Descriptive stats ------

  insert_msg('Descriptive stats')

  bcg_protein$stats <- bcg_protein$data %>%
    column_to_rownames('sample_id') %>%
    fast_num_stats(split_factor = 'clust_id')

# ANOVA --------

  insert_msg('ANOVA')

  bcg_protein$anova <- bcg_protein$data %>%
    test_anova(split_fct = 'clust_id',
               variables = bcg_protein$lexicon$variable,
               adj_method = 'BH',
               .parallel = FALSE)

  ## formatting the results and identification of significant
  ## effects in ANOVA

  bcg_protein$anova <- bcg_protein$anova$anova %>%
    re_adjust(method = 'none') %>%
    mutate(variable = response,
           plot_cap = paste('\u03B7\u00B2 =', signif(effect_size, 2)),
           plot_cap = paste(plot_cap, significance, sep = ', '))

  bcg_protein$anova_significant  <- bcg_protein$anova %>%
    filter(p_adjusted < 0.05,
           effect_size >= 0.14) %>%
    .$variable

# Deviation from the cohort mean -------

  insert_msg('Deviation from the cocohort mean')

  bcg_protein$dev_test <-
    avg_deviation(bcg_protein$data,
                  split_fct = 'clust_id',
                  variables = bcg_protein$lexicon$variable)

  ## formatting, indicating the significant effects:
  ## significance in ANOVA and significant difference of expression
  ## in the cluster as compared with the cohort mean

  bcg_protein$dev_test <- bcg_protein$dev_test %>%
    mutate(anova_significant = ifelse(variable %in% bcg_protein$anova_significant,
                                      'yes', 'no'),
           regulation = ifelse(p_adjusted >= 0.05 | anova_significant == 'no',
                               'ns',
                               ifelse(deviation_center > 0,
                                      'upregulated',
                                      ifelse(deviation_center < 0,
                                             'downregulated', 'ns'))),
           regulation = factor(regulation,
                               c('upregulated', 'downregulated', 'ns')))

# Differentially regulated proteins -------

  insert_msg('Differentially regulated proteins')

  bcg_protein$significant <- bcg_protein$dev_test %>%
    filter(regulation %in% c('upregulated', 'downregulated')) %>%
    blast(clust_id, regulation) %>%
    map(~.x$variable)

# Volcano plots -------

  insert_msg('Volcano plots')

  bcg_protein$volcano_plots <-
    list(data = bcg_protein$dev_test %>%
           filter(anova_significant == 'yes') %>%
           blast(clust_id),
         plot_title = paste('RPPA proteins, cluster',
                            bcg_protein$clust_levels,
                            ', ', globals$cohort_labs["tcga"])) %>%
    pmap(plot_volcano,
         regulation_variable = 'deviation_center',
         p_variable = 'p_adjusted',
         regulation_level = 0,
         label_variable = 'variable',
         top_regulated = 7,
         label_type = 'text',
         txt_size = 2,
         txt_face = 'plain',
         cust_theme = globals$common_theme,
         x_lab = expression('log'[2] * ' fold-regulation vs cohort mean'),
         y_lab = expression('-log'[10] * ' pFDR'))  %>%
    map(~.x +
          labs(subtitle = .x$labels$tag) +
          theme(plot.tag = element_blank()))

# Forest plots for the top regulated genes in the clusters ------

  insert_msg('Top DGE Forest plots')

  bcg_protein$top_plots <-
    list(data = bcg_protein$dev_test %>%
           filter(regulation %in% c('upregulated', 'downregulated')) %>%
           blast(clust_id),
         plot_title = paste('Proteins, cluster',
                            bcg_protein$clust_levels,
                            ', ', globals$cohort_labs["tcga"])) %>%
    pmap(plot_top,
         regulation_variable = 'deviation_center',
         label_variable = 'variable',
         p_variable = 'p_adjusted',
         lower_ci_variable = 'lower_ci',
         upper_ci_variable = 'upper_ci',
         top_regulated = 10,
         fill_title = 'Regulation\nvs cohort mean',
         x_lab = expression('log'[2] * ' fold-regulation vs cohort mean'),
         cust_theme = globals$common_theme)

# Protein classification and heat maps of expression Z scores ------

  insert_msg('Classification and heat maps of expression Z scores')

  bcg_protein$classification <-
    classify(bcg_protein$data,
             variables = reduce(bcg_protein$significant, union),
             split_fct = 'clust_id') %>%
    .$classification

  ## differentially expressed proteins are presented

  bcg_protein$hm_plot <-
    heat_map(bcg_protein$data,
             variables = reduce(bcg_protein$significant, union),
             split_fct = 'clust_id',
             normalize = TRUE,
             variable_classification = bcg_protein$classification,
             plot_title = paste('Protein expression,',
                                globals$cohort_labs["tcga"]),
             x_lab = 'caner_sample',
             cust_theme = globals$common_theme +
               theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.y = element_blank()),
             limits = c(-3, 3),
             midpoint = 0,
             oob = scales::squish) +
    guides(y = guide_axis(n.dodge = 3))

# Box plots for the top protein markers of the clusters ------

  insert_msg('Box plots')

  ## proteins to be plotted and their stats
  ## appending with the ER and AR, EGFR and ERBBs,
  ## which are interesting as well

  bcg_protein$top_markers <- bcg_protein$classification %>%
    group_by(clust_id) %>%
    top_n(n = 5, delta_auc) %>%
    .$variable

  bcg_protein$top_markers <-
    c(bcg_protein$top_markers,
      'NRAS|N-Ras',
      'PIK3CA |PI3K-p110-alpha',
      'PIK3R1 PIK3R2|PI3K-p85',
      'AKT1 AKT2 AKT3|Akt_pS473',
      'AKT1 AKT2 AKT3|Akt_pT308',
      'MTOR|mTOR_pS2448',
      'RICTOR|Rictor',
      'RICTOR|Rictor_pT1135',
      'AR|AR',
      'ESR1|ER-alpha',
      'ESR1|ER-alpha_pS118',
      'EGFR|EGFR',
      'EGFR|EGFR_pY1068',
      'EGFR|EGFR_pY1173',
      'ERBB2|HER2',
      'ERBB2|HER2_pY1248',
      'ERBB3|HER3',
      'ERBB3|HER3_pY1298',
      'BRAF|B-Raf',
      'RAF1|C-Raf',
      'RAF1|C-Raf_pS338',
      'MAP2K1|MEK1',
      'MAP2K1|MEK1_pS217_S221',
      'MAPK1|ERK2',
      'MAPK1 MAPK3|MAPK_pT202_Y204',
      'BCL2L11|Bim',
      'BCL2L1|Bcl-xL',
      'BCL2|Bcl-2',
      'BID|Bid',
      'BAX|Bax',
      'BAD|Bad_pS112',
      'PTEN|PTEN',
      'TSC1|TSC1',
      'TSC2|Tuberin',
      'TSC2|Tuberin_pT1462')

  bcg_protein$top_captions <- bcg_protein$anova %>%
    filter(variable %in% bcg_protein$top_markers) %>%
    mutate(variable = factor(variable, bcg_protein$top_markers)) %>%
    arrange(variable) %>%
    .$plot_cap

  ## box plots

  bcg_protein$box_plots <-
    list(variable = bcg_protein$top_markers,
         plot_title = bcg_protein$top_markers %>%
           stri_replace_all(fixed = '_', replacement = '') %>%
           paste(globals$cohort_labs["tcga"], sep = ', '),
         plot_subtitle = bcg_protein$top_captions) %>%

    pmap(plot_variable,
         bcg_protein$data,
         split_factor = 'clust_id',
         type = 'box',
         cust_theme = globals$common_theme,
         x_lab = 'cluster',
         y_lab = expression('log'[2] * ' protein expression'),
         x_n_labs = TRUE) %>%
    map(~.x +
          scale_fill_manual(values = globals$cluster_colors)) %>%
    set_names(bcg_protein$top_markers)

# Plots for the ERBB - AKT - mTOR pathway -------

  insert_msg('Plot panel for ERBB signaling pathways')

  ## variables of interest

  bcg_protein$detail_vars <-
    list(receptors = c(#'KIT|c-Kit',
                       'EGFR|EGFR',
                       'EGFR|EGFR_pY1068',
                       'EGFR|EGFR_pY1173',
                       'ERBB2|HER2',
                       'ERBB2|HER2_pY1248',
                       'ERBB3|HER3',
                       'ERBB3|HER3_pY1298'),
         SRC = c('SRC|Src',
                 'SRC|Src_pY416',
                 'SRC|Src_pY527'),
         MAPK = c('NRAS|N-Ras',
                  'BRAF|B-Raf',
                  'RAF1|C-Raf',
                  'RAF1|C-Raf_pS338',
                  'MAP2K1|MEK1',
                  'MAP2K1|MEK1_pS217_S221',
                  'MAPK1|ERK2',
                  'MAPK1 MAPK3|MAPK_pT202_Y204'),
         `PI3K/AKT/mTOR` = c('PIK3CA |PI3K-p110-alpha',
                             'PIK3R1 PIK3R2|PI3K-p85',
                             'PTEN|PTEN',
                             'AKT1 AKT2 AKT3|Akt',
                             'AKT1 AKT2 AKT3|Akt_pS473',
                             'AKT1 AKT2 AKT3|Akt_pT308',
                             'TSC1|TSC1',
                             'TSC2|Tuberin',
                             'TSC2|Tuberin_pT1462',
                             'RICTOR|Rictor',
                             'RICTOR|Rictor_pT1135',
                             'MTOR|mTOR',
                             'MTOR|mTOR_pS2448'),
         apoptosis = c('BCL2L11|Bim',
                       'BAX|Bax',
                       'BID|Bid',
                       'BCL2L1|Bcl-xL',
                       'BCL2|Bcl-2',
                       'BAD|Bad_pS112')) %>%
    map(~tibble(variable = .x)) %>%
    compress(names_to = 'class') %>%
    mutate(class = factor(class,
                          c('receptors',
                            'SRC',
                            'MAPK',
                            'PI3K/AKT/mTOR',
                            'apoptosis')))

  ## n numbers of samples in the clusters to be presented
  ## in plot facets

  bcg_protein$n_numbers <- bcg_protein$data %>%
    count(clust_id)

  bcg_protein$n_numbers <-
    map2_chr(bcg_protein$n_numbers[[1]],
             bcg_protein$n_numbers[[2]],
             paste, sep = '\nn = ') %>%
    set_names(bcg_protein$n_numbers[[1]])

  ## and regulation estimates, p values, and effect sizes

  bcg_protein$detail_data <-
    list(test = bcg_protein$dev_test[c("clust_id",
                                       "variable",
                                       "deviation_center",
                                       "lower_ci",
                                       "upper_ci",
                                       "p_adjusted",
                                       "regulation")],
         anova = bcg_protein$anova[c("variable",
                                     "p_adjusted",
                                     "effect_size")] %>%
           set_names(c('variable', 'p_anova', 'eta_sq'))) %>%
    reduce(left_join, by = 'variable') %>%
    filter(variable %in% bcg_protein$detail_vars$variable) %>%
    left_join(bcg_protein$detail_vars,
              by = 'variable') %>%
    mutate(variable = factor(variable, rev(bcg_protein$detail_vars$variable)))

  ## plot panel

  bcg_protein$detail_plot <- bcg_protein$detail_data %>%
    ggplot(aes(x = deviation_center,
               y = variable,
               color = regulation,
               fill = regulation)) +
    facet_grid(class ~ clust_id,
               scales = 'free_y',
               space = 'free_y',
               labeller = labeller(.cols = bcg_protein$n_numbers)) +
    geom_vline(xintercept = 0,
               linetype = 'dashed') +
    geom_segment(aes(x = 0,
                     y = variable,
                     xend = deviation_center,
                     yend = variable)) +
    geom_point(aes(size = -log10(p_adjusted)),
               shape = 16) +
    scale_fill_manual(values = c(upregulated = 'firebrick',
                                 downregulated = 'steelblue',
                                 ns = 'gray70'),
                      name = 'Regulation\nvs cohort mean') +
    scale_color_manual(values = c(upregulated = 'firebrick',
                                  downregulated = 'steelblue',
                                  ns = 'gray70'),
                       name = 'Regulation\nvs cohort mean') +
    scale_size_area(max_size = 4,
                    name = expression('-log'[10] * ' pFDR')) +
    globals$common_theme +
    theme(axis.title.y = element_blank()) +
    labs(title = paste('Regulation of GF pathway proteins,',
                       globals$cohort_labs["tcga"]),
         x = expression('log'[2] * ' fold-regulation'))

# END ------

  bcg_protein$data <- NULL
  bcg_protein$lexicon <- NULL
  bcg_protein$top_captions <- NULL
  bcg_protein$top_markers <- NULL
  bcg_protein$detail_data <- NULL

  bcg_protein <- compact(bcg_protein)

  insert_tail()
