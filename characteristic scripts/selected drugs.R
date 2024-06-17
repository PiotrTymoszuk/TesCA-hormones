# Detailed visualization for selected drug classes:
#
# 1) Drugs acting via estrogen receptor.
#
# 2) Drugs acting via SRC, KIT, MAPKs, and FGFRs, MEK/ERK
#
# 3) Drugs acting via CHEK1/2, ATR, CDK2, WEE1, PLK, TOP - i.e. cell cycle
# interference
#
# 4) Apoptosis modulators acting via BCL, MCL, PARP-family members
#
# Those drugs are restricted to substances found to be differ in response
# between the clusters both in the TCGA and GSE99420 cohort

  insert_head()

# container --------

  bcg_drdetails <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  bcg_drdetails$drug_class_colors <-
    c(`sex hormones` = 'aquamarine4',
      `growth factor` = 'orangered3',
      `cell cycle/checkpoint` = 'plum4',
      apoptosis = 'gray40',
      epigenetics = 'bisque',
      `transcription/translation` = 'steelblue',
      `DNA synthesis` = 'darkslategray2',
      `adhesion/cytoskeleton` = 'goldenrod2')

# search regular expressions -------

  insert_msg('Regular expressions for search of the molecular targets')

  ## the regular expressions will be used for search within the molecular
  ## targets, moas (mechanism of actions), and pathways

  bcg_drdetails$regex <-
    list(`sex hormones` = c('ESR', 'AR'),
         `growth factor` = c('KIT', 'KDR', 'SRC', 'FLT', 'MEK',
                             'ERK', 'FGFR', 'MET', 'BRAF', 'ABL',
                             'PI3', 'PDGFR', 'MAPK', 'MEK', 'RAF\\s+',
                             'EGFR', 'ERBB', 'STAT(3|5)', 'RET', 'AKT',
                             'KRAS', 'NIK', 'NTRK2', 'ERK',
                             'VEGFR', 'YES1', 'LCK', 'EPHA2',
                             'c-FGR', 'EPHB4', 'Ephrins', 'PKC',
                             'TGFB1', 'BMP', '(m|M)TOR', 'IGF', 'HER\\d{1}',
                             'JAK\\d{1}', 'SYK', 'LCK', 'FYN',
                             'RTK\\s{1}signaling'),
         `cell cycle/checkpoint` = c('CDK', 'CHEK', 'ATR', 'WEE', 'PLK',
                                     'WRN', 'NPM1', 'CDC25', 'UAF1', 'USP1',
                                     'AURK', 'MDM\\d{1}', 'p53', 'ATM',
                                     'GADD34', 'P53', 'TP53', 'KIF11', 'GADD34',
                                     '(C|c)ell\\s{1}cycle', '(M|m)itosis'),
         apoptosis = c('BCL', 'MCL', 'PARP', 'PPM1D', 'CASP3', 'Bax', 'BIRC',
                       'IAP', 'SMAC', 'DIABLO',
                       'FAS', 'PIM3', 'BFL1', 'DAPK3', 'WIP1', 'Bcl', 'TRAIL'),
         epigenetics = c('HDAC', 'PORCN', 'DOT1L', 'EP300', 'G9A', 'LSD1',
                         'BRD\\d+', 'EHMT\\d+', 'Polybromo', 'SMARCA\\d+',
                         '(C|c)hromatin', '(D|d)emethyl'),
         `transcription/translation` = c('XPO1', 'RNA', 'EIF\\d+', 'eEF',
                                         'CLK4', 'TAF1', 'eIF\\d+', 'EF\\d+'),
         `DNA synthesis` = c('DHFR', 'TERT', 'TOP',
                             'Anti-metabolite', 'Antimetabolite',
                             'DNA', 'Mutant', 'Tankyrase', 'Telomerase',
                             'TYMS', '(g|Genome)', 'DHFR'),
         `adhesion/cytoskeleton` = c('KSP', 'tubule', 'MRCKB',
                                     'PPK', 'RAC\\d+', 'ROCK2',
                                     '(C|c)ytoskeleton', 'centrosome')) %>%
    map_chr(paste, collapse = '|') %>%
    map_chr(~paste0('^(', .x, ')'))

# Drugs of interest ------

  insert_msg('Drugs of interest')

  ## common significant drugs, significant in the TCGA and GSE99420
  ## cohorts

  bcg_drdetails$common_drugs <- bcg_drugs$common_significant %>%
    map(unlist) %>%
    map(unname) %>%
    map(unique)

  ## selection of drugs for the targets of interest, that were found significant
  ## in the TCGA and GSE99420 cohort

  bcg_drdetails$lexicon <- drugs$lexicons %>%
    map(select, variable, drug_name, targets,
        any_of(c('moa', 'pathway_name')))

  for(i in names(bcg_drdetails$lexicon)) {

    bcg_drdetails$variables[[i]]$targets <- bcg_drdetails$regex %>%
      map(function(rex) bcg_drdetails$lexicon[[i]] %>%
            filter(map_lgl(targets, ~any(stri_detect(.x, regex = rex)))))

    bcg_drdetails$variables[[i]]$moa <- bcg_drdetails$regex %>%
      map(~reglook(select(bcg_drdetails$lexicon[[i]],
                          variable, any_of(c('moa', 'pathway_name'))),
                   regex = .x))

    bcg_drdetails$variables[[i]] <- bcg_drdetails$variables[[i]] %>%
      map(map, ~.x[c('variable')]) %>%
      map(map,
          filter, variable %in% bcg_drdetails$common_drugs[[i]]) %>%
      map_dfr(compress,
              names_to = 'drug_class') %>%
      mutate(drug_class = factor(drug_class, names(bcg_drdetails$regex))) %>%
      filter(!duplicated(variable))

  }

# ANOVA results to be presented in box plots and heat maps --------

  insert_msg('ANOVA results to be shown in the plots')

  for(i in names(bcg_drdetails$variables)) {

    bcg_drdetails$anova[[i]] <- bcg_drugs$anova[[i]] %>%
      map(~left_join(bcg_drdetails$variables[[i]],
                     .x,
                     by = 'variable')) %>%
      map(re_adjust, method = 'none') %>%
      map(mutate,
          eff_size = paste('\u03B7\u00B2 =', signif(effect_size, 2)),
          plot_cap = paste(eff_size, significance, sep = ', ')) %>%
      map(select,
          drug_class,
          variable,
          drug_name,
          plot_cap)

  }

# Drug response predictions for the compounds of interest --------

  insert_msg('Drug response predictions')

  for(i in names(bcg_drdetails$variables)) {

    bcg_drdetails$data[[i]] <- drugs$predictions[[i]] %>%
      map(select, sample_id, all_of(bcg_drdetails$variables[[i]]$variable)) %>%
      map2(bcg_globals$assignment, .,
           inner_join,
           by = 'sample_id')

  }

# Box plots for single drugs -------

  insert_msg('Box plots for single drugs')

  for(i in names(bcg_drdetails$anova)) {

    for(j in names(bcg_drdetails$anova[[i]])) {

      bcg_drdetails$plots[[i]][[j]] <-
        list(variable = bcg_drdetails$anova[[i]][[j]]$variable,
             plot_title = bcg_drdetails$anova[[i]][[j]]$drug_name %>%
               paste(globals$cohort_labs[j], sep = ', '),
             plot_subtitle = bcg_drdetails$anova[[i]][[j]]$plot_cap) %>%
        pmap(plot_variable,
             bcg_drdetails$data[[i]][[j]],
             split_factor = 'clust_id',
             type = 'box',
             cust_theme = globals$common_theme,
             x_lab = 'cluster',
             y_lab = paste(globals$drug_exp_labs[i],
                           globals$drug_unit_labs[i],
                           sep = '-trained predictions, '),
             x_n_labs = TRUE) %>%
        map(~.x +
              scale_fill_manual(values = globals$cluster_colors)) %>%
        set_names(bcg_drdetails$anova[[i]][[j]]$variable)

    }

  }

# Classification of drugs by cluster specificity -------

  insert_msg('Drug classification in the TCGA cohort')

  bcg_drdetails$classification_tcga <-
    list(data = map(bcg_drdetails$data, ~.x$tcga),
         variables = map(bcg_drdetails$variables, ~.x$variable)) %>%
    pmap(classify,
         split_fct = 'clust_id') %>%
    map(~.x$classification) %>%
    map2(.,
         map(bcg_drdetails$variables,
             ~.x[c('variable', 'drug_class')]),
         left_join, by = 'variable')

# Heat map plots -------

  insert_msg('Heat map plots')

  for(i in names(bcg_drdetails$data)) {

    bcg_drdetails$hm_plots[[i]] <-
      list(data = bcg_drdetails$data[[i]],
           plot_title = globals$drug_exp_labs[i] %>%
             paste(globals$cohort_labs[names(bcg_drdetails$data[[i]])],
                   sep = '-trained predictions, ')) %>%
      pmap(heat_map,
           split_fct = 'clust_id',
           variables = bcg_drdetails$variables[[i]]$variable,
           normalize = TRUE,
           variable_classification = bcg_drdetails$classification_tcga[[i]],
           cust_theme = globals$common_theme +
             theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.y = element_blank()),
           midpoint = 0,
           limits = c(-3, 3),
           oob = scales::squish,
           x_lab = 'cancer sample')

  }

  ## Y axis labels

  bcg_drdetails$hm_plots$ctrp2 <-
    bcg_drdetails$hm_plots$ctrp2 %>%
    map(~.x +
          scale_y_discrete(labels = function(x) exchange(x,
                                                         drugs$lexicons$ctrp2,
                                                         value = 'drug_name')))

  bcg_drdetails$hm_plots$gdsc <-
    bcg_drdetails$hm_plots$gdsc %>%
    map(~.x +
          scale_y_discrete(labels = function(x) exchange(x,
                                                         drugs$lexicons$gdsc,
                                                         value = 'drug_name')))

# Rug heat maps coding for the drug class -------

  insert_msg('Rug heat maps with the drug class')

  bcg_drdetails$rug_hm_plots <- bcg_drdetails$classification_tcga %>%
    map(~ggplot(.x,
                aes(x = 'class',
                    y = reorder(variable, delta_auc),
                    fill = drug_class)) +
          facet_grid(clust_id ~ .,
                     scales = 'free',
                     space = 'free') +
          geom_tile(color = 'black') +
          scale_fill_manual(values = bcg_drdetails$drug_class_colors,
                            labels = function(x) stri_replace(x,
                                                              fixed = '/',
                                                              replacement = '\n'),
                            name = 'Drug\nclassification') +
          globals$common_theme +
          theme(axis.title = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.line = element_blank(),
                axis.ticks = element_blank(),
                panel.grid.major = element_blank()))

# END -------

  rm(i, j)

  insert_tail()
