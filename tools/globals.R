# Project globals

# container -----

  globals <- list()

# tools --------

  library(readxl)

# graphic globals ------

  globals$common_text <- element_text(size = 8,
                                      face = 'plain',
                                      color = 'black')

  globals$common_margin <- ggplot2::margin(t = 5,
                                           l = 4,
                                           r = 2,
                                           unit = 'mm')

  globals$common_theme <- theme_classic() +
    theme(axis.text = globals$common_text,
          axis.title = globals$common_text,
          plot.title = element_text(size = 8,
                                    face = 'bold'),
          plot.subtitle = globals$common_text,
          plot.tag = element_text(size = 8,
                                  face = 'plain',
                                  color = 'black',
                                  hjust = 0,
                                  vjust = 1),
          plot.tag.position = 'bottom',
          legend.text = globals$common_text,
          legend.title = globals$common_text,
          strip.text = globals$common_text,
          strip.background = element_rect(fill = 'gray95',
                                          color = 'gray80'),
          plot.margin = globals$common_margin,
          panel.grid.major = element_line(color = 'gray90'))

# Genes of interest -------

  globals$gene_lexicon <-
    c('GNRH1' = 'pituitary',
      'GNRH2' = 'pituitary',
      'PRL' = 'pituitary',
      'CGA' = 'pituitary',
      'FSHB' = 'pituitary',
      'LHB' = 'pituitary',
      'POMC' = 'pituitary',

      'STAR' = 'steroid',
      'STARD3' = 'steroid',
      'STARD3NL' = 'steroid',
      'STARD4' = 'steroid',
      #'STARD6' = 'steroid', ## not available for all cohorts
      'TSPO' = 'steroid',
      'TSPOAP1' = 'steroid',
      'CYP11A1' = 'steroid',
      'CYP17A1' = 'steroid',
      'FDX1' = 'steroid',
      'FDX2' = 'steroid',
      'FDXR' = 'steroid',
      'HSD3B1' = 'steroid',
      'HSD3B2' = 'steroid',
      'SERPINA6' = 'steroid',

      'CYP11B1' = 'adrenal',
      'CYP11B2' = 'adrenal',
      'CYP21A2' = 'adrenal',
      'HSD11B1' = 'adrenal',
      'HSD11B2' = 'adrenal',


      'HSD17B1' = 'gonadal',
      'HSD17B2' = 'gonadal',
      'HSD17B3' = 'gonadal',
      'HSD17B11' = 'gonadal',
      'HSD17B12' = 'gonadal',
      'HSD17B14' = 'gonadal',
      'CYP19A1' = 'gonadal',
      'SRD5A1' = 'gonadal',
      'SRD5A2' = 'gonadal',
      'SRD5A3' = 'gonadal',
      'SHBG' = 'gonadal') %>%
    compress(names_to = 'gene_symbol',
             values_to = 'class') %>%
    mutate(class = factor(class,
                          c('pituitary', 'steroid', 'adrenal', 'gonadal')),
           color = fct_recode(class,
                              gray60 = 'pituitary',
                              orangered3 = 'steroid',
                              steelblue = 'adrenal',
                              aquamarine3 = 'gonadal'),
           color = as.character(color))

  globals$genes <- globals$gene_lexicon$gene_symbol

  globals$gene_class_colors <- globals$gene_lexicon %>%
    filter(!duplicated(class))

  globals$gene_class_colors <-
    set_names(globals$gene_class_colors$color,
              as.character(globals$gene_class_colors$class))

# Cohort labels, colors and expressions -------

  globals$analysis_cohorts <- c('tcga', 'gse3218', 'gse99420')

  globals$cohort_expr <- globals$analysis_cohorts %>%
    map_chr(~paste(.x, .x, sep = ' = ')) %>%
    paste(collapse = ', ') %>%
    paste0('list(', ., ')') %>%
    parse_expr

  globals$cohort_labs <-
    c('tcga' = 'TCGA',
      'gse3218' = 'GSE3218',
      'gse99420' = 'GSE99420')

  globals$cohort_colors <-
    c('tcga' = 'indianred3',
      'gse3218' = 'steelblue',
      'gse99420' = 'steelblue')

# Tissue histology colors -------

  globals$histo_colors <- c('normal' = 'cornsilk',
                            'seminoma' = 'darkseagreen',
                            'NSGCT' = 'firebrick',
                            'not assigned' = 'gray80')

  globals$histo_icd_colors <- c('seminoma' = 'darkseagreen',
                                'SEM' = 'darkseagreen',
                                'germinal mixed histology' = 'firebrick3',
                                'MGCT' = 'firebrick3',
                                'embryonal carcinoma' = 'plum4',
                                'EMBCA' = 'plum4',
                                'teratoma' = 'gray60',
                                'TT' = 'gray60',
                                'yolk sac cancer' = 'bisque2',
                                'TYST' = 'bisque2',
                                'choriocarcinoma' = 'darkolivegreen4',
                                'TCCA' = 'darkolivegreen4')

# clinical variable lexicon --------

  globals$clinic_lexicon <- read_xlsx('./data/variable_lexicon.xlsx')

# Cluster colors and labels ------

  globals$cluster_colors <-
    c('#1' = 'gray60',
      '#2' = 'darkolivegreen4',
      '#3' =  'orangered3',
      '#4' = 'steelblue')

  globals$cluster_hex_colors <-
    c('#1' = '#999999',
      '#2' = '#6e8b3d',
      '#3' =  '#cd3700',
      '#4' = '#4682b4')

# labels for the drug response experiments --------

  globals$drug_exp_labs <- c(ctrp2 = 'CTRP2',
                             gdsc = 'GDSC1/2')

  globals$drug_unit_labs <- c(ctrp2 = 'AUC',
                              gdsc = 'log IC50 [µM]')

# Protein category colors ------

  globals$protein_colors <-
    c('cell cycle' = 'aquamarine4',
      'DNA repair/apoptosis' = 'gray60',
      'ECM/adhesion/cytoskeleton' = 'hotpink3',
      'GF receptors' = 'orangered3',
      'inflammation' = 'brown1',
      'MAPK signaling' = 'steelblue',
      'metabolism/energy' = 'orange1',
      'non-receptor kinases' = 'plum4',
      'other' = 'cornsilk',
      'other signaling' = 'bisque3',
      'PI3K/AKT/mTOR' = 'darkslategrey',
      'RNA turnover and translation' = 'gold',
      'sex hormones' = 'lightcoral',
      'transcription/chromatin' = 'darkorange3',
      'WNT/NOTCH signaling' = 'cyan4')

# END -------
