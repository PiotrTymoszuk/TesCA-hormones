# Detailed analysis of differential gene expression for histological subtypes
# defined by ICD. Teratomas are merged.
# The analysis is done for the TCGA and GSE3218 cohorts with ComBat-corrected
# expression estimates for the genes passing variability and minimal expression
# crietria.

  insert_head()

# container -----

  expl_icd <- list()

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis data ------

  insert_msg('Analysis data')

  expl_icd$variables <- expl_dist$top_variables

  expl_icd$data$clinic <- list(tcga = tcga, gse3218 = gse3218) %>%
    map(~.x$clinic) %>%
    map(filter, tissue == 'tumor') %>%
    map(~.x[c('sample_id', 'histology_icd')])

  expl_icd$data$expression <- combat$expression[c("tcga", "gse3218")] %>%
    map(~.x[c('sample_id', expl_icd$variables)])

  expl_icd$data <-
    map2(expl_icd$data[[1]],
         expl_icd$data[[2]],
         inner_join, by = 'sample_id') %>%
    map(~filter(.x, complete.cases(.x))) %>%
    map(mutate,
        histology_icd = droplevels(histology_icd))

  ## collapsing the teratomas

  expl_icd$data <- expl_icd$data %>%
    map(mutate,
        histology_icd = fct_collapse(histology_icd,
                                     TT = c('benign teratoma',
                                            'malignant teratoma',
                                            'teratocarcinoma'),
                                     TYST = c('yolk sac cancer'),
                                     MGCT = c('germinal mixed histology'),
                                     EMBCA = c('embryonal carcinoma'),
                                     SEM = c('seminoma')),
        histology_icd = fct_relevel(histology_icd,
                                    'SEM',
                                    'MGCT',
                                    'EMBCA',
                                    'TT',
                                    'TYST'))

# Descriptive stats -------

  insert_msg('Descriptive stats')

  expl_icd$stats <- expl_icd$data %>%
    map(fast_num_stats,
        split_factor = 'histology_icd')

# Tests --------

  insert_msg('Tests')

  expl_icd$test <- expl_icd$data %>%
    future_map(compare_variables,
               variables = expl_icd$variables,
               split_factor = 'histology_icd',
               what = 'eff_size',
               types = 'kruskal_etasq',
               exact = FALSE,
               ci = FALSE,
               pub_styled = TRUE,
               .options = furrr_options(seed = TRUE)) %>%
    map(mutate,
        plot_cap = paste(eff_size, significance, sep = ', '))

  ## significant effects

  expl_icd$significant <- expl_icd$test %>%
    map(filter, p_adjusted < 0.05) %>%
    map(~.x$variable)

  expl_icd$common_significant <- expl_icd$significant %>%
    reduce(intersect)

# Box plots for single variables --------

  insert_msg('Box plots for single variables')

  for(i in names(expl_icd$data)) {

    expl_icd$plots[[i]] <-
      list(variable = expl_icd$test[[i]]$variable,
           plot_title = expl_icd$test[[i]]$variable %>%
             html_italic %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = expl_icd$test[[i]]$plot_cap) %>%
      pmap(plot_variable,
           expl_icd$data[[i]],
           split_factor = 'histology_icd',
           type = 'box',
           cust_theme = globals$common_theme,
           x_n_labs = TRUE,
           y_lab = expression('log'[2] * ' expression')) %>%
      map(~.x +
            scale_fill_manual(values = globals$histo_icd_colors) +
            guides(fill = 'none') +
            theme(plot.title = element_markdown())) %>%
      set_names(expl_icd$test[[i]]$variable)

  }

# Result table -------

  insert_msg('Result table')

  expl_icd$result_tbl <-
    map2(expl_icd$stats,
         map(expl_icd$test,
             ~.x[c('variable', 'significance', 'eff_size')]),
         left_join, by = 'variable') %>%
    map(format_tbl,
        lexicon = NULL,
        rm_complete = TRUE) %>%
    map2(., expl_icd$data,
         ~set_names(.x,
                    'Variable',
                    levels(.y$histology_icd),
                    'Significance',
                    'Effect size'))

# Heat maps for the regulated genes -------

  insert_msg('Heat maps')

  ## gene classification, developed in the TCGA cohort
  ## and appended with the gene classification used later for rug heat maps

  expl_icd$classification <- expl_icd$data$tcga %>%
    classify(variables = expl_icd$significant$tcga,
             split_fct = 'histology_icd') %>%
    .$classification

  expl_icd$classification <- expl_icd$classification %>%
    mutate(gene_symbol = variable) %>%
    left_join(globals$gene_lexicon[c('gene_symbol', 'class')],
              by = 'gene_symbol')

  ## heat maps

  expl_icd$hm_plots <-
    list(data = expl_icd$data,
         plot_title = globals$cohort_labs[names(expl_icd$data)]) %>%
    pmap(heat_map,
         variables = expl_icd$significant$tcga,
         split_fct = 'histology_icd',
         normalize = TRUE,
         midpoint = 0,
         limits = c(-3, 3),
         oob = scales:::squish,
         cust_theme = globals$common_theme +
           theme(axis.text.y = element_text(face = 'italic'),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_blank()),
         variable_classification = expl_icd$classification,
         x_lab = 'cancer sample')

# Rug heat map with the color coded gene classification -------

  insert_msg('Rug heat map')

  expl_icd$rug_hm_plot <- expl_icd$classification %>%
    ggplot(aes(x = 'class',
               y = reorder(variable, delta_auc),
               fill = class)) +
    geom_tile(color = 'black') +
    facet_grid(histology_icd ~ .,
               scales = 'free',
               space = 'free') +
    scale_fill_manual(values = globals$gene_class_colors,
                      name = 'Gene\nclassification') +
    theme_void() +
    theme(strip.text = element_blank(),
          legend.text = globals$common_text,
          legend.title = globals$common_text)


# END ------

  rm(i)

  plan('sequential')

  expl_icd$data <- NULL
  expl_icd$n_numbers <- NULL
  expl_icd$variables <- NULL

  expl_icd <- compact(expl_icd)

  insert_tail()
