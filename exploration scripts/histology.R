# Differences in expression of genes of interest between the main
# histologies (normal, seminoma, NSGCT).
#
# Non-parametric tests are used for comparison of expression between
# the histologies, because some genes have a clearly bi-modal expression pattern.
# The analysis is done with Combat-corrected expression estimates for the genes
# passing variability and expression criteria

  insert_head()

# container ------

  expl_histo <- list()

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis data ------

  insert_msg('Analysis data')

  expl_histo$variables <- expl_dist$top_variables

  expl_histo$data$clinic <- globals$cohort_expr %>%
    eval %>%
    map(~.x$clinic[c('sample_id', 'histology')])

  expl_histo$data$expression <- combat$expression %>%
    map(~.x[c('sample_id', expl_histo$variables)])

  expl_histo$data <-
    map2(expl_histo$data[[1]],
         expl_histo$data[[2]],
         inner_join, by = 'sample_id') %>%
    map(~filter(.x, complete.cases(.x))) %>%
    map(mutate,
        histology = droplevels(histology)) %>%
    map(column_to_rownames, 'sample_id')

# Descriptive stats --------

  insert_msg('Descriptive stats')

  expl_histo$stats <- expl_histo$data %>%
    map(fast_num_stats,
        split_factor = 'histology')

# Tests ------

  insert_msg('Tests')

  expl_histo$test$test <- expl_histo$data %>%
    future_map(compare_variables,
               variables = expl_histo$variables,
               split_factor = 'histology',
               what = 'eff_size',
               types = 'wilcoxon_r',
               exact = FALSE,
               ci = FALSE,
               adj_method = 'BH',
               pub_styled = FALSE,
               .options = furrr_options(seed = TRUE)) %>%
    map(mutate,
        eff_size = paste('r =', signif(estimate, 2)),
        plot_cap = paste(eff_size, significance, sep = ', '))

  ## computing difference in medians as an additional
  ## effect size measure to be presented in volcano plots
  ## appending the test results with gene classification information
  ## which will be displayed in the volcano plots as well

  expl_histo$test$delta_medians <- expl_histo$data %>%
    map(blast, histology, .skip = TRUE) %>%
    map(map, colMedians) %>%
    map(~map2_dbl(.x[['NSGCT']],
                  .x[['seminoma']],
                  `-`)) %>%
    map(compress,
        names_to = 'variable',
        values_to = 'delta_median')

  expl_histo$test <-
    map2(expl_histo$test$test,
         expl_histo$test$delta_medians,
         left_join, by = 'variable') %>%
    map(mutate,
        gene_symbol = variable,
        volcano_lab = ifelse(p_adjusted < 0.05 & estimate >= 0.3,
                             gene_symbol, NA),
        regulation = ifelse(p_adjusted >= 0.05, 'ns',
                            ifelse(delta_median > 0,
                                   'upregulated',
                                   ifelse(delta_median < 0,
                                          'downregulated', 'ns'))),
        regulation = factor(regulation,
                            c('upregulated', 'downregulated', 'ns'))) %>%
    map(left_join,
        globals$gene_lexicon[, c("gene_symbol", "class")],
        by = 'gene_symbol')

  ## significance

  expl_histo$significant <- expl_histo$test %>%
    map(filter, regulation %in% c('upregulated', 'downregulated')) %>%
    map(blast, regulation) %>%
    transpose %>%
    map(map, ~.x$variable)

  ## common significant effects shared by the TCGA and GSE99420 cohorts

  expl_histo$common_significant <- expl_histo$significant %>%
    map(~.x[c("tcga", "gse99420")]) %>%
    map(reduce, intersect)

# Numbers of significantly up- and downregulated genes ------

  insert_msg('Numbers of significantly regulated genes')

  ## they will be presented in captions of volcano plots

  expl_histo$n_numbers <- expl_histo$significant %>%
    transpose %>%
    map(map_dbl, length) %>%
    map(~map2_chr(names(.x), .x, paste, sep = ': n = ')) %>%
    map(paste, collapse = ', ')

# plots for single variables -------

  insert_msg('Box plots')

  for(i in names(expl_histo$data)) {

    expl_histo$plots[[i]] <-
      list(variable = expl_histo$test[[i]]$variable,
           plot_title = expl_histo$test[[i]]$variable %>%
             html_italic %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = expl_histo$test[[i]]$plot_cap) %>%
      pmap(plot_variable,
           expl_histo$data[[i]],
           split_factor = 'histology',
           type = 'box',
           x_n_labs = TRUE,
           cust_theme = globals$common_theme,
           y_lab = expression('log'[2] * ' expression')) %>%
      map(~.x +
            scale_fill_manual(values = globals$histo_colors) +
            guides(fill = 'none') +
            theme(plot.title = element_markdown(),
                  axis.title.x = element_blank())) %>%
      set_names(expl_histo$test[[i]]$variable)

  }

# Result table ---------

  insert_msg('Result table')

  expl_histo$result_tbl <-
    map2(expl_histo$stats,
         map(expl_histo$test,
             ~.x[c('variable', 'significance', 'eff_size')]),
         left_join, by = 'variable') %>%
    map(format_tbl,
        lexicon = NULL,
        rm_complete = TRUE) %>%
    map(set_names,
        c('Variable',
          levels(expl_histo$data[[1]]$histology),
          'Significance', 'Effect size'))

# Heat maps for the common significant effects -------

  insert_msg('Heat maps')

  ## ready-to-use Y axis labels: significant effects highlighted
  ## in bold

  expl_histo$axis_labs <- expl_histo$test %>%
    map(mutate,
        #axis_lab = html_italic(variable),
        axis_lab = paste(html_italic(variable),
                         plot_cap,
                         sep = '<br>'),
        axis_lab = ifelse(p_adjusted < 0.05,
                          html_bold(axis_lab),
                          axis_lab)) %>%
    map(~set_names(.x$axis_lab, .x$variable))

  ## common gene classification, developed in the TCGA cohort

  expl_histo$classification <- expl_histo$data$tcga %>%
    classify(variables = reduce(expl_histo$common_significant, union),
             split_fct = 'histology')

  ## heat maps

  expl_histo$hm_plots <-
    list(data = expl_histo$data,
         plot_title = globals$cohort_labs[names(expl_histo$data)]) %>%
    pmap(heat_map,
         variables = reduce(expl_histo$common_significant, union),
         split_fct = 'histology',
         normalize = TRUE,
         midpoint = 0,
         limits = c(-3, 3),
         oob = scales:::squish,
         cust_theme = globals$common_theme +
           theme(axis.text.y = element_markdown(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_blank()),
         variable_classification = expl_histo$classification$classification,
         x_lab = 'Cancer sample')

  expl_histo$hm_plots <-
    list(x = expl_histo$hm_plots,
         y = expl_histo$axis_labs) %>%
    pmap(function(x, y) x +
           scale_y_discrete(labels = y))

# Volcano plots -------

  insert_msg('Volcano plots')

  expl_histo$volcano_plots <-
    list(x = expl_histo$test,
         y = globals$cohort_labs[names(expl_histo$test)],
         z = expl_histo$n_numbers) %>%
    pmap(function(x, y, z) x %>%
           ggplot(aes(x = delta_median,
                      y = -log10(p_adjusted),
                      fill = class,
                      alpha = regulation,
                      size = estimate)) +
           geom_vline(xintercept = 0,
                      linetype = 'dashed') +
           geom_hline(yintercept = -log10(0.05),
                      linetype = 'dashed') +
           geom_point(shape = 21) +
           geom_text_repel(aes(label = volcano_lab),
                           size = 2,
                           fontface = 'italic',
                           show.legend = FALSE) +
           scale_fill_manual(values = globals$gene_class_colors,
                             name = 'Gene\nclassification') +
           scale_size_area(max_size = 5,
                           limits = c(0, 1),
                           name = 'Effect size\nr') +
           scale_alpha_manual(values = c('upregulated' = 1,
                                         'downregulated' = 1,
                                         'ns' = 0.25)) +
           globals$common_theme +
           labs(title = y,
                subtitle = z,
                x = expression('log'[2] * ' fold-regulation, NSGCT vs seminoma'),
                y = expression('-log'[10] * ' pFDR')))

# END -------

  expl_histo$data <- NULL
  expl_histo$n_numbers <- NULL
  expl_histo$variables <- NULL
  expl_histo$test_types <- NULL
  expl_histo$classifiction <- NULL

  expl_histo <- compact(expl_histo)

  plan('sequential')

  rm(i)

  insert_tail()
