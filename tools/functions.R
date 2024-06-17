# Functions

# tools ------

  library(tidyverse)
  library(trafo)
  library(stringi)

  library(clustTools)
  library(igraph)
  library(ggnetwork)

  library(survival)
  library(survminer)

  library(ranger)

# formatting of tables with descriptive stats ------

  format_stats <- function(x) {

    tbl <- x %>%
      reduce(left_join, by = 'variable') %>%
      set_names(c('variable', names(x)))

  }

  format_tbl <- function(x,
                         lexicon = globals$clinic_lexicon,
                         rm_complete = FALSE, ...) {

    ## common format of the descriptive statistic tables

    x <- x %>%
      map_dfc(stri_replace,
              regex = '^Mean.*\\nMedian\\s{1}=\\s{1}',
              replacement = '') %>%
      map_dfc(stri_replace,
              fixed = 'Range',
              replacement = 'range') %>%
      map_dfc(stri_replace,
              fixed = 'Complete',
              replacement = 'complete') %>%
      map_dfc(stri_replace,
              regex = '^no:.*\\nyes:\\s{1}',
              replacement = '')

    if(rm_complete) {

      x <- x %>%
        map_dfc(stri_replace,
                regex = '\\ncomplete.*$',
                replacement = '')

    }

    if(!is.null(lexicon)) {

      x <- x %>%
        mutate(variable = exchange(variable,
                                   dict = lexicon, ...))

    }

    x

  }

  fast_num_stats <- function(data, split_factor = 'clust_id') {

    ## fast calculation of medians, interquartile ranges and ranges

    data <- data[names(data) != 'sample_id']

    stats <- data %>%
      blast(all_of(split_factor), .skip = TRUE) %>%
      map(~tibble(variable = names(.x),
                  medians = colMedians(.x),
                  lower_quart = colQuantiles(.x, 0.25),
                  upper_quart = colQuantiles(.x, 0.75),
                  min = colMins(.x),
                  max = colMax(.x))) %>%
      map(mutate,
          tbl_cell = paste0(signif(medians, 2),
                            ' [IQR: ', signif(lower_quart, 2),
                            ' - ', signif(upper_quart, 2),']\nrange: ',
                            signif(min, 2), ' - ', signif(max, 2))) %>%
      map(~.x[c('variable', 'tbl_cell')]) %>%
      reduce(left_join, by = 'variable') %>%
      set_names(c('variable', levels(data[[split_factor]])))

    ## appending with the n numbers

    n_numbers <- count(data, .data[[split_factor]]) %>%
      column_to_rownames(split_factor) %>%
      t %>%
      as_tibble %>%
      mutate(variable = 'Samples, N') %>%
      relocate(variable)

    full_rbind(n_numbers, stats)

  }

# Machine learning --------

  tune_rf <- function(data,
                      formula,
                      tune_grid, ...) {

    ## tunes a Random Forest model by minimizing the OOB prediction error

    ## construction of the models ---------

    tune_grid <- tune_grid %>%
      mutate(model_id = paste0('rf_', 1:nrow(.)))

    tune_models <- tune_grid %>%
      select(-model_id) %>%
      pmap(ranger,
           formula = formula,
           data = data, ...) %>%
      set_names(tune_grid$model_id)

    ## OOB stats ----------

    oob_errors <- tune_models %>%
      map_dbl(~.x$prediction.error)

    tune_stats <- tune_grid %>%
      mutate(oob_error = oob_errors) %>%
      relocate(model_id, oob_error) %>%
      as_tibble

    best_tune <- tune_stats %>%
      filter(oob_error == min(oob_error, na.rm = TRUE))

    best_tune <- best_tune[1, ]

    tune_stats <- tune_stats %>%
      mutate(best = ifelse(model_id == best_tune$model_id,
                           'yes', 'no'))

    list(stats = tune_stats,
         best_tune = best_tune)

  }

  plot_ranger_tuning <- function(tune_stats,
                                 plot_title = NULL,
                                 plot_subtitle = NULL,
                                 x_lab = '# variables per random tree, mtry',
                                 y_lab = 'Classification error, OOB',
                                 split_color = c(gini = 'cornflowerblue',
                                                 extratrees = 'orangered3',
                                                 hellinger = 'darkolivegreen4',
                                                 variance = 'darkolivegreen4',
                                                 maxstat = 'steelblue'),
                                 split_labels = c(gini = 'Gini index',
                                                  extratrees = 'ExtraTrees',
                                                  hellinger = 'Hellinger',
                                                  variance = 'Variance',
                                                  maxstat = 'MaxStat'),
                                 split_title = 'Split rule',
                                 line_alpha = 0.5,
                                 point_alpha = 0.75,
                                 split_by_node_size = FALSE) {


    ## plots results of tuning of a Random Forest classifier or a regressor

    ## plot subtitle --------

    if(is.null(plot_subtitle)) {

      sub_tbl <- tune_stats %>%
        filter(best == 'yes')

      min_error <- paste('Minimal error =', signif(sub_tbl$oob_error, 2))

      sub_tbl <- sub_tbl %>%
        select(-model_id, -oob_error, -best)

      tune_params <- sub_tbl %>%
        as.list %>%
        map2_chr(names(.), .,
                 paste, sep = ' = ') %>%
        paste(collapse = ', ')

      plot_subtitle <- paste(min_error, tune_params, sep = '\n')

    }

    ## plotting -------

    tune_plot <- tune_stats %>%
      ggplot(aes(x = mtry,
                 y = oob_error,
                 fill = splitrule,
                 color = splitrule))

    if(split_by_node_size) {

      tune_plot <- tune_plot +
        facet_grid(. ~ factor(min.node.size),
                   scales = 'free',
                   space = 'free',
                   labeller = as_labeller(function(x) paste('min.node =', x)))

    }

    tune_plot <- tune_plot +
      geom_path(alpha = line_alpha) +
      geom_point(shape = 16,
                 size = 2,
                 alpha = point_alpha) +
      geom_smooth(se = FALSE,
                  show.legend = FALSE)  +
      scale_color_manual(values = split_color,
                         labels = split_labels,
                         name = split_title) +
      scale_fill_manual(values = split_color,
                        labels = split_labels,
                        name = split_title) +
      globals$common_theme +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab)

    tune_plot

  }

# Co-expression graphs --------

  simil_graph <- function(data,
                          method = 'spearman',
                          corr_cutoff = 0.3,
                          abs_cutoff = FALSE,
                          mode = 'undirected',
                          weighted = TRUE,
                          attr_lexicon = globals$gene_lexicon, ...) {

    ## The function creates first a similarity matrix with
    ## correlation coefficients.
    ## In this matrix, all correlation coefficients with
    ## abs(corr_coeff) < corr_cutoff are set to zero.
    ## Finally, a graph object is created via `graph_from_adjacency_matrix()`
    ##
    ## If the arguments `class_attribute` and `color_attribute`

    ## correlation matrix: scaled to the 0 - 1 range

    corr_mtx <- cor(data, method = method)

    if(abs_cutoff) {

      corr_mtx <- ifelse(abs(corr_mtx) < corr_cutoff, 0, corr_mtx)

      corr_mtx <- 0.5 * (1 + corr_mtx)

    } else {

      corr_mtx <- ifelse(corr_mtx < corr_cutoff, 0, corr_mtx)

    }

    graph_obj <-
      graph_from_adjacency_matrix(adjmatrix = corr_mtx,
                                  mode = mode,
                                  weighted = weighted, ...)

    if(is.null(attr_lexicon)) return(graph_obj)

    attr_names <- names(attr_lexicon)[-1]

    vertice_names <- V(graph_obj)$name

    for(i in attr_names) {

      attr_values <- exchange(vertice_names,
                              attr_lexicon,
                              key = names(attr_lexicon)[1],
                              value = i) %>%
        unname %>%
        as.character

      graph_obj <- set.vertex.attribute(graph_obj,
                                        name = i,
                                        value = attr_values)

    }

    graph_obj

  }

  plot_simil_graph <- function(graph_object,
                               node_color_var = 'class',
                               weighting_order = 1,
                               point_size = 2,
                               point_rim_color = NULL,
                               txt_size = 2.75,
                               node_txt_color = 'black',
                               node_txt_face = 'italic',
                               label_edges = FALSE,
                               node_palette = globals$gene_class_colors,
                               node_title = 'Gene\nclassification',
                               edge_title = "Spearman's\n\u03C1",
                               cust_theme = theme_void(),
                               plot_title = NULL,
                               plot_subtitle = NULL,
                               select_nodes = NULL,
                               label_fun = identity, ...) {

    ## visualized a network of gene co-expression stored in a graph object
    ## node and node label color corresponds to gene classification
    ## edge width and alpha represent similarity, i.e. correlation coefficient

    ## plot labels --------

    if(is.null(plot_subtitle)) {

      node_n <- length(V(graph_object))

      edge_n <- length(E(graph_object))

      plot_subtitle <- paste0('genes: n = ', node_n,
                              ', edges: n = ', edge_n)

    }

    ## plotting -------

    graph_plot <- graph_object %>%
      ggplot(aes(x = x,
                 y = y,
                 xend = xend,
                 yend = yend)) +
      geom_edges(aes(alpha = weight^weighting_order,
                     linewidth = weight^weighting_order))

    if(!is.null(point_rim_color)) {

      graph_plot <- graph_plot +
        geom_nodes(aes(fill = .data[[node_color_var]]),
                   color = point_rim_color,
                   size = point_size)

    } else {

      graph_plot <- graph_plot +
        geom_nodes(aes(color = .data[[node_color_var]],
                       fill = .data[[node_color_var]]),
                   size = point_size)

    }

    if(is.null(select_nodes)) {

      graph_plot <- graph_plot +
        geom_nodelabel_repel(aes(label = name,
                                 fill = .data[[node_color_var]]),
                             color = node_txt_color,
                             size = txt_size,
                             fontface = node_txt_face,
                             show.legend = FALSE,
                             box.padding = unit(0.09, "lines"),
                             point.padding = unit(1e-06, "lines"), ...)

    } else {

      graph_plot <- graph_plot +
        geom_nodelabel_repel(aes(label = ifelse(name %in% select_nodes,
                                                label_fun(name), NA),
                                 fill = .data[[node_color_var]]),
                             color = node_txt_color,
                             size = txt_size,
                             fontface = node_txt_face,
                             show.legend = FALSE,
                             box.padding = unit(0.09, "lines"),
                             point.padding = unit(1e-06, "lines"), ...)

    }



    if(label_edges) {

      graph_plot <- graph_plot +
        geom_edgelabel_repel(aes(label = signif(weight, 2)),
                             size = txt_size,
                             show.legend = FALSE)

    }

    graph_plot +
      scale_fill_manual(values = node_palette,
                        name = node_title) +
      scale_color_manual(values = node_palette,
                         name = node_title) +
      scale_linewidth(limits = c(0, 1),
                      range = c(0.2, 1.2),
                      name = edge_title) +
      scale_alpha_continuous(limits = c(0, 1),
                             range = c(0.2, 1),
                             name = edge_title) +
      cust_theme +
      labs(title = plot_title,
           subtitle = plot_subtitle)

  }

  get_graph_stats <- function(graph_object, ...) {

    ## retrieves basic graph stats from a graph object

    hub_sc <- hub_score(graph_object)

    tibble(variable = V(graph_object)$name,
           class = V(graph_object)$class,
           community = V(graph_object)$comm_id,
           degree = degree(graph_object),
           betweenness = betweenness(graph_object, ...),
           hub_score = hub_sc$vector,
           einegvalue = hub_sc$value)

  }

  community_go <- function(assignment_data) {

    ## performs Go enrichment analysis for protein communities

    ## Entrez IDs for the corre

    proteins <- assignment_data %>%
      mutate(variable = stri_replace(variable,
                                     regex = '\\|.*$',
                                     replacement = '')) %>%
      blast(comm_id) %>%
      map(~.x$variable) %>%
      map(map,
          stri_split_fixed,
          pattern = ' ') %>%
      map(unlist) %>%
      map(unique) %>%
      map(unname) %>%
      map(~.x[!is.na(.x)]) %>%
      map(~.x[.x != '']) %>%
      map(~tibble(protein_symbol = .x))

    entrez_ids <- proteins %>%
      map(mutate,
          entrez_id = mapIds(org.Hs.eg.db,
                             keys = protein_symbol,
                             keytype = 'SYMBOL',
                             column = 'ENTREZID')) %>%
      map(mutate,
          entrez_id = ifelse(is.na(entrez_id) | entrez_id == '',
                             mapIds(org.Hs.eg.db,
                                    keys = protein_symbol,
                                    keytype = 'ALIAS',
                                    column = 'ENTREZID'),
                             entrez_id)) %>%
      map(~filter(.x, complete.cases(.x)))

    universe <- entrez_ids %>%
      map(~.x$entrez_id) %>%
      reduce(union)

    entrez_ids <- entrez_ids[names(entrez_ids) != 'other']

    ## enrichment analysis

    enrichment <- entrez_ids %>%
      map(~.x$entrez_id) %>%
      future_map(~GOana(.x,
                        universe = universe,
                        ontology = 'BP',
                        adj_method = 'BH'))

    enrichment

  }

# survival analysis -------

  make_survfit <- function(formula, data) {

    ## generates a list of surv_fit objects for the general differences
    ## in survival and differences in survival between cluster pairs

    clust_pairs <- combn(levels(data$clust_id), m = 2, simplify = FALSE)

    clust_pairs <- clust_pairs %>%
      set_names(map_chr(clust_pairs, paste, collapse = '|'))

    clust_data <- clust_pairs %>%
      map(~filter(data, clust_id %in% .x)) %>%
      c(list(global = data), .)

    clust_data %>%
      map(surv_fit, formula = formula)

  }

# Drug target analysis -------

  extract_append <- function(vector,
                             target_lst,
                             regex_lexicon) {

    ## looks for a regular expression in a vector of the same
    ## length as `target_list` and appends with the list elements
    ## with a string defined by names of the regex_lexicon

    helper <- function(x, y, regex, value) {

      if(all(is.na(x))) return(y)

      if(any(stri_detect(x, regex = regex))) return(c(y, value))

      return(y)

    }

    for(i in seq_along(regex_lexicon)) {

      target_lst <-
        map2(vector,
             target_lst,
             helper,
             regex = regex_lexicon[i],
             value = names(regex_lexicon)[i])

    }

    target_lst %>%
      map(function(x) if(length(x) > 1) x[!is.na(x)] else x) %>%
      map(unique)

  }

  target_extractor <- function(data,
                               target_column = 'targets',
                               pathway_column = 'moa',
                               regex_lexicon = c('DNA' = 'DNA')) {

    data %>%
      mutate(!!target_column := extract_append(.data[[pathway_column]],
                                               .data[[target_column]],
                                               regex_lexicon))


  }

  test_targets <- function(stats) {

    ## Fisher's exact test for drug target enrichment

    ## contingency matrices

    ctg_mtx <- map2(stats$n,
                    stats$n_total,
                    ~rbind(c(.x, .y - .x),
                           c(0, .y))) %>%
      set_names(stats$target)

    tests <- map(ctg_mtx, fisher.test)

    p_values <- map_dbl(tests, ~.x$p.value)

    stats %>%
      mutate(p_value = p_values) %>%
      re_adjust(method = 'BH')


  }


# varia ------

  space_evenly <- function(x) {

    ## numbers of NA and non-NA elements

    non_na_elm <- x[!is.na(x)]
    na_elm <- x[is.na(x)]

    ## numbers and lengths of NA stretches between the non-NA
    ## elements

    num_na_seg <- length(non_na_elm) - 1

    na_seg_len <- ceiling(length(na_elm)/num_na_seg)

    if(na_seg_len == 0) na_seg_len <- 1

    ## segments of consecutive non-NA elements with the appropriate
    ## numbers of NAs attached to them, merging

    out_vec <- non_na_elm[-length(non_na_elm)] %>%
      map(~c(.x, rep(NA, na_seg_len))) %>%
      reduce(c)

    out_vec <- c(out_vec, non_na_elm[length(non_na_elm)])

    ## in case the output vector is shorter as the input vector,
    ## NAs are attached to it

    if(length(out_vec) < length(x)) {

      padding_len <- length(x) - length(out_vec)

      for(i in 1:padding_len) {

        out_vec <- c(NA, out_vec)

        if(length(out_vec) == length(x)) break

        out_vec <- c(out_vec, NA)

        if(length(out_vec) == length(x)) break

      }

    }

    if(length(out_vec) > length(x)) {

      na_indexes <- 1:length(out_vec)

      na_indexes <- na_indexes[is.na(out_vec)]

      n_indexes <- sample(na_indexes,
                          size = length(na_indexes),
                          replace = FALSE)

      for(i in na_indexes) {

        out_vec <- out_vec[-i]

        if(length(out_vec) == length(x)) break

      }

    }

    out_vec

  }

  stri_capitalize <- function(x) {

    start_x <- stri_extract(x, regex = '.{1}')

    trail_x <- stri_replace(x, regex = '.{1}', replacement = '')

    paste0(toupper(start_x), trail_x)

  }

  make_ft_heat_map <- function(ft_hm,
                               sample_rug,
                               gene_rug,
                               rel_widths = c(0.94, 0.06),
                               rel_heights = c(0.92, 0.08),
                               bottom_panel_dims = c(0.305, 1, 0.11)) {

    ## stitches heat maps for gene expression
    ## with rug color-coded panels for samples and genes

    upper_panel <-
      plot_grid(ft_hm +
                  theme(strip.text.y = element_blank(),
                        strip.background.y = element_blank(),
                        legend.position = 'none',
                        axis.title.x = element_blank()),
                gene_rug +
                  theme(strip.text.y = element_blank(),
                        strip.background.y = element_blank(),
                        axis.text.x = element_blank(),
                        axis.text.y = element_blank(),
                        axis.title = element_blank(),
                        axis.ticks = element_blank(),
                        panel.grid.major = element_blank(),
                        legend.position = 'none',
                        plot.margin = ggplot2::margin()),
                ncol = 2,
                align = 'h',
                axis = 'tblr',
                rel_widths = rel_widths)

    bottom_panel <-
      plot_grid(ggdraw(),
                sample_rug +
                  theme(strip.text.x = element_blank(),
                        strip.background.x = element_blank(),
                        axis.text.x = element_blank(),
                        axis.text.y = element_blank(),
                        axis.title.y = element_blank(),
                        axis.ticks = element_blank(),
                        panel.grid.major = element_blank(),
                        legend.position = 'none',
                        plot.margin = ggplot2::margin()) +
                  labs(x = 'cancer sample'),
                ggdraw(),
                ncol = 3,
                rel_widths = bottom_panel_dims)

    plot_grid(upper_panel,
              bottom_panel,
              nrow = 2,
              rel_heights = rel_heights)

  }

# Labellers and markdown utilities -------

  subsystem_labeller <- function(x) {

    x %>%
      stri_replace(regex = '(A|a)ndrogen',
                   replacement = 'A') %>%
      stri_replace(regex = '(A|e)strogen',
                   replacement = 'E') %>%
      stri_replace(fixed = 'AE',
                   replacement = 'A/E') %>%
      stri_replace(fixed = 'synthesis and metabolism',
                   replacement = 'metabolism') %>%
      stri_replace(regex = '(M|m)ethionine',
                   replacement = 'Met') %>%
      stri_replace(regex = '(C|c)ysteine',
                   replacement = 'Cys') %>%
      stri_replace(fixed = 'golgi',
                   replacement = 'Golgi') %>%
      stri_replace(fixed = 'Oxidative phosphorylation',
                   replacement = 'OxPhos') %>%
      stri_replace(fixed = 'Fatty acid synthesis',
                   replacement = 'FASynth') %>%
      stri_replace(fixed = 'Fatty acid oxidation',
                   replacement = 'FAOx') %>%
      stri_replace_last(fixed = ' ',
                        replacement = '\n') %>%
      stri_replace(fixed = '\nmetabolism',
                   replacement = '') %>%
      stri_replace(fixed = '\ndetoxification',
                   replacement = 'detox')

  }

  metab_labeller <- function(x) {

    ## labeller for metabolic reactions

    lex <-
      c('R_CYOOm3' = 'III, CYOOm3',
        'R_CYOOm2' = 'III, CYOOm2',
        'R_NADH2_u10m' = 'I, NADH2 u10m',
        'R_ATPS4m' = 'IV, ATPS4m',

        'R_ICDHyrm' = 'IDH2, mitochondria',
        'R_ICDHxm' = 'IDH3, mitochondria',
        'R_ICDHy' = 'IDH1',
        'R_ICDHyp' = 'IDH1, peroxisome',
        'R_r0081' = 'Alanine transaminase',

        'R_HSD17B1' = 'HSD17B1, testicular',
        'R_HSD17B2r' = 'HSD17B2',
        'R_HSD17B8r' = 'HSD17B8',
        'R_DHEASULT' = 'DHEA sulfotransferase',
        'R_P45017A1r' = 'CYP17A1, pregnolone',
        'R_P45017A2r' = 'CYP17A1, DHEA',
        'R_P45017A3r' = 'CYP17A1, 17a-OH-progesterone',
        'R_P45017A4r' = 'CYP17A1, Androst-4-ene-3,17-dione',
        'R_HSD11B1r' = 'HSD11B1',
        'R_P45019A1r' = 'Aromatase, estrone',
        'R_P45019A2r' = 'Aromatase, estradiol',
        'R_P45011B11m' = 'CYP11B1, corticosterone',
        'R_P45011B12m' = 'CYP11B1, cortisol',
        'R_RE2235R' = 'Monooxygenase, 16a-OH-estrone',
        'R_RE3013R' = 'Monooxygenase, 2-OH-estradiol-17b',
        'R_HSD11B2r' = 'HSD11B2',
        'R_5ADTSTSTERONESULT' = '5a-DHT sulfotransferase',
        'R_TSTSTERONESULT' = 'testosterone sulfotransferase',
        'R_P45011A1m' = 'CYP11A1',
        'R_P4503A7r' = 'CYP3A7, testosterone',
        'R_P45021A1r' = 'CYP21A2, progesterone',
        'R_P45021A2r' = 'CYP21A2, 17a-OH-progesterone',
        'R_SR5AR2r' = 'SRD5A1/2, Androst-4-ene-3,17-dione',
        'R_SR5ARr' = 'SRD5A1/2, testosterone',
        'R_HSD17B3r' = 'HSD17B3')

    x <- ifelse(x %in% names(lex),
                lex[x],
                ifelse(is.na(annotate_bigg(x)),
                       stri_replace(x, regex = '^R_', replacement = ''),
                       annotate_bigg(x)))

    x

  }

  steroid_classifier <- function(x) {

    ## classifies steroid metabolism enzymes

    lex <-
      c('R_HSD17B2r' = 'androgens/testosterone',
        'R_HSD17B1' = 'estrogens/progesteron',
        'R_HSD17B8r' = 'estrogens/progesteron',
        'R_DHEASULT' = 'steroids',
        'R_P45017A1r' = 'estrogens/progesteron',
        'R_P45017A2r' = 'steroids',
        'R_P45017A3r' = 'steroids',
        'R_P45017A4r' = 'androgens/testosterone',
        'R_RE1096M' = 'steroids',
        'R_RE1096R' = 'steroids',
        'R_RE1134M' = 'steroids',
        'R_RE1134R' = 'steroids',
        'R_RE2768M' = 'steroids',
        'R_RE2768R' = 'steroids',
        'R_HSD11B1r' = 'corticosteroids',
        'R_P45019A1r' = 'estrogens/progesteron',
        'R_P45019A2r' = 'estrogens/progesteron',
        'R_P45011B11m' = 'corticosteroids',
        'R_P45011B12m' = 'corticosteroids',
        'R_RE2235R' = 'estrogens/progesteron',
        'R_RE3013R' = 'estrogens/progesteron',
        'R_HSD11B2r' = 'corticosteroids',
        'R_5ADTSTSTERONESULT' = 'androgens/testosterone',
        'R_TSTSTERONESULT' = 'androgens/testosterone',
        'R_P45011A1m' = 'steroids',
        'R_PRGNLONESULT' = 'steroids',
        'R_P4503A7r' = 'androgens/testosterone',
        'R_P45021A1r' = 'estrogens/progesteron',
        'R_P45021A2r' = 'steroids',
        'R_SR5AR2r' = 'androgens/testosterone',
        'R_SR5ARr' = 'androgens/testosterone',
        'R_RE1635R' = 'estrogens/progesteron',
        'R_RE1815R' = 'estrogens/progesteron',
        'R_RE1816R' = 'estrogens/progesteron',
        'R_RE1817R' = 'estrogens/progesteron',
        'R_RE2318R' = 'estrogens/progesteron',
        'R_RE2319R' = 'estrogens/progesteron',
        'R_RE1635X' = 'estrogens/progesteron',
        'R_RE1815X' = 'estrogens/progesteron',
        'R_RE1816X' = 'estrogens/progesteron',
        'R_RE1817X' = 'estrogens/progesteron',
        'R_RE2318X' = 'estrogens/progesteron',
        'R_RE2319X' = 'estrogens/progesteron',
        'R_RE1635M' = 'estrogens/progesteron',
        'R_RE1815M' = 'estrogens/progesteron',
        'R_RE1816M' = 'estrogens/progesteron',
        'R_RE1817M' = 'estrogens/progesteron',
        'R_RE2318M' = 'estrogens/progesteron',
        'R_RE2319M' = 'estrogens/progesteron',
        'R_HSD17B3r' = 'androgens/testosterone')

    factor(lex[x],
           c('steroids',
             'corticosteroids',
             'estrogens/progesteron',
             'androgens/testosterone'))

  }

  protein_labeller <-
    function(x) stri_replace(x, regex = '^.*\\|', replacement = '')

  compact_tbl <- function(x) {

    stopifnot(is.data.frame(x))

    x %>%
      map_dfc(stri_replace, fixed = ' [', replacement = '\n[') %>%
      map_dfc(stri_replace, regex = '^ns\\s{1}', replacement = 'ns\n')

  }

# END -------
