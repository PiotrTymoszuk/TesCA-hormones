# Plots for the results of clustering of the common significantly
# enriched GO terms

  insert_head()

# container ------

  bcg_gocplots <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  ## clustering objects and cluster assignment

  bcg_gocplots[c("clust_objects", "assignment")] <-
    bcg_goclust[c("clust_objects", "assignment")]

  ## descriptions of the clusters

  bcg_gocplots$description[['#1']] <-
    list('1' = 'RNA processing',
         '2' = 'organ development',
         '3' = c('reproduction',
                  'proliferation',
                  'ECM',
                  'mesenchymal cell'))

  bcg_gocplots$description[['#2']] <-
    list('1' = c('collagen, ECM',
                 'LDL/HDL',
                 'TGFB/BMP/WNT signaling',
                 'adhesion, juctions',
                 'motility'),
         '2' = 'transcription',
         '3' = c('organ development',
                 'angiogenesis'),
         '4' = c('epithelium',
                 'connective tissue',
                 'neuronal development'))

  bcg_gocplots$description[['#3']] <-
    list('1' = c('organ development',
                 'sex differentiation'),
         '2' = c('mesenchymal cells',
                 'hemostasis/coagulation',
                 'BMP/SMAD/TGFB signaling',
                 'NOTCH/WNT signaling',
                 'FGFR signaling',
                 'wound healing',
                 'adhesion,junctions, motility'),
         '3' = c('RNA metabolism'))

  bcg_gocplots$description <- bcg_gocplots$description %>%
    map(map, paste, collapse = '\n')

# Plots of the MDS spaces -------

  insert_msg('Plots of the MDS spaces')

  bcg_gocplots$mds_plots <-
    list(x = bcg_gocplots$clust_objects,
         y = bcg_gocplots$description,
         w = paste('Common GO terms, cluster',
                   names(bcg_gocplots$clust_objects)),
         z = map(bcg_gocplots$clust_objects, nobs)) %>%
    pmap(function(x, y, w, z) x %>%
           plot(type = 'data',
                cust_theme = globals$common_theme +
                  theme(plot.tag = element_blank())) +
           scale_fill_viridis_d(labels = y) +
           labs(title = w,
                subtitle = paste('GO terms: n =', z$observations)))

# END ------

  bcg_gocplots$clust_objects <- NULL
  bcg_gocplots$assignment <- NULL

  bcg_gocplots <- compact(bcg_gocplots)

  insert_tail()
