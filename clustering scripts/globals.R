# Globals of the clustering analysis


  insert_head()

# container ------

  clust_globals <- list()

# data -------

  insert_msg('Analysis data')

  clust_globals$variables <- expl_dist$top_variables

  clust_globals$data <- combat$expression %>%
    map(select,
        sample_id, all_of(clust_globals$variables)) %>%
    map(column_to_rownames, 'sample_id')

# END -----

  insert_tail()
