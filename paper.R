# Scripts for generation of paper figures, tables, supplements, and
# rendering of the manuscript parts and the supplementary material.

# tools --------

  library(tidyverse)
  library(rlang)
  library(trafo)
  library(stringi)

  library(ggrepel)
  library(figur)
  library(ggtext)

  library(AnnotationDbi)
  library(org.Hs.eg.db)

  library(cowplot)
  library(figur)

  library(flextable)
  library(rmarkdown)
  library(bookdown)
  library(knitr)

  library(soucer)

  insert_head()

  select <- dplyr::select
  reduce <- purrr::reduce
  explore <- exda::explore
  set_rownames <- trafo::set_rownames
  map <- purrr::map
  rename <- dplyr::rename
  components <- generics::components
  transpose <- purrr::transpose

  c('./tools/globals.R',
    './tools/functions.R') %>%
    source_all(message = TRUE, crash = TRUE)

# figures and tables ---------

  insert_msg('Figures and tables')

  c('./paper scripts/figures.R',
    './paper scripts/supplementary_figures.R',
    './paper scripts/tables.R',
    './paper scripts/supplementary_tables.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Rendering manuscript parts -------

  insert_msg('Rendering manuscript parts')

  c('./paper scripts/links.R',
    './paper scripts/render.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END -----

  insert_tail()
