# Renders the manuscript and supplements

  insert_head()

# Bibliography -----

  insert_msg('Reading the bibliography')

  tesca_bib <- read_bib('./paper/markdown/tesca.bib') %>%
    as_mdbib

# Rendering the manuscript parts ------

  insert_msg('Rendering the manuscript parts')

  c('./paper/markdown/figures_tables_methods.Rmd',
    './paper/markdown/supplementary_material.Rmd') %>%
    walk(render,
         output_format = word_document2(number_sections = FALSE,
                                        reference_docx = 'ms_template.docx'),
         output_dir = './paper')

# END ----

  insert_tail()
