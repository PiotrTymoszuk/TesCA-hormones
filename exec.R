# Launches the entire pipeline

  library(soucer)

  print(source_all(c('import.R',
                     'exploration.R',
                     'clustering.R',
                     'characteristic.R',
                     'paper.R'),
                   message = TRUE, crash = TRUE))
