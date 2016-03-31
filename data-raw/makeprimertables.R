#load the degeneracy tables
source("degeneracy_tables.R")

#save them as package data
devtools::use_data(degenS,  overwrite = TRUE)
devtools::use_data(degenSZ, overwrite = TRUE)
devtools::use_data(degenZ,  overwrite = TRUE)
