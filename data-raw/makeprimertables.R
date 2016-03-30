#load the degeneracy tables
source("degeneracy_tables.R")

#save them as package data
devtools::use_data(degenS)
devtools::use_data(degenSZ)
devtools::use_data(degenZ)
