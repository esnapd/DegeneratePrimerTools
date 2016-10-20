library(DegeneratePrimerTools)
library(dplyr)

ahbafna <- system.file("sequences/AHBA.fna",package="DegeneratePrimerTools")
AHBAseqs <- Biostrings::readDNAStringSet(ahbafna)

# run alignment, create tree, run DEGEPRIME on the MSA
AHBA_degeprimer <- degeprimer(AHBAseqs) %>%
  run_alignment() %>%
  build_tree() %>%
  design_primers( maxdegeneracies=c(1, 10, 20), number_iterations=2) %>%
  add_primerpair(name="primerpair1", fpos=455, fdeg=10, rpos=617, rdeg=10) %>%
  add_primerpair(name="primerpair2", fpos=928, fdeg=10, rpos=1133, rdeg=10)


AHBA_degeprimer

plot_GC(AHBA_degeprimer)

plot_msa(AHBA_degeprimer@msa)

msa2 <- add_primers_to_MSA(AHBA_degeprimer, max.mismatch = 3)

plot_msa(msa2)

devtools::use_data(AHBA_degeprimer)
