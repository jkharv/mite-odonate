#Author: Jake Harvey (jakekharvey@gmail.com)

library(picante)
library(tidyverse)

calc_indices <- function(file_path){  
  
  bin_network <- read_csv(file_path)
  
  # Filtering out species that we don't have in our phylogeny.
  bin_network <- bin_network %>%
    filter(odonate_spp != "Enallagma_annexum") %>%
    filter(odonate_spp != "Lestes_forcipatus") %>%
    filter(odonate_spp != "Lestes_inaequalis") %>%
    filter(odonate_spp != "Lestes_vigilax")
  
  phylo_path <- "datasets_primary/phylogeny/FinalOdonatetree.tre"
  phylo <- read.nexus(phylo_path)
  
  #~~~~~~~~~~~~~~~~~~Resource Range~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Number of links normalized by number of possible links.
  # Non-phylogenetic specialization index.
  resource_range <- function(bin_network, mite){
    
    mite <- pull(bin_network, mite)
    R <- length(mite)
    r <- sum(mite) # Will only work on the binary network.
    (R - r)/(R - 1)
  }
  
  partial_rr <- Vectorize(function(m){resource_range(bin_network, m)})
  
  mites <- tibble(mite = colnames(bin_network)[-1]) %>%
    mutate(resource_range = partial_rr(mite))
  
  #~~~~~~~~~~~~~~~~~~phylo Specialization index~~~~~~~~~~~~~~~~~
  
  x <- bin_network %>%
    pivot_longer(2:last_col()) %>%
    rename(mite = name) %>%
    arrange(mite) %>%
    relocate(mite, value, odonate_spp) %>%
    sample2matrix()
  
  simpson <- as.data.frame(raoD(x, phylo)$Dkk) %>%
    rownames_to_column() %>%
    rename(simpson = "raoD(x, phylo)$Dkk", mite = rowname)
  
  faiths <- pd(x, phylo) %>%
    rownames_to_column() %>%
    rename(faiths = "PD", mite = rowname, num_host= SR)
  
  div <- full_join(faiths, simpson, by = "mite")
  mites <- full_join(mites, div, by = "mite")
  
  # File names are getting out of control.
  if(grepl("otu90", file_path)){
    write_csv(mites, "datasets_derived/mite_phylo_scale_otu90.csv")
    return()
  }
  if(grepl("otu97", file_path)){
    write_csv(mites, "datasets_derived/mite_phylo_scale_otu97.csv")
    return()
  }
  write_csv(mites, "datasets_derived/mite_phylo_scale_asv.csv")
}

calc_indices("datasets_derived/mite_sequences_annotated_bin_network.csv")
calc_indices("datasets_derived/mite_sequences_otu90_annotated_bin_network.csv")
calc_indices("datasets_derived/mite_sequences_otu97_annotated_bin_network.csv")
