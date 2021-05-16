# Author: Jacob Harvey - jakekharvey@gmail.com
# Summaries from the mite pov.

library(tidyverse)

mite_summary <- function(network_path, indices_path){
  
  # network <- read_csv("datasets_derived/mite_sequences_annotated_bin_network.csv")
  # mites <- read_csv("datasets_derived/mite_phylo_scale_asv.csv")

  network <- read_csv(network_path)
  mites <- read_csv(indices_path)
  
  # Get mite statistics at the host-species level.
  get_rr <- Vectorize(function(mite){mites$resource_range[which(mites$mite == mite)]})
  get_hn <- Vectorize(function(mite){mites$num_host[which(mites$mite == mite)]})
  get_faith <- Vectorize(function(mite){mites$faiths[which(mites$mite == mite)]})
  get_simp <- Vectorize(function(mite){mites$simpson[which(mites$mite == mite)]})
  
  mite_interactions <- network %>% 
    pivot_longer(2:last_col()) %>%
    filter(value > 0) %>%
    rename(mite = name) %>%
    select(odonate_spp, mite) %>%
    mutate(resource_range = get_rr(mite)) %>%
    mutate(num_host = get_hn(mite)) %>%
    mutate(faiths = get_faith(mite)) %>%
    mutate(simpson = get_simp(mite)) %>%
    select(odonate_spp, mite, resource_range, num_host,
           faiths, simpson) %>%
    rename(species = odonate_spp)
  
  mite_interactions <- mite_interactions %>%
    # Spec/gen based on phylo scale.
    mutate(specialist    = ifelse(simpson <= 100, 1, 0)) %>%
    mutate(generalist    = ifelse(simpson >  100, 1, 0)) %>%
    mutate(specialist_50 = ifelse(simpson <= 50, 1, 0)) %>%
    mutate(generalist_50 = ifelse(simpson >  50, 1, 0)) %>%
    mutate(specialist_25 = ifelse(simpson <= 25, 1, 0)) %>%
    mutate(generalist_25 = ifelse(simpson >  25, 1, 0)) %>%
    # Spec/gen based on nhost.  
    mutate(specialist_nhost    = ifelse(num_host <= 2, 1, 0)) %>%
    mutate(generalist_nhost    = ifelse(num_host >  2, 1, 0)) %>%
    mutate(specialist_nhost_5  = ifelse(num_host <= 5, 1, 0)) %>%
    mutate(generalist_nhost_5  = ifelse(num_host >  5, 1, 0)) %>%
    mutate(specialist_nhost_10 = ifelse(num_host <= 10, 1, 0)) %>%
    mutate(generalist_nhost_10 = ifelse(num_host >  10, 1, 0))
  
  if(grepl("otu90", indices_path)){
    write_csv(mite_interactions, "datasets_derived/mite_summaries_otu90.csv")
    return()
  }
  if(grepl("otu97", indices_path)){
    write_csv(mite_interactions, "datasets_derived/mite_summaries_otu97.csv")
    return()
  }
  if(grepl("asv", indices_path)){
    write_csv(mite_interactions, "datasets_derived/mite_summaries_asv.csv")
    return()
  }
}

mite_summary("datasets_derived/mite_sequences_annotated_bin_network.csv",
             "datasets_derived/mite_phylo_scale_asv.csv")
mite_summary("datasets_derived/mite_sequences_otu90_annotated_bin_network.csv",
             "datasets_derived/mite_phylo_scale_otu90.csv")
mite_summary("datasets_derived/mite_sequences_otu97_annotated_bin_network.csv",
             "datasets_derived/mite_phylo_scale_otu97.csv")
