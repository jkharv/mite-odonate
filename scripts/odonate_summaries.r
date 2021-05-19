# Author: Jake Harvey -- jakekharvey@gmail.com

# Pull together stuff to make a nice odonate species level dataset on mite
# community characteristics.

library(tidyverse)
library(lubridate)
library(ape)
library(phytools)

odonate_summary <- function(mites_path){
  
  # For some reason, readr has trouble finding the column type when it begins in a bunch
  # of NAs. It guesses col_logical and then trips up when the values are not logical.
  odonates <- read_csv("datasets_derived/odonates_merged.csv",
                       col_types = cols(immune_response = col_double()))
  mites <- read_csv(mites_path)
  
  # The 2015 data was collected using a standard sampling protocol.
  abundances <- odonates %>%
    filter(year(date) == 2015) %>%
    select(genus_species) %>%
    group_by(genus_species) %>%
    tally() %>%
    rename(abundance = n) %>%
    rename(species = genus_species)
  
  masses <- odonates %>%
    drop_na(mass) %>%
    select(c(genus_species, mass)) %>%
    mutate(n = 1) %>%
    group_by(genus_species) %>%
    summarise_all(sum) %>%
    mutate(mass = mass/n) %>%
    select(c(genus_species, mass)) %>%
    rename(species = genus_species)
  
  immune_responses <- odonates %>%
    drop_na(immune_response) %>%
    select(c(genus_species, immune_response)) %>%
    group_by(genus_species) %>%
    summarise(across(everything(), mean)) %>%
    rename(species = genus_species)
  
  mite_spec <- mites %>%
    group_by(species) %>%
    select(7:last_col()) %>%
    summarise(across(.cols = everything(), sum))
  
  # Join these dataframes together
  odonates <- list(abundances, masses, mite_spec, immune_responses) %>%
    reduce(function(x, y) full_join(x, y, by = "species")) 
  
  # Add a column for suborder which will be useful later for plots.
  phylo <- read.nexus("datasets_primary/phylogeny/FinalOdonatetree.tre")
  zygoptera_mrca <- getMRCA(phylo, c("Lestes_congener", "Enallagma_hageni"))
  zygoptera <- getDescendants(phylo, zygoptera_mrca)
  zygoptera <- keep.tip(phylo, zygoptera)
  
  odonates <- mutate(odonates, suborder = if_else(species %in% zygoptera$tip.label, 
                                                  "Zygoptera", "Anisoptera"))
  
  if(grepl("otu90", mites_path)){
    write_csv(odonates, "datasets_derived/odonate_summaries_otu90.csv")
    return()
  }
  if(grepl("otu97", mites_path)){
    write_csv(odonates, "datasets_derived/odonate_summaries_otu97.csv")
    return()
  }
  if(grepl("asv", mites_path)){
    write_csv(odonates, "datasets_derived/odonate_summaries_asv.csv")
    return()
  }
}

odonate_summary("datasets_derived/mite_summaries_otu90.csv")
odonate_summary("datasets_derived/mite_summaries_otu97.csv")
odonate_summary("datasets_derived/mite_summaries_asv.csv")
