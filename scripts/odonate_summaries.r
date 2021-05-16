# Author: Jake Harvey -- jakekharvey@gmail.com

# Pull together stuff to make a nice odonate species level dataset on mite
# community characteristics.

library(tidyverse)
library(lubridate)
library(ape)
library(phytools)

odonate_summary <- function(mites_path){
  
  odonates <- read_csv("datasets_derived/odonates_merged.csv")
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
  
  mite_spec <- mites %>%
    group_by(species) %>%
    select(7:last_col()) %>%
    summarise(across(.cols = everything(), sum))
  
  # Join these dataframes together
  odonates <- list(abundances, masses, mite_spec) %>%
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
