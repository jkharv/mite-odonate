# Author: Jake Harvey -- jakekharvey@gmail.com

# Pull together stuff to make a nice odonate species level dataset on mite
# community characteristics.

library(tidyverse)
library(lubridate)
library(ape)
library(phytools)

odonates <- read_csv("datasets_derived/odonates_merged.csv")
mites <- read_csv("datasets_derived/mite_summaries.csv")

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
  select(c(species, simpson, num_host)) %>%
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
  mutate(generalist_nhost_10 = ifelse(num_host >  10, 1, 0)) %>%
  group_by(species) %>%
  summarise_all(sum) %>%
  select(c(species, specialist, generalist, specialist_25, generalist_25, 
           specialist_50, generalist_50, specialist_nhost, generalist_nhost,
           specialist_nhost_5, generalist_nhost_5, specialist_nhost_10, 
           generalist_nhost_10))

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

write_csv(odonates, "datasets_derived/odonate_summaries.csv")
