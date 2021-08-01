# Author: Jake Harvey -- jakekharvey@gmail.com

# This script preps the data for the co-phylogenetic models we will run to 
# predict odonate - mite interactions.

library(tidyverse)
library(lubridate)
library(ape)
library(phytools)
library(mgcv)

# For some reason, readr has trouble finding the column type when it begins in a
# bunch of NAs. It guesses col_logical and then trips up when the values are not
# logical.
odonates <- read_csv("datasets_derived/odonates_merged.csv",
                     col_types = cols(melanization = col_double(),
                                      control = col_double()))
odonates <- odonates %>%
  unite(genus_species, c(genus, species), sep = "_", remove = FALSE)

# We just need a field to tell if a species has parasite data or no for running
# our models later on. We need to filter out species without parasite data.
network <- read_csv("datasets_derived/bin_network.csv")

network <- network %>%
  pivot_longer(2:last_col()) %>%
  select(odonate_spp, value) %>%
  group_by(odonate_spp) %>%
  summarise(parasite_interactions = sum(value)) %>%
  rename(species = odonate_spp)

# The 2015 data was collected using a standard sampling protocol.
# Estimate site-level relative abundance for each species at each site then
# average each estimate from every site at the species level.
abundances <- odonates %>%
  filter(year(date) == 2015) %>%
  select(genus_species, site) %>%
  group_by(site, genus_species) %>%
  tally() %>%
  group_by(site) %>%
  mutate(abundance = n/sum(n)) %>%
  group_by(genus_species) %>%
  select(genus_species, abundance) %>%
  summarize(abundance = mean(abundance)) %>%
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

# The lighting setup was different in 2019 and 2020. We'll need to correct the
# melanization values for this real quick using a random intercept model and
# extracting the residuals to use as immune response values.
od <- odonates %>%
  drop_na(melanization)

# From Maggie's code
yr <- as.factor(year(od$date))

# Control is the darkness of the unmelanized side of the implant.
immune_response <- od$melanization - od$control
ir <- lme4::lmer(immune_response ~ (1|yr), data = od)

immune_responses <- od %>%
  add_column(immune_response = immune_response) %>%
  select(c(genus_species, immune_response)) %>%
  group_by(genus_species) %>%
  summarise(across(everything(), mean)) %>%
  rename(species = genus_species)

# Join these dataframes together
odonates <- list(abundances, masses, immune_responses, network) %>%
  reduce(function(x, y) full_join(x, y, by = "species")) 

# Add a column for suborder which will be useful later for plots.
phylo <- read.nexus("datasets_primary/phylogeny/FinalOdonatetree.tre")
zygoptera_mrca <- getMRCA(phylo, c("Lestes_congener", "Enallagma_hageni"))
zygoptera <- getDescendants(phylo, zygoptera_mrca)
zygoptera <- keep.tip(phylo, zygoptera)

odonates <- mutate(odonates, suborder = if_else(species %in% zygoptera$tip.label, 
                                                "Zygoptera", "Anisoptera"))

write_csv(odonates, "datasets_derived/odonate_summaries.csv")

