# Author: Jacob Harvey - jakekharvey@gmail.com
# Summaries from the mite pov.

library(tidyverse)

network <- read_csv("datasets_derived/bin_network.csv")
mites <- read_csv("datasets_derived/mite_phylo_scale.csv")
odonates <- read_csv("datasets_derived/odonate_summaries.csv")

# Get mite statistics at the host-species level.
get_scale <- Vectorize(function(mite){mites$combo_phylo_scale[which(mites$mite == mite)]})
get_rr <- Vectorize(function(mite){mites$resource_range[which(mites$mite == mite)]})
get_hn <- Vectorize(function(mite){mites$num_host[which(mites$mite == mite)]})

mite_points <- network %>% 
  pivot_longer(2:last_col()) %>%
  filter(value > 0) %>%
  rename(mite = name) %>%
  select(odonate_spp, mite) %>%
  mutate(mite_scale = get_scale(mite)) %>%
  mutate(resource_range = get_rr(mite)) %>%
  mutate(num_host = get_hn(mite)) %>%
  select(odonate_spp, mite_scale, resource_range, num_host) %>%
  rename(species = odonate_spp)

odonates <- select(odonates, species, prevalence, mass)

mites <- inner_join(mite_points, odonates, by = "species")

write_csv(mites, "datasets_derived/mite_summaries.csv")

