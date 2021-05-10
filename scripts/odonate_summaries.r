# Author: Jacob Harvey - jakekharvey@gmail.com

library(tidyverse)
library(ape)
library(phytools)

a2015 <- read_csv("datasets_primary/2015_data.csv", na = c("", "NA", "N/A", "n/a"))
a2019 <- read_csv("datasets_primary/2019_data.csv", na = c("", "NA", "N/A", "n/a"))
a2020 <- read_csv("datasets_primary/2020_data.csv", na = c("", "NA", "N/A", "n/a"))
masses <- read_csv("datasets_primary/2019-2020-mass.csv", na = c("", "NA", "N/A", "n/a"))
network <- read_csv("datasets_derived/bin_network.csv")
mites <- read_csv("datasets_derived/mite_phylo_scale.csv")

a2015 <- a2015 %>%
  mutate(species = paste(Genus, Species, sep = "_")) %>%
  rename(Mite_Abundance = "Mite Load") %>%
  rename(mass = "Mass (g)") %>%
  filter(species != "Lestes_sp" & species != "Aeshna_sp") %>%
  drop_na(Species) %>%
  select(c(species, Mite_Abundance, mass))

# Julie used old taxonomy without Gomphus split up. I'm using Gomphus in the analysis.
# Gomphurus, Hylogomphus, Phanogomphus, and Stenogomphurus -> Gomphus
a2020 <-  a2020 %>%
  mutate_if(is.character, str_replace_all, 
            pattern = "Gomphurus|Hylogomphus|Phanogomphus|Stenogomphurus", 
            replacement = "Gomphus") %>%
  mutate(species = paste(Genus, Species, sep = "_")) %>%
  drop_na(Species) %>%
  select(species, Mite_Abundance)

a2019 <-  a2019 %>%
  mutate_if(is.character, str_replace_all, 
            pattern = "Gomphurus|Hylogomphus|Phanogomphus|Stenogomphurus", 
            replacement = "Gomphus") %>%
  mutate(species = paste(Genus, Species, sep = "_")) %>%
  drop_na(Species) %>%
  select(species, Mite_Abundance)

masses <-  masses %>%
  mutate_if(is.character, str_replace_all, 
            pattern = "Gomphurus|Hylogomphus|Phanogomphus|Stenogomphurus", 
            replacement = "Gomphus") %>%
  mutate(species = paste(Genus, Species, sep = "_")) %>%
  rename(mass = "Weight(g)") %>%
  drop_na(Species, mass) %>%
  select(species, mass)

# Julie's data used a standard sampling protocol.
abundances <- a2015 %>%
  group_by(species) %>%
  tally() %>%
  rename(abundance = n)

all_data <- union_all(union_all(a2015, a2019), union_all(a2020, masses))

prevalences <- all_data %>%
  drop_na(Mite_Abundance) %>%
  mutate(abundance = 1) %>%
  group_by(species) %>%
  summarise_all(sum) %>%
  mutate(prevalence = Mite_Abundance/abundance) %>%
  select(c(species, prevalence))

masses <- all_data %>%
  drop_na(mass) %>%
  mutate(abundance = 1) %>%
  group_by(species) %>%
  summarise_all(sum) %>%
  mutate(mass = mass/abundance) %>%
  select(c(species, mass))
  
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
  select(odonate_spp, mite_scale, resource_range, num_host)

specialist <- mite_points %>%
  mutate(specialist = ifelse(mite_scale <= 100, 1, 0)) %>%
  mutate(generalist = ifelse(mite_scale > 100, 1, 0)) %>%
  mutate(specialist_50 = ifelse(mite_scale <= 50, 1, 0)) %>%
  mutate(generalist_50 = ifelse(mite_scale > 50, 1, 0)) %>%
  mutate(specialist_25 = ifelse(mite_scale <= 25, 1, 0)) %>%
  mutate(generalist_25 = ifelse(mite_scale > 25, 1, 0)) %>%
  select(odonate_spp, specialist, generalist, specialist_25, 
         generalist_25, specialist_50, generalist_50) %>%
  group_by(odonate_spp) %>%
  summarise_all(sum)  %>%
  rename(species = odonate_spp) %>%
  mutate(spec_ratio = specialist/generalist)

spp_avg <- mite_points %>%
  group_by(odonate_spp) %>%
  summarise_all(c(mean = mean, sd = sd)) %>%
  rename(species = odonate_spp)

odonates <- list(prevalences, masses, abundances, specialist, spp_avg) %>%
            reduce(function(x, y) full_join(x, y, by = "species")) 
                     
# Add col for suborder, useful for plots later
phylo <- read.nexus("datasets_primary/phylogeny/FinalOdonatetree.tre")
zygoptera_mrca <- getMRCA(phylo, c("Lestes_congener", "Enallagma_hageni"))
zygoptera <- getDescendants(phylo, zygoptera_mrca)
zygoptera <- keep.tip(phylo, zygoptera)

odonates <- mutate(odonates, suborder = if_else(species %in% zygoptera$tip.label, 
                                                "Zygoptera", "Anisoptera"))

write_csv(odonates, "datasets_derived/odonate_summaries.csv")
