# Author: Jake Harvey - jakekharvey@gmail.com

library(tidyverse)
library(mgcv)
library(MRFtools)
library(ape)
library(lmtest)

a201920 <- read_csv("datasets_primary/2019-2020-mass.csv", na = c("n/a", "", "NA"))
a2015 <- read_csv("datasets_primary/2015_data.csv", na = c("n/a", "", "NA"))

a201920 <- a201920 %>%
  rename(mass = "Weight(g)") %>%
  unite(species, c("Genus", "Species")) %>%
  select(c(mass, species)) %>%
  drop_na(c(mass, species))

a2015 <- a2015 %>%
  rename(mass = "Mass (g)") %>%
  unite(species, c("Genus", "Species")) %>%
  select(c(mass, species)) %>%
  drop_na(c(mass, species))

masses <- rbind(a2015, a201920)
masses$species <- factor(masses$species)

phylo <- read.nexus("datasets_primary/phylogeny/FinalOdonatetree.tre")
# Phylogeny is missing many species. I have mite data for four of the missing.
keep <- intersect(phylo$tip.label, masses$species)
masses <- filter(masses, species %in% keep)
phylo <- keep.tip(phylo, keep)

#-------------------------------------------------------------------------------
# Is mass phylogenetically conserved?
#

pmatrix <- mrf_penalty(phylo)

# Mass is highly phylogenetically conserved. ~ 80% deviance explained. 
mass_phylo_model <- gam(log(mass) ~ s(species, bs = 'mrf', xt = 
                                 list(penalty = pmatrix), k = 65), 
                        data = masses, method = "REML", family = gaussian)
#summary(mass_phylo_model)
#gam.check(mass_phylo_model)

#-------------------------------------------------------------------------------
# Is abundance phylogenetically conserved?
#

odonates <- read_csv("datasets_derived/odonate_summaries.csv") %>%
  drop_na(abundance)
odonates$species <- factor(odonates$species)

keep <- intersect(phylo$tip.label, odonates$species)
odonates <- filter(odonates, species %in% keep)
abun_phylo <- keep.tip(phylo, keep)
abun_pmatrix <- mrf_penalty(abun_phylo)

abun_phylo_model <- gam(log(abundance) ~ s(species, bs = 'mrf', xt = 
                                      list(penalty = abun_pmatrix)), 
                        data = odonates, method = "REML", family = gaussian)
#summary(abun_phylo_model)
#appraise(abun_phylo_model)

#-------------------------------------------------------------------------------
# Is mite scale biased by number of hosts (sampling effort)?
#

mites <- read_csv("datasets_derived/mite_summaries.csv") %>%
  select(mite, faiths, simpson, num_host, mite_scale) %>%
  distinct()

n_model <- gam(mite_scale ~ log(num_host), data = mites, method = "REML", 
               family = gaussian)
#summary(n_model)
#appraise(n_model)
