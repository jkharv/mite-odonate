# Author: Jake Harvey - jakekharvey@gmail.com

library(tidyverse)
library(mgcv)
library(MRFtools)
library(ape)

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

pmatrix <- mrf_penalty(phylo)

# Mass is highly phylogenetically conserved. ~ 80% deviance explained. 
mass_phylo_model <- gam(mass ~ s(species, bs = 'mrf', xt = 
                                 list(penalty = pmatrix), k = 65), 
                        data = masses, method = "REML", family = gaussian)
summary(phylo_model)

mites <- read_csv("datasets_derived/mite_summaries.csv") %>%
  select(mite, faiths, simpson, num_host, mite_scale) %>%
  distinct()

n_model <- gam(mite_scale ~ num_host, data = mites, method = "REML", 
               family = gaussian)
summary(n_model)

pdps_model <- gam(mite_scale ~ s(faiths), data = mites, method = "REML",
                  family = gaussian)
summary(pdps_model)

summary(lm(mite_scale ~ simpson, data = mites))


pd_plot <- ggplot(mites, aes(x = faiths, y = mite_scale)) +
  geom_point() +
  scale_x_log10()
print(pd_plot)





