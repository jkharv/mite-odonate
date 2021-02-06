# Author: Jacob Harvey - jakekharvey@gmail.com  

library(MRFtools)
library(ape)
library(tidyverse)
library(mgcv)

odonates <- read.csv("datasets_derived/odonate_summaries.csv")

odonates <- odonates %>%
  filter(mite_scale > 0)

odonates$species = factor(odonates$species)

phylo <- read.nexus("datasets_primary/phylogeny/FinalOdonatetree.tre")
# Phylogeny is missing many species. I have mite data for four of the missing.
keep <- intersect(phylo$tip.label, odonates$species)
odonates <- filter(odonates, species %in% keep)
phylo <- keep.tip(phylo, keep)
plot(phylo)

pmatrix <- mrf_penalty(phylo)

x <- gam(mite_scale ~ s(species, bs = 'mrf', xt = list(penalty = pmatrix)), data = odonates,
         method = "REML", family = gaussian)

summary(x)
