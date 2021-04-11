# Author: Jacob Harvey - jakekharvey@gmail.com  

library(MRFtools)
library(ape)
library(tidyverse)
library(mgcv)
library(car)
library(gratia)

odonates <- read_csv("datasets_derived/odonate_summaries.csv")

odonates <- odonates %>% 
  drop_na(mite_scale_mean) %>%
  filter(generalist > 0) # no div zero
odonates$species <- factor(odonates$species)

phylo <- read.nexus("datasets_primary/phylogeny/FinalOdonatetree.tre")
# Phylogeny is missing many species. I have mite data for four of the missing.
keep <- intersect(phylo$tip.label, odonates$species)
odonates <- filter(odonates, species %in% keep)
phylo <- keep.tip(phylo, keep)

#-------------------------------------------------------------------------------
# Is there phylogenetic signal in ratio of specialists on a host?
#

pmatrix <- mrf_penalty(phylo)
spec_gen <- cbind(odonates$specialist, odonates$generalist)

phylo_model <- gam(spec_gen ~ s(species, bs = 'mrf', 
                                  xt = list(penalty = pmatrix), k = 10), 
                   data = odonates, method = "REML", family = binomial)
#summary(phylo_model)

#appraise(phylo_model)

#-------------------------------------------------------------------------------
# Does mass predict the proportion of specialists on an odonate?
#

mass_model <- gam(spec_gen ~ log(mass), 
         data = odonates, method = "REML", family = binomial)
#summary(mass_model)

#appraise(mass_model)

#-------------------------------------------------------------------------------
# Is specialization predicted by host abundance?
#

odonates <- drop_na(odonates, abundance)
spec_gen <- cbind(odonates$specialist, odonates$generalist)

abun_model <- gam(spec_gen ~ abundance, 
                  data = odonates, method = "REML", family = binomial)
#summary(abun_model)

#appraise(abun_model)
