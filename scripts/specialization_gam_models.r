# Author: Jacob Harvey - jakekharvey@gmail.com  

library(MRFtools)
library(ape)
library(tidyverse)
library(mgcv)
library(car)
library(gratia)

odonates <- read_csv("datasets_derived/odonate_summaries_asv.csv")
                                              # _otu90 _otu97 asv
odonates <- odonates %>% 
  drop_na(specialist) %>%
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

# Sensitivity analysis for phylogenetic v non-phylogenetic indices.
# Still significnant, but larger p and lower R sq.
spec_gen_nhost <- cbind(odonates$specialist_nhost, odonates$generalist_nhost)
phylo_model_nhost <- gam(spec_gen ~ s(species, bs = 'mrf', 
                                xt = list(penalty = pmatrix), k = 10), 
                   data = odonates, method = "REML", family = binomial)
# summary(phylo_model)
# appraise(phylo_model)

#-------------------------------------------------------------------------------
# Does mass predict the proportion of specialists on an odonate?
#

mass_model <- gam(spec_gen ~ log(mass), 
         data = odonates, method = "REML", family = binomial)
#summary(mass_model)
#appraise(mass_model)

# Sensitivity analysis
# Becomes marginally significant. V low R sq.
mass_model_nhost <- gam(spec_gen_nhost ~ log(mass), 
                  data = odonates, method = "REML", family = binomial)
# summary(mass_model_nhost)
# appraise(mass_model_nhost)

#-------------------------------------------------------------------------------
# Is specialization predicted by host abundance?
#

odonates_abun <- drop_na(odonates, abundance)
spec_gen_abun <- cbind(odonates_abun$specialist, odonates_abun$generalist)

abun_model <- gam(spec_gen_abun ~ log(abundance), 
                  data = odonates_abun, method = "REML", family = binomial)
#summary(abun_model)
#appraise(abun_model)

# Sensitivity analysis
# Still significant, higher R sq ????? (slightly tho, basically the same)
spec_gen_abun_nhost <- cbind(odonates_abun$specialist_nhost, odonates_abun$generalist_nhost)
abun_model_nhost <- gam(spec_gen_abun_nhost ~ log(abundance), 
                        data = odonates_abun, method = "REML", family = binomial)
# summary(abun_model_nhost)
# appraise(abun_model_nhost)
