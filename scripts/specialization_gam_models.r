# Author: Jacob Harvey - jakekharvey@gmail.com  

library(MRFtools)
library(ape)
library(tidyverse)
library(mgcv)
library(car)
library(gratia)
library(phytools)

#odonates <- read_csv("datasets_derived/odonate_summaries_asv.csv")
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

#-------------------------------------------------------------------------------
# Is specialization predicted by host immune response?
#

odonates_immune <- drop_na(odonates, immune_response)
spec_gen_immune <- cbind(odonates_immune$specialist, odonates_immune$generalist)

immune_model <- gam(spec_gen_immune ~ log(immune_response), 
                  data = odonates_immune, method = "REML", family = binomial)
# summary(immune_model)
# appraise(immune_model)

# Sensitivity analysis
# Marginal p, half the deviance explained.
# This model is real heteroscedastic tho.
spec_gen_immune_nhost <- cbind(odonates_immune$specialist_nhost, 
                               odonates_immune$generalist_nhost)
immune_model_nhost <- gam(spec_gen_immune_nhost ~ log(immune_response), 
                          data = odonates_immune, method = "REML", family = binomial)
# summary(immune_model_nhost)
# appraise(immune_model_nhost)

#-------------------------------------------------------------------------------
# Does this pattern hold at different phylogenetic scales? Doubles as sensitivity
# tests for sampling biases between different taxa.
#

# Do it only on Zygoptera.
zygoptera <- odonates %>% 
  filter(suborder == "Zygoptera") %>%
  drop_na(immune_response) %>%
  drop_na(abundance)

zygoptera_mrca <- getMRCA(phylo, c("Lestes_congener", "Enallagma_hageni"))
z_phylo <- getDescendants(phylo, zygoptera_mrca)
z_phylo <- keep.tip(phylo, z_phylo)
z_phylo <- keep.tip(z_phylo, as.character(zygoptera$species))

z_pmatrix <- mrf_penalty(z_phylo)
z_spec_gen <- cbind(zygoptera$specialist, zygoptera$generalist)

zygoptera_phylo_model <- gam(z_spec_gen ~ s(species, bs = 'mrf', 
                                          xt = list(penalty = z_pmatrix), k = 10), 
                             data = zygoptera, method = "REML", family = binomial)
# Still a 'lil significant, high deviance explained too.
summary(zygoptera_phylo_model)
#appraise(zygoptera_phylo_model)

zygoptera_immune_model <- gam(z_spec_gen ~ log(immune_response),
                             data = zygoptera, method = "REML", family = binomial)
# Still a 'lil significant, high deviance explained too.
summary(zygoptera_immune_model)
#appraise(zygoptera_phylo_model)

zygoptera_abun_model <- gam(z_spec_gen ~ log(abundance), data = zygoptera, 
                            method = "REML", family = binomial)
summary(zygoptera_abun_model)
#appraise(zygoptera_abun_model)

# Do it only on Anisoptera.
anisoptera <- odonates %>% 
  filter(suborder == "Anisoptera") %>%
  drop_na(immune_response) %>%
  drop_na(abundance)

anisoptera_mrca <- getMRCA(phylo, c("Gomphus_spicatus", "Sympetrum_costiferum"))
a_phylo <- getDescendants(phylo, anisoptera_mrca)
a_phylo <- keep.tip(phylo, a_phylo)
a_phylo <- keep.tip(a_phylo, as.character(anisoptera$species))

a_pmatrix <- mrf_penalty(a_phylo)
a_spec_gen <- cbind(anisoptera$specialist, anisoptera$generalist)

anisoptera_phylo_model <- gam(a_spec_gen ~ s(species, bs = 'mrf', 
                                            xt = list(penalty = a_pmatrix), k = 3), 
                             data = anisoptera, method = "REML", family = binomial)
# Sample size just really seems too low here.
summary(anisoptera_phylo_model)
#appraise(anisoptera_phylo_model) 

anisoptera_immune_model <- gam(a_spec_gen ~ log(immune_response), data = anisoptera,
                               method = "REML", family = binomial)
summary(anisoptera_immune_model)

anisoptera_abun_model <- gam(a_spec_gen ~ log(abundance), data = anisoptera,
                               method = "REML", family = binomial)
summary(anisoptera_abun_model)

#
# Sample size of 8 is just not high enough to do anything with Anisoptera. Besides we
# only really have interactions for the smaller dragonflies, We'd need some Aeshnidae
# if we wanted dragonflies which actually contrast against eachother. 
#
