# Author: Jacob Harvey - jakekharvey@gmail.com  
#
# Models looking at the distribution of mite specialization on odonate hosts.

library(MRFtools)
library(ape)
library(tidyverse)
library(mgcv)
library(gratia)
library(phytools)

odonates <- read_csv("datasets_derived/odonate_summaries_predicted.csv") %>%
  drop_na(mpd_mean)
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

odonates <- mutate(odonates, median_p5_diff = ses_median - ses_5)

phylo_model <- gam(median_p5_diff ~ s(species, bs = 'mrf', 
                                  xt = list(penalty = pmatrix), k = 10),
                   data = odonates, method = "REML", family = gaussian)
#summary(phylo_model)
#appraise(phylo_model)

#-------------------------------------------------------------------------------
# Is specialization predicted by host immune response?
#

immune_model <- gam(median_p5_diff ~ log(immune_response), 
                  data = odonates, method = "REML", family = gaussian)
#summary(immune_model)
#appraise(immune_model)

#-------------------------------------------------------------------------------
# Is specialization predicted by host relative abundance?
#

abundance_model <- gam(median_p5_diff ~ log(abundance), 
                    data = odonates, method = "REML", family = gaussian)
#summary(abundance_model)
#appraise(abundance_model)

#-------------------------------------------------------------------------------
# Does this pattern hold at different phylogenetic scales? Doubles as sensitivity
# tests for sampling biases between different taxa.
#

# Do it only on Zygoptera.
zygoptera <- odonates %>% 
  filter(suborder == "Zygoptera")

zygoptera_mrca <- getMRCA(phylo, c("Lestes_congener", "Enallagma_hageni"))
z_phylo <- getDescendants(phylo, zygoptera_mrca)
z_phylo <- keep.tip(phylo, z_phylo)
z_phylo <- keep.tip(z_phylo, as.character(zygoptera$species))

z_pmatrix <- mrf_penalty(z_phylo)

zygoptera_phylo_model <- gam(median_p5_diff ~ s(species, bs = 'mrf', 
                                          xt = list(penalty = z_pmatrix), k = 10), 
                             data = zygoptera, method = "REML", family = gaussian)
#summary(zygoptera_phylo_model)
#appraise(zygoptera_phylo_model)

zygoptera_immune_model <- gam(median_p5_diff ~ log(immune_response),
                             data = zygoptera, method = "REML", family = gaussian)
#summary(zygoptera_immune_model)
#appraise(zygoptera_phylo_model)


# Do it only on Anisoptera.
#------------------------------------------------------------------------------
anisoptera <- odonates %>% 
  filter(suborder == "Anisoptera")

anisoptera_mrca <- getMRCA(phylo, c("Gomphus_spicatus", "Leucorrhinia_intacta"))
a_phylo <- getDescendants(phylo, anisoptera_mrca)
a_phylo <- keep.tip(phylo, a_phylo)
a_phylo <- keep.tip(a_phylo, as.character(anisoptera$species))

a_pmatrix <- mrf_penalty(a_phylo)

anisoptera_phylo_model <- gam(median_p5_diff ~ s(species, bs = 'mrf', 
                                            xt = list(penalty = a_pmatrix), k = 3), 
                             data = anisoptera, method = "REML", family = gaussian)
#summary(anisoptera_phylo_model)
#appraise(anisoptera_phylo_model) 

# Nope, But we've got a tiny sample of Anisoptera alone (n = 8)
anisoptera_immune_model <- gam(median_p5_diff ~ log(immune_response), data = anisoptera,
                               method = "REML", family = gaussian)
#summary(anisoptera_immune_model)
#appraise(anisoptera_immune_model)
