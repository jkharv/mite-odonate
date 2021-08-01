# Author: Jake Harvey -- jakekharvey@gmail.com

library(tidyverse)
library(mgcv)
library(gratia)
library(MRFtools)
library(ape)

mite_phylo <- read.tree("datasets_derived/mite_tree.tre")
odonate_phylo <- read.nexus("datasets_primary/phylogeny/FinalOdonatetree.tre")

edges <- read_csv("datasets_derived/edge_list.csv")
odonates <- read_csv("datasets_derived/odonate_summaries.csv")

# Filter out odonate species which we don't have immune response data for.
# Select only the cols we need for this analysis. Also remove odonates we have 
# no parasitism records for.
odonates <- odonates %>%
  drop_na(immune_response, parasite_interactions, abundance) %>%
  select(c(species, abundance, immune_response))

# Keep only mites which are both in the dataset and the phylogeny.
keep <- intersect(edges$mite, mite_phylo$tip.label)
edges <- filter(edges, mite %in% keep)
mite_phylo <- keep.tip(mite_phylo, keep)

# Keep only odonates which are both in the dataset and the phylogeny. 
keep <- intersect(odonates$species, odonate_phylo$tip.label)
odonates <- filter(odonates, species %in% keep)
odonate_phylo <- keep.tip(odonate_phylo, keep)

# Not all the mite-level data is relevant for this analysis. The phylo-based
# specialization measured would be double counting the phylogeny.
# Also this still references some of the odonates we just removed.
edges <- edges %>%
  select(odonate_spp, mite, interaction) %>%
  filter(odonate_spp %in% keep) %>%
  rename(species = odonate_spp)

# Join the odonate data to the edge-list, so we can use both mite-level and 
# odonate-level predictors in the model.
edges <- left_join(edges, odonates, by = "species")

# Mites and odonate names need to be factors
edges$species <- factor(edges$species)
edges$mite <- factor(edges$mite)

# Set up penalty matrices for the mite and odonate phylogenies.
m_pmatrix <- mrf_penalty(mite_phylo)
o_pmatrix <- mrf_penalty(odonate_phylo)

co_phylo_model <- gam(interaction ~ s(species, bs = 'mrf', 
                                      xt = list(penalty = o_pmatrix), k = 10) +
                                    s(mite, bs = 'mrf', 
                                      xt = list(penalty = m_pmatrix), k = 10) +
                                    ti(species, mite, bs = c('mrf', 'mrf'), 
                                       k = c(10, 10),
                                       xt = list(list(penalty = o_pmatrix),
                                                 list(penalty = m_pmatrix))) +
                                    log(immune_response) + 
                                    log(abundance),
                   data = edges, method = "REML", family = binomial)
#summary(co_phylo_model)
#appraise(co_phylo_model)

immune_model <- gam(interaction ~ s(mite, bs = 'mrf', 
                                    xt = list(penalty = m_pmatrix), k = 10) +
                                  log(immune_response),
                    data = edges, method = "REML", family = binomial)
#summary(immune_model)
#appraise(immune_model)

immune_cons <- gam(log(immune_response) ~ s(species, bs = 'mrf', 
                                             xt = list(penalty = o_pmatrix), k = 10),
                    data = edges, method = "REML", family = gaussian)
#summary(immune_cons)
#appraise(immune_cons)

#-------------------------------------------------------------------------------
#                                                                              |
#                 Make some predictions based on the model                     |
#                                                                              |
#-------------------------------------------------------------------------------

edges <- edges %>%
  unite("int", c("species", "mite"), remove = FALSE) %>%
  distinct(int, .keep_all = TRUE) %>% # Only need to predict an interaction once
  select(species, mite, interaction, immune_response, abundance)

predicted <- predict(co_phylo_model, edges)
edges <- add_column(edges, predicted)

logit_to_prob <- function(logit){
  
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

edges <- mutate(edges, prob  = logit_to_prob(predicted))
edges <- select(edges, c(species, mite, prob))

adj_matrix <- pivot_wider(edges, names_from = mite, 
                          values_from = prob)
write_csv(adj_matrix, "datasets_derived/predicted_network.csv")

