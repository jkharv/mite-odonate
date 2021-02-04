# Author: Jake Harvey jakekharvey@gmail.com

library(tidyverse)
library(ggstance)
library(ggtree)
library(ape)

mites <- read_csv("datasets_derived/mite_phylo_scale.csv")
network <-read_csv("datasets_derived/bin_network.csv")
phylo <- read.nexus("datasets_primary/phylogeny/FinalOdonatetree.tre")

# Filtering out species that we don't have in our phylogeny.
network <- network %>%
  filter(odonate_spp != "Enallagma_annexum") %>%
  filter(odonate_spp != "Lestes_forcipatus") %>%
  filter(odonate_spp != "Lestes_inaequalis") %>%
  filter(odonate_spp != "Lestes_vigilax")

# Getting the data into a plotable format. One point for each mite on each Odonate.
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
  mutate(num_host = get_hn(mite))

trimmed_phylo <- keep.tip(phylo, network$odonate_spp)

tree <- ggtree(trimmed_phylo) +
  geom_tiplab() + 
  xlim_tree(450)

scale_plot <- facet_plot(tree, panel = "Phylogenetic Scale", data = mite_points,
                         geom = geom_point, aes(x = mite_scale)) +
                         xlim_expand(c(25, 100), panel = "Phylogenetic Scale") + 
                         theme_tree2()
print(scale_plot)

rr_plot <- facet_plot(tree, panel = "Phylogenetic Scale", data = mite_points,
                      geom = geom_point, aes(x = resource_range)) +
                      xlim_expand(c(0, 1), panel = "Phylogenetic Scale") + 
                      theme_tree2()
print(rr_plot)

hr_plot <- facet_plot(tree, panel = "Phylogenetic Scale", data = mite_points,
                      geom = geom_point, aes(x = num_host)) +
                      xlim_expand(c(0, 1), panel = "Phylogenetic Scale") + 
                      theme_tree2()
print(hr_plot)




facet_widths(facet_plot, c(3,1))