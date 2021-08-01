# Author: Jake Harvey jakekharvey@gmail.com

library(tidyverse)
library(ggstance)
library(ggtree)
library(ape)
library(ggridges)

mites <- read_csv("datasets_derived/mite_summaries_predicted.csv")
phylo <- read.nexus("datasets_primary/phylogeny/FinalOdonatetree.tre")

# I have data for some odonates which aren't in the phylogeny.
# These are left out of all the analysis.
keep <- intersect(phylo$tip.label, mites$species)
mites <- filter(mites, species %in% keep)
phylo <- keep.tip(phylo, keep)

tree <- ggtree(phylo, ladderize = TRUE) +
  geom_tiplab() + 
  xlim_tree(500) +
  theme_tree2()

scale_plot <- facet_plot(tree, panel = "Phylogenetic Scale", data = mites,
                         geom_density_ridges, mapping = 
                           aes(group = label, x = ses_mpd * interaction))

ggsave("figures/pscale_phylo_plot.svg", scale_plot, width = 15, height = 8)
